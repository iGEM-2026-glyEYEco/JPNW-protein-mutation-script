import os
import re
import shutil

from Bio import SeqIO
from itertools import product

def generate_mutations(input_file_name, input_positions, output_file_name, max_mutations=None, position_mutation_options=None, chain="A", pdb_offset=0):
    """
    Parameters (inputs):
    - input_file_name: The name of the FASTA file containing the scaffold protein sequence.
    - input_positions: A list of positions (each indexed from 1) in the protein sequence where the residue will be mutated
    - output_file_name: The name of the output individual_list.txt file for FoldX.
    - max_mutations (optional): The maximum number of positions allowed to be mutated. If None, all possible mutations will be generated.
    - position_mutation_options (optional): A dictionary mapping positions to lists of allowed amino acid options.
    - chain (optional): PDB chain identifier (default "A").
    - pdb_offset (optional): Added to sequence position to get PDB residue number (default 0).

    Outputs:
    - A FoldX individual_list.txt file where each line encodes one variant.
      Format: <WT_residue><chain><PDB_number><mutant_residue>,...;
      e.g. single: RA144G;   combinatorial: RA144G,KB78R;
    """
    positions = [position - 1 for position in input_positions]  # Convert to 0-based indexing
    amino_acids = list("ACDEFGHIKLMNPQRSTVWY")

    # Read the scaffold protein sequence
    scaffold_record = SeqIO.read(input_file_name, "fasta")
    scaffold_sequence = str(scaffold_record.seq)

    # Generate mutation options per position
    mutation_options = []
    for pos_1based, pos_0based in zip(input_positions, positions):
        if scaffold_sequence[pos_0based] == "-":  # Account for gaps in the scaffold sequence (no mutations at gaps)
            mutation_options.append(["-"])
        elif position_mutation_options is not None and pos_1based in position_mutation_options:
            mutation_options.append(position_mutation_options[pos_1based])
        else:
            mutation_options.append(amino_acids)

    # List of all possible combinations of mutations at the specified positions
    all_combinations = product(*mutation_options)

    lines_written = 0
    with open(output_file_name, "w") as out:
        for combo in all_combinations:
            mutations = []
            mutation_count = 0

            for pos_1based, pos_0based, amino_acid in zip(input_positions, positions, combo):
                wt_aa = scaffold_sequence[pos_0based]
                if wt_aa == "-":
                    continue  # Skip gaps

                if amino_acid != wt_aa:
                    mutation_count += 1
                    pdb_number = pos_1based + pdb_offset
                    mutations.append(f"{wt_aa}{chain}{pdb_number}{amino_acid}")

            # Skip this combination if it exceeds the maximum allowed mutations
            if max_mutations is not None and mutation_count > max_mutations:
                continue

            # Skip identity (no mutations) — nothing to write to FoldX
            if not mutations:
                continue

            out.write(",".join(mutations) + ";\n")
            lines_written += 1

    print(f"individual_list.txt written to '{output_file_name}'")
    print(f"  Variants written: {lines_written}")
    return lines_written


def _parse_mutation_name(line, multi_chain):
    """
    Convert a raw individual_list.txt line into a human-readable name.

    e.g. "NA143A,VB78C;"  →  "N143A_V78C"  (single chain)
                          →  "NA143A_VB78C" (multi-chain)
    """
    line = line.strip().rstrip(";")
    tokens = line.split(",")
    parts = []
    for token in tokens:
        m = re.fullmatch(r"([A-Z])([A-Z])(\d+)([A-Z])", token.strip())
        if not m:
            parts.append(token.strip())
            continue
        wt_aa, chain, pdb_num, mut_aa = m.groups()
        if multi_chain:
            parts.append(f"{wt_aa}{chain}{pdb_num}{mut_aa}")
        else:
            parts.append(f"{wt_aa}{pdb_num}{mut_aa}")
    return "_".join(parts)


def rename_foldx_outputs(individual_list_file, foldx_output_dir, pdb_base_name):
    """
    Rename FoldX BuildModel output files to use human-readable mutation names.

    After running FoldX BuildModel, outputs are numbered sequentially
    (e.g. 1A3K_Repair_1.pdb, 1A3K_Repair_2.pdb ...). This function:
      1. Reads individual_list.txt to build a list of mutation names
      2. Renames each PDB output file to include the mutation name
         (1A3K_Repair_1.pdb  →  1A3K_Repair_N143A.pdb)
      3. Adds a MutationName column to every .fxout file so each row is
         labelled with the actual mutation instead of just a number.

    """
    # --- Read and parse individual_list.txt ---------------------------------
    with open(individual_list_file) as f:
        raw_lines = [l.strip() for l in f if l.strip()]

    # Detect multi-chain: more than one distinct chain letter used
    all_chains = set(re.findall(r"[A-Z](?=[A-Z]\d+[A-Z])", " ".join(raw_lines)))
    multi_chain = len(all_chains) > 1

    mutation_names = [_parse_mutation_name(line, multi_chain) for line in raw_lines]

    # --- Rename PDB files ---------------------------------------------------
    renamed_pdbs = 0
    for i, name in enumerate(mutation_names, start=1):
        old_path = os.path.join(foldx_output_dir, f"{pdb_base_name}_{i}.pdb")
        new_path = os.path.join(foldx_output_dir, f"{pdb_base_name}_{name}.pdb")
        if os.path.exists(old_path):
            shutil.move(old_path, new_path)
            renamed_pdbs += 1

    print(f"PDB files renamed: {renamed_pdbs}")

    # --- Annotate .fxout files ----------------------------------------------
    fxout_pattern = re.compile(r"^(.*BuildModel.*\.fxout)$", re.IGNORECASE)
    annotated_fxout = 0

    for fname in os.listdir(foldx_output_dir):
        if not fxout_pattern.match(fname):
            continue

        fxout_path = os.path.join(foldx_output_dir, fname)
        with open(fxout_path) as f:
            lines = f.readlines()

        new_lines = []
        data_row_index = 0  # counts non-header, non-blank data rows
        for line in lines:
            stripped = line.rstrip("\n")

            # Header lines start with "Pdb" or are blank — pass through unchanged
            if not stripped or stripped.startswith("Pdb"):
                new_lines.append(line)
                continue

            # Data row: insert mutation name as a new first column
            if data_row_index < len(mutation_names):
                mut_label = mutation_names[data_row_index]
            else:
                mut_label = f"mutant_{data_row_index + 1}"

            new_lines.append(f"{mut_label}\t{stripped}\n")
            data_row_index += 1

        # Also update the header to include the new column
        for j, line in enumerate(new_lines):
            if line.startswith("Pdb"):
                new_lines[j] = f"MutationName\t{line}"
                break

        with open(fxout_path, "w") as f:
            f.writelines(new_lines)

        annotated_fxout += 1

    print(f".fxout files annotated: {annotated_fxout}")


if __name__ == "__main__":
    # ----- CHANGE INPUTS BELOW IF NEEDED -----
    input_file_name = "./test_inputs/scaffold_0.fasta"
    output_file_name = "../foldxMac/individual_list.txt"

    positions = [144]  # List of positions to mutate (1-based indexing)
    chain = "A"            # PDB chain identifier
    pdb_offset = 0         # Add to sequence position to get PDB residue number
    max_mutations = None   # Optional: Maximum number of positions allowed to be mutated (set to None to allow all mutations)

    # Optional: specify amino acid options for specific positions (set to None to use all 20 amino acids for all positions)
    # For any position not specified in position_mutation_options, the default amino acid options (all 20 amino acids) will be used.
    position_mutation_options = {
        144: list("ACDEFGHIKLMNPQRSTVWY"),  # Position 1 can mutate to any amino acid
    }
        
    
    # -----------------------------------------

    generate_mutations(input_file_name, positions, output_file_name, max_mutations, position_mutation_options, chain, pdb_offset)

    # ----- OPTIONAL: run AFTER FoldX BuildModel has finished -----
    # Uncomment and fill in the values below to rename PDB outputs and
    # annotate .fxout files with human-readable mutation names.
    #
    rename_foldx_outputs(
         individual_list_file = output_file_name,   # same file written above
         foldx_output_dir     = "../foldxMac",       # directory where FoldX wrote its outputs
         pdb_base_name        = "1A3K_Repair",       # base name FoldX used, e.g. "1A3K_Repair"
    )