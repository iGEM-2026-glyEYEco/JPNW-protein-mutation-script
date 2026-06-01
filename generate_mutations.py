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
        2. Renames each PDB output file to include the mutation name (1A3K_Repair_1.pdb  →  1A3K_Repair_N143A.pdb)
        3. Adds a MutationName column to every .fxout file so each row is labelled with the actual mutation instead of just a number.
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


def split_individual_list(input_file, output_dir,batch_size=100):
    """
    Splits a FoldX individual_list file into smaller batch files.

    Parameters:
    - input_file: Path to the FoldX individual_list file
    - output_dir: Directory where the batch files will be saved
    - batch_size: Number of lines (mutants) per batch file
    """
    os.makedirs(output_dir, exist_ok=True)

    with open(input_file) as f:
        lines = [line for line in f if line.strip()]

    # Get base filename without extension
    base_name = os.path.splitext(os.path.basename(input_file))[0]

    batch_files = []

    for i in range(0, len(lines), batch_size):
        batch_num = i // batch_size + 1

        batch_filename = os.path.join(
            output_dir,
            f"{base_name}_batch_{batch_num:03d}.txt"
        )

        with open(batch_filename, "w") as out:
            out.writelines(lines[i:i + batch_size])

        batch_files.append(batch_filename)

    print(f"Created {len(batch_files)} batch files")
    for bf in batch_files:
        print(f"  {bf}")
    
    return batch_files


def merge_pdbs(batch_dirs, merged_dir):
    """
    Merges PDB files from multiple batch directories into a single directory.
    Skips WT files and duplicates.

    Parameters:
    - batch_dirs: List of batch directories containing PDB outputs
    - merged_dir: Directory where merged PDB files will be saved
    """
    os.makedirs(merged_dir, exist_ok=True)

    copied = 0

    for batch_dir in batch_dirs:
        for fname in os.listdir(batch_dir):
            if fname.startswith("WT_"):
                continue

            if not fname.endswith(".pdb"):
                continue

            src = os.path.join(batch_dir, fname)
            dst = os.path.join(merged_dir, fname)

            if os.path.exists(dst):
                print(f"Duplicate skipped: {fname}")
                continue

            shutil.copy2(src, dst)
            copied += 1

    print(f"Copied {copied} mutant PDB files")


if __name__ == "__main__":
    # ------------------------------------------------------------------
    # 1. CHANGE INPUTS IF NEEDED
    # ------------------------------------------------------------------

    input_file_name = "./fasta_inputs/galectin3.fasta"
    output_file_name = "./foldx_inputs/individual_list_set5.txt"  # File name Must begin with individual_list

    positions = [52, 53]  # List of positions to mutate (1-based indexing)
    chain = "A"            # PDB chain identifier
    pdb_offset = 0         # Add to sequence position to get PDB residue number
    max_mutations = None   # Optional: Maximum number of positions allowed to be mutated (set to None to allow all mutations)

    # Optional: specify amino acid options for specific positions (set to None to use all 20 amino acids for all positions)
    # For any position not specified in position_mutation_options, the default amino acid options (all 20 amino acids) will be used.
    position_mutation_options = None   



    # ------------------------------------------------------------------
    # 2. GENERATE MUTATIONS (write individual_list file for FoldX)
    # Uncomment the line below to run.
    # Once finished, then comment it out again.
    # ------------------------------------------------------------------
    generate_mutations(input_file_name, positions, output_file_name, max_mutations, position_mutation_options, chain, pdb_offset)
    


    # ------------------------------------------------------------------
    # 3. Optional: If you have many mutants, split the individual_list.txt into smaller batches for FoldX. 
    # Uncomment the line below and modify any parameters if needed. 
    # Once finished, comment it out again.
    # ------------------------------------------------------------------
    # batch_files = split_individual_list(output_file_name, output_dir="./foldx_inputs/set5/", batch_size=100)



    # ------------------------------------------------------------------
    # 4. Run FoldX BuildModel in your terminal using the generated individual_list.txt file(s). 
    # This will produce PDB outputs and .fxout files in the specified output directory
    # ------------------------------------------------------------------



    # ------------------------------------------------------------------
    # 5. Optional: run AFTER FoldX BuildModel has finished
    # Uncomment and fill in the values below to rename PDB outputs and annotate .fxout files with human-readable mutation names.
    # Once finished, comment it out again.
    # ------------------------------------------------------------------

    # *** IF YOU DID NOT PERFORM STEP 4 (NO BATCHES)
    # rename_foldx_outputs(
    #     individual_list_file = "./foldx_inputs/set5/individual_list_set5_batch_002.txt",
    #     foldx_output_dir     = "./foldx_outputs/set5/batch_002",       # directory where FoldX wrote its outputs
    #     pdb_base_name        = "scaffold_galectin3_Repair",       # base name FoldX used, e.g. "1A3K_Repair"
    # )

    # *** IF YOU PERFORMED STEP 4 (IF YOU HAVE MULTIPLE BATCHES):
    # for batch_num in range(1, 5):
    #     rename_foldx_outputs(
    #         individual_list_file=f"./foldx_inputs/set5/individual_list_set5_batch_{batch_num:03d}.txt",
    #         foldx_output_dir=f"./foldx_outputs/set5/batch_{batch_num:03d}",
    #         pdb_base_name="scaffold_galectin3_Repair"
    # )



    # ------------------------------------------------------------------
    # 6. Optional: If you split into batches (performed step 4), uncomment the code below to merge the renamed PDB outputs from all batches into a single directory.
    # Change any parameters if needed.
    # One finished, comment it out again.
    # ------------------------------------------------------------------
    # merge_pdbs(
    #     batch_dirs=[
    #         "./foldx_outputs/set5/batch_001",
    #         "./foldx_outputs/set5/batch_002",
    #         "./foldx_outputs/set5/batch_003",
    #         "./foldx_outputs/set5/batch_004",
    #     ],
    #     merged_dir="./foldx_outputs/set5/merged"
    # )