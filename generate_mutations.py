from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from itertools import product

def generate_mutations(input_file_name, input_positions, output_file_name):
    """
    Parameters (inputs):
    - input_file_name: The name of the FASTA file containing the scaffold protein sequence.
    - input_positions: A list of positions (each indexed from 1) in the protein sequence where the residue will be mutated
    - output_file_name: The name of the output FASTA file that will contain all the mutated sequences.

    Outputs:
    - A FASTA file containing all the possible mutations, where each mutated sequence has a unique name.
    """
    positions = [position - 1 for position in input_positions]  # Convert to 0-based indexing
    amino_acids = list("ACDEFGHIKLMNPQRSTVWY")

    # Read the scaffold protein sequence
    scaffold_record = SeqIO.read(input_file_name, "fasta")
    scaffold_sequence = str(scaffold_record.seq)

    # Generate all possible mutation combinations at the specified positions
    mutation_options = [amino_acids] * len(positions)  # Create a list of lists; each inner list is amino acid options for each position
    all_combinations = product(*mutation_options)  # List of all possible combinations of mutations at the specified positions
    
    # Apply mutations to the scaffold sequence and create SeqRecords for each mutant
    mutant_records = []
    for index, combo in enumerate(all_combinations):
        mutated_sequence = list(scaffold_sequence)  # Scaffold sequence (that will be mutated) converted to a list for mutability

        for position, amino_acid in zip(positions, combo):
            mutated_sequence[position] = amino_acid  # Mutate the specified position
        mutated_sequence = "".join(mutated_sequence)  # Convert back to string
        
        # Create a SeqRecord for the mutant sequence
        mutant_record = SeqRecord(
            Seq(mutated_sequence),
            id=f"scaffold_0_mutant_{index}",
            description=""
        )
        mutant_records.append(mutant_record)
    
    # Write all mutated sequences to a FASTA file
    SeqIO.write(mutant_records, output_file_name, "fasta")


if __name__ == "__main__":
    # Change inputs if needed
    input_file_name = "scaffold_0.fasta"
    output_file_name = "mutants_test.fasta"
    positions = [1]  # List of positions to mutate (1-based indexing)

    generate_mutations(input_file_name, positions, output_file_name)