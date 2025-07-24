"""
parse_patent_to_fasta.py
Command-line tool to parse amino acid sequences from patent text files and convert them to FASTA format.

Overview:
    - Reads patent sequence files containing amino acid sequences in 3-letter codes.
    - Converts 3-letter amino acid codes to 1-letter codes using a standard mapping.
    - Handles invalid or ambiguous codes by prompting the user for correction or skipping.
    - Outputs a FASTA file with the parsed and converted sequences.

Arguments:
    --input (str): Path to the input file containing patent sequences (in 3-letter code format).
    --output (str): Path to the output FASTA file to be generated.

Requirements:
    - Input file must be present at the specified path.
    - Python packages: re, argparse.

Usage Example:
    python parse_patent_to_fasta.py --input patent_sequences.txt --output patent_sequences.fasta

Outputs:
    <output_fasta>   # FASTA file containing converted amino acid sequences

Author: Nadia
"""
import re
import argparse

# 3-letter to 1-letter amino acid codes
AMINO_ACID_MAP = {
    'Ala': 'A', 'Arg': 'R', 'Asn': 'N', 'Asp': 'D', 'Cys': 'C',
    'Gln': 'Q', 'Glu': 'E', 'Gly': 'G', 'His': 'H', 'Ile': 'I',
    'Leu': 'L', 'Lys': 'K', 'Met': 'M', 'Phe': 'F', 'Pro': 'P',
    'Ser': 'S', 'Thr': 'T', 'Trp': 'W', 'Tyr': 'Y', 'Val': 'V',
    'Sec': 'U', 'Pyl': 'O', 'Asx': 'B', 'Glx': 'Z', 'Xaa': 'X',
    'Ter': '*', 'Trp*': 'W'
}

def prompt_for_correction(invalid_code):
    """Prompt the user to correct an invalid amino acid code.
    Args:
        invalid_code (str): The invalid 3-letter amino acid code.
    Returns:
        str: The corrected 1-letter code or None if the user chooses to skip.
    This function prompts the user to enter a valid 3-letter code or skip the invalid code.
    It checks the input against the AMINO_ACID_MAP and returns the corresponding 1-letter code.
    If the input is not valid, it continues to prompt until a valid code is entered
    or the user chooses to skip."""
    while True:
        user_input = input(f"Invalid amino acid code '{invalid_code}'. Enter a valid 3-letter code or type 'skip': ").strip()
        if user_input.lower() == 'skip':
            return None
        elif user_input in AMINO_ACID_MAP:
            return AMINO_ACID_MAP[user_input]
        else:
            print("Invalid input. Try again.")

def parse_sequences(input_file, output_fasta):
    """Parse patent sequences from a text file and convert to FASTA format.
    Args:
        input_file (str): Path to the input file containing patent sequences.
        output_fasta (str): Path to the output FASTA file.
    """
    with open(input_file, 'r') as infile:
        lines = infile.readlines()

    sequences = []
    current_seq = []
    current_id = None
    in_sequence_block = False

    for line in lines:
        line = line.strip()
        if line.startswith('SEQ ID NO'):
            if current_id and current_seq:
                sequences.append((current_id, current_seq))
                current_seq = []
            current_id = line.replace('SEQ ID NO', '').strip()
            in_sequence_block = False

        elif line.startswith('SEQUENCE:'):
            in_sequence_block = True
            content = line.split(':', 1)[1].strip()
            if content:
                current_seq.extend(content.split())

        elif in_sequence_block:
            if re.match(r'^SEQ ID NO', line):
                in_sequence_block = False
                continue
            if line.startswith('---') or not line:
                continue
            current_seq.extend(line.split())

    if current_id and current_seq:
        sequences.append((current_id, current_seq))

    with open(output_fasta, 'w') as fasta_out:
        for seq_id, amino_acids in sequences:
            one_letter_seq = ""
            for aa in amino_acids:
                if aa.isdigit():
                    # Automatically skip numeric entries
                    print(f"Skipping numeric code '{aa}' in SEQ ID NO {seq_id}")
                    continue
                elif aa in AMINO_ACID_MAP:
                    one_letter_seq += AMINO_ACID_MAP[aa]
                else:
                    corrected = prompt_for_correction(aa)
                    if corrected:
                        one_letter_seq += corrected
                    else:
                        print(f"Skipping invalid code '{aa}' in SEQ ID NO {seq_id}")

            fasta_out.write(f">SEQ_ID_{seq_id}\n")
            fasta_out.write(one_letter_seq + '\n')

    print(f"\n✅ FASTA file saved as '{output_fasta}'")

if __name__ == "__main__":
    """ Command-line interface for parsing patent amino acid sequences and converting to FASTA format.
    Usage:
        python parse_patent_to_fasta.py --input <input_file> --output <output_fasta>
    """ 
    parser = argparse.ArgumentParser(
        description="Parse patent amino acid sequences and convert to FASTA format.",
        usage="python parse_patent_to_fasta.py --input <input_file> --output <output_fasta>"
    )
    parser.add_argument("--input", required=True, help="Input file containing patent sequences")
    parser.add_argument("--output", required=True, help="Output FASTA file")
    args = parser.parse_args()

    parse_sequences(args.input, args.output)
