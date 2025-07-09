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
    while True:
        user_input = input(f"Invalid amino acid code '{invalid_code}'. Enter a valid 3-letter code or type 'skip': ").strip()
        if user_input.lower() == 'skip':
            return None
        elif user_input in AMINO_ACID_MAP:
            return AMINO_ACID_MAP[user_input]
        else:
            print("Invalid input. Try again.")

def parse_sequences(input_file, output_fasta):
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
