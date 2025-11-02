import pandas as pd

def update_csv_with_fasta(csv_path, fasta_path, output_csv_path):
    # --- Step 1: Read the CSV ---
    df = pd.read_csv(csv_path)

    # Ensure the nucleotide_sequence column exists
    if 'nucleotide_sequence' not in df.columns:
        df['nucleotide_sequence'] = ''

    # --- Step 2: Parse the FASTA ---
    fasta_dict = {}
    with open(fasta_path, 'r') as f:
        accession = None
        seq_lines = []
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                # Save previous record
                if accession and seq_lines:
                    fasta_dict[accession] = ''.join(seq_lines)
                accession = line[1:].split()[0]  # take first token as accession
                seq_lines = []
            else:
                seq_lines.append(line)
        # Save the last one
        if accession and seq_lines:
            fasta_dict[accession] = ''.join(seq_lines)

    # --- Step 3: Update the CSV dataframe ---
    for i, row in df.iterrows():
        acc = row['uniprot_accession']
        if acc in fasta_dict:
            df.at[i, 'nucleotide_sequence'] = fasta_dict[acc]

    # --- Step 4: Save the updated CSV ---
    df.to_csv(output_csv_path, index=False)
    print(f"✅ Updated CSV saved to {output_csv_path}")


# Example usage:
update_csv_with_fasta('data/S.aureus_analysis/staphylococcus_aureus_compiled_proteins.csv', 'data/S.aureus_analysis/reverse_translated_antigen_sequences.fasta', 'data/S.aureus_analysis/staphylococcus_aureus_compiled_proteins_corrected.csv')

