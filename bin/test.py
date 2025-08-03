from Bio import AlignIO
import os

def find_bad_clustal_files(dir_path):
    print(f"Checking files in {dir_path} for Clustal alignment issues...")
    if not os.path.isdir(dir_path):
        print(f"[ERROR] The directory {dir_path} does not exist.")
        return
    for filename in os.listdir(dir_path):
        if not filename.endswith("_aligned.fasta"):
            continue
        filepath = os.path.join(dir_path, filename)
        try:
            AlignIO.read(filepath, "clustal")
        except Exception as e:
            print(f"[BAD] {filename} — {type(e).__name__}: {e}")
        else:
            print(f"[OK]  {filename}")

find_bad_clustal_files("data/S.aureus/epitope_outputs/mafft_rate4site/")
