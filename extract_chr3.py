from Bio import SeqIO
from pathlib import Path

def extract_chromosome(input_file, output_file, chr_id="NC_001135.5"):
    """Extract chromosome III from the genome file"""
    with open(input_file, 'r') as fin, open(output_file, 'w') as fout:
        for record in SeqIO.parse(fin, 'fasta'):
            if chr_id in record.id:
                SeqIO.write(record, fout, 'fasta')
                print(f"Extracted {record.id}, length: {len(record.seq)} bp")
                return True
    return False

if __name__ == "__main__":
    genome_file = "r64/GCF_000146045.2_R64_genomic.fasta"
    output_file = "chr3.fasta"
    
    if extract_chromosome(genome_file, output_file):
        print(f"Successfully extracted chromosome III to {output_file}")
    else:
        print("Failed to find chromosome III") 