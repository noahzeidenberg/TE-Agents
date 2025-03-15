import os
import requests
from pathlib import Path

# Create consensus directory if it doesn't exist
Path("consensus").mkdir(exist_ok=True)

# SGD TE consensus sequences
te_sequences = {
    "Ty1": "https://sgd-dev-upload.s3.amazonaws.com/S000007345/S000007345.fsa",
    "Ty2": "https://sgd-dev-upload.s3.amazonaws.com/S000007346/S000007346.fsa",
    "Ty3": "https://sgd-dev-upload.s3.amazonaws.com/S000007347/S000007347.fsa",
    "Ty4": "https://sgd-dev-upload.s3.amazonaws.com/S000007348/S000007348.fsa",
    "Ty5": "https://sgd-dev-upload.s3.amazonaws.com/S000007349/S000007349.fsa"
}

def download_sequence(name, url):
    """Download a sequence file from SGD"""
    print(f"Downloading {name} consensus sequence...")
    response = requests.get(url)
    if response.status_code == 200:
        output_file = f"consensus/{name}.fasta"
        with open(output_file, 'w') as f:
            f.write(response.text)
        print(f"Saved to {output_file}")
    else:
        print(f"Failed to download {name}: {response.status_code}")

def main():
    for name, url in te_sequences.items():
        download_sequence(name, url)
    print("Done downloading consensus sequences.")

if __name__ == "__main__":
    main() 