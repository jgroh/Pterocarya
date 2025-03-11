import random
import argparse

# Dictionary mapping IUPAC ambiguous codes (uppercase and lowercase) to possible nucleotides
IUPAC_CODES = {
    'R': ['A', 'G'], 'r': ['A', 'G'],
    'Y': ['C', 'T'], 'y': ['C', 'T'],
    'S': ['G', 'C'], 's': ['G', 'C'],
    'W': ['A', 'T'], 'w': ['A', 'T'],
    'K': ['G', 'T'], 'k': ['G', 'T'],
    'M': ['A', 'C'], 'm': ['A', 'C'],
    'B': ['C', 'G', 'T'], 'b': ['C', 'G', 'T'],
    'D': ['A', 'G', 'T'], 'd': ['A', 'G', 'T'],
    'H': ['A', 'C', 'T'], 'h': ['A', 'C', 'T'],
    'V': ['A', 'C', 'G'], 'v': ['A', 'C', 'G'],
    'N': ['A', 'C', 'G', 'T'], 'n': ['A', 'C', 'G', 'T']
}

def replace_ambiguous(sequence):
    """
    Replace IUPAC ambiguous codes in the sequence with one of the corresponding nucleotides at random.
    """
    result = []
    for base in sequence:
        if base in IUPAC_CODES:
            result.append(random.choice(IUPAC_CODES[base]))
        else:
            result.append(base)
    return ''.join(result)

def process_fasta(input_file):
    """
    Process a FASTA file, replacing IUPAC ambiguous codes in sequences
    while preserving headers, and write to standard output.
    """
    with open(input_file, 'r') as infile:
        for line in infile:
            if line.startswith(">"):
                # Write the header as is
                print(line.strip())
            else:
                # Process and write the sequence
                clean_sequence = replace_ambiguous(line.strip())
                print(clean_sequence)

if __name__ == "__main__":
    # Argument parser setup
    parser = argparse.ArgumentParser(description="Replace IUPAC ambiguous codes in a FASTA file.")
    parser.add_argument("input_file", help="Input FASTA file")
    
    # Parse arguments
    args = parser.parse_args()
    
    # Process the FASTA file and output to standard out
    process_fasta(args.input_file)
