
pip install matplotlib



import matplotlib.pyplot as plt
from collections import Counter

# Reverse Complement Tool
def reverse_complement(dna_seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    dna_seq = dna_seq.upper().replace(' ', '')  # Remove spaces and ensure uppercase
    try:
        return ''.join(complement[base] for base in dna_seq[::-1])
    except KeyError as e:
        return f"Error: Invalid character '{e.args[0]}' found in the DNA sequence."

# GC Content Calculator
def gc_content(dna_seq):
    dna_seq = dna_seq.upper().replace(' ', '')  # Clean input
    gc_count = dna_seq.count('G') + dna_seq.count('C')
    return (gc_count / len(dna_seq)) * 100 if len(dna_seq) > 0 else 0

# Transcription Simulator (DNA to RNA)
def transcribe(dna_seq):
    return dna_seq.upper().replace('T', 'U')

# RNA to Protein Translation Tool
genetic_code = {
    'AUG': 'M', 'UGG': 'W', 'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
    # Add rest of the codons
    'UAA': 'Stop', 'UAG': 'Stop', 'UGA': 'Stop'
}

def translate_rna_to_protein(rna_seq):
    protein_seq = []
    for i in range(0, len(rna_seq), 3):
        codon = rna_seq[i:i+3]
        if len(codon) != 3:
            continue
        amino_acid = genetic_code.get(codon, '-')
        if amino_acid == 'Stop':
            break
        protein_seq.append(amino_acid)
    return ''.join(protein_seq)

# Codon Usage Calculator
def codon_usage(dna_seq):
    dna_seq = dna_seq.upper().replace(' ', '')  # Clean input
    codon_count = Counter([dna_seq[i:i+3] for i in range(0, len(dna_seq)-2, 3)])
    return codon_count

#def plot_codon_usage(codon_count):
#    labels, values = zip(*codon_count.items())
#    plt.bar(labels, values)
#    plt.xlabel('Codons')
#    plt.ylabel('Frequency')
#    plt.title('Codon Usage')
#    plt.xticks(rotation=90)
#    plt.show()

def plot_codon_usage(codon_count):
    labels, values = zip(*codon_count.items())
    
    # Increase figure size for better readability
    plt.figure(figsize=(10, 6))  # Adjust the size as per your requirement
    
    plt.bar(labels, values)
    plt.xlabel('Codons')
    plt.ylabel('Frequency')
    plt.title('Codon Usage')
    
    # Rotate x-ticks, adjust alignment, and reduce font size
    plt.xticks(rotation=90, ha='right', fontsize=8)
    
    plt.tight_layout()  # Automatically adjust subplot parameters for a nice fit
    plt.show()




# Open Reading Frame (ORF) Finder
def find_orfs(dna_seq):
    start_codon = 'ATG'
    stop_codons = ['TAA', 'TAG', 'TGA']
    orfs = []
    for i in range(len(dna_seq)):
        if dna_seq[i:i+3] == start_codon:
            for j in range(i, len(dna_seq), 3):
                if dna_seq[j:j+3] in stop_codons:
                    orfs.append(dna_seq[i:j+3])
                    break
    return orfs

# k-mer Frequency Counter
def kmer_frequency(dna_seq, k=3):
    dna_seq = dna_seq.upper().replace(' ', '')  # Clean input
    kmers = [dna_seq[i:i+k] for i in range(len(dna_seq) - k + 1)]
    return Counter(kmers)

# Motif Search
def search_motif(dna_seq, motif):
    dna_seq = dna_seq.upper().replace(' ', '')  # Clean input
    positions = []
    for i in range(len(dna_seq) - len(motif) + 1):
        if dna_seq[i:i+len(motif)] == motif:
            positions.append(i)
    return positions

# Main Pipeline
def dna_analysis_pipeline(dna_seq, motif):
    dna_seq = dna_seq.upper().replace(' ', '')  # Clean input
    
    print(f"Original DNA Sequence: {dna_seq}\n")
    
    # Reverse Complement
    rev_complement = reverse_complement(dna_seq)
    print(f"Reverse Complement: {rev_complement}\n")
    
    # GC Content
    gc_percentage = gc_content(dna_seq)
    print(f"GC Content: {gc_percentage:.2f}%\n")
    
    # Transcription
    rna_seq = transcribe(dna_seq)
    print(f"RNA Sequence: {rna_seq}\n")
    
    # Translation
    protein_seq = translate_rna_to_protein(rna_seq)
    print(f"Protein Sequence: {protein_seq}\n")
    
    # Codon Usage
    codon_count = codon_usage(dna_seq)
    print("Codon Usage Frequencies:")
    print(codon_count)
    plot_codon_usage(codon_count)
    
    # ORF Finder
    orfs = find_orfs(dna_seq)
    print(f"Open Reading Frames (ORFs): {orfs}\n")
    
    # k-mer Frequency
    kmer_freq = kmer_frequency(dna_seq)
    print(f"k-mer Frequencies: {kmer_freq}\n")
    
    # Motif Search
    motif_positions = search_motif(dna_seq, motif)
    print(f"Motif '{motif}' found at positions: {motif_positions}\n")


if __name__ == "__main__":
    dna_seq = input("Enter your DNA sequence: ").strip()
    motif = input("Enter the motif to search for: ").strip()
    dna_analysis_pipeline(dna_seq, motif)
