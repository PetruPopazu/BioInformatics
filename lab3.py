from collections import Counter

GENETIC_CODE = {
    'UUU': 'Phe (F)', 'UUC': 'Phe (F)', 'UUA': 'Leu (L)', 'UUG': 'Leu (L)',
    'UCU': 'Ser (S)', 'UCC': 'Ser (S)', 'UCA': 'Ser (S)', 'UCG': 'Ser (S)',
    'UAU': 'Tyr (Y)', 'UAC': 'Tyr (Y)', 'UAA': 'Stop (*)', 'UAG': 'Stop (*)',
    'UGU': 'Cys (C)', 'UGC': 'Cys (C)', 'UGA': 'Stop (*)', 'UGG': 'Trp (W)',
    'CUU': 'Leu (L)', 'CUC': 'Leu (L)', 'CUA': 'Leu (L)', 'CUG': 'Leu (L)',
    'CCU': 'Pro (P)', 'CCC': 'Pro (P)', 'CCA': 'Pro (P)', 'CCG': 'Pro (P)',
    'CAU': 'His (H)', 'CAC': 'His (H)', 'CAA': 'Gln (Q)', 'CAG': 'Gln (Q)',
    'CGU': 'Arg (R)', 'CGC': 'Arg (R)', 'CGA': 'Arg (R)', 'CGG': 'Arg (R)',
    'AUU': 'Ile (I)', 'AUC': 'Ile (I)', 'AUA': 'Ile (I)', 'AUG': 'Met (M)',
    'ACU': 'Thr (T)', 'ACC': 'Thr (T)', 'ACA': 'Thr (T)', 'ACG': 'Thr (T)',
    'AAU': 'Asn (N)', 'AAC': 'Asn (N)', 'AAA': 'Lys (K)', 'AAG': 'Lys (K)',
    'AGU': 'Ser (S)', 'AGC': 'Ser (S)', 'AGA': 'Arg (R)', 'AGG': 'Arg (R)',
    'GUU': 'Val (V)', 'GUC': 'Val (V)', 'GUA': 'Val (V)', 'GUG': 'Val (V)',
    'GCU': 'Ala (A)', 'GCC': 'Ala (A)', 'GCA': 'Ala (A)', 'GCG': 'Ala (A)',
    'GAU': 'Asp (D)', 'GAC': 'Asp (D)', 'GAA': 'Glu (E)', 'GAG': 'Glu (E)',
    'GGU': 'Gly (G)', 'GGC': 'Gly (G)', 'GGA': 'Gly (G)', 'GGG': 'Gly (G)'
}

def read_fasta(file_path):
    sequence_parts = []
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                continue
            sequence_parts.append(line.strip().upper())
    return "".join(sequence_parts)

def extract_codons(sequence, step=3):
    seq = "".join([base for base in sequence if base in "ACGT"])
    
    codons = []
    for i in range(0, len(seq) - 2, step):
        codon_dna = seq[i:i+3]
        if len(codon_dna) == 3:
            codons.append(codon_dna)
            
    return Counter(codons)

def dna_to_rna(dna_codon):
    return dna_codon.replace('T', 'U')


fasta_file = "genome.fasta"
    
try:
    genome_seq = read_fasta(fasta_file)
except FileNotFoundError:
    print(f"Error: The file '{fasta_file}' was not found. Please check the path.")


codon_counts = extract_codons(genome_seq, step=3)
    
if not codon_counts:
    print("No valid codons found in the sequence.")

print("Top 10 Highest Frequency Codons:")
print("-" * 65)
print(f"{'Rank':<6}{'DNA Codon':<12}{'RNA Codon':<12}{'Frequency':<12}{'Amino Acid':<15}")
print("-" * 65)

for rank, (dna_codon, freq) in enumerate(codon_counts.most_common(10), 1):
    rna_codon = dna_to_rna(dna_codon)
    amino_acid = GENETIC_CODE.get(rna_codon, 'Unknown')
    print(f"{rank:<6}{dna_codon:<12}{rna_codon:<12}{freq:<12}{amino_acid:<15}")
print("-" * 65)