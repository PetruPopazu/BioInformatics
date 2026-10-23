import random

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

def read_fasta(file_path, max_bases=20000):
    sequence_parts = []
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                continue
            sequence_parts.append(line.strip().upper())
            
    raw_seq = "".join(sequence_parts)
    clean_seq = "".join([base for base in raw_seq if base in "ACGT"])
    return clean_seq[:max_bases]

def generate_overlapping_samples(genome, a=100, b=150, min_coverage=10):
    n = len(genome)
    coverage = [0] * n
    reads = []
    
    print(f"[+] Sampling genome (Target length: {n:,} bp) into reads of length [{a}, {b}]...")
    
    while any(c < min_coverage for c in coverage):
        read_len = random.randint(a, min(b, n))
        start_pos = random.randint(0, n - read_len)
        end_pos = start_pos + read_len
        
        read = genome[start_pos:end_pos]
        reads.append(read)
        
        for i in range(start_pos, end_pos):
            coverage[i] += 1
            
    random.shuffle(reads)
    print(f"[+] Generated {len(reads):,} reads. Guaranteed min coverage depth >= {min_coverage}x.")
    return reads

def assemble_genome(reads, min_overlap=10):
    print("[+] Rebuilding genome from read collection...")
    contigs = set(reads)
    
    while True:
        merged_any = False
        prefix_map = {}
        
        for c in contigs:
            if len(c) >= min_overlap:
                p = c[:min_overlap]
                prefix_map.setdefault(p, []).append(c)
                
        for c1 in list(contigs):
            if c1 not in contigs:
                continue
                
            max_check = min(len(c1), 150)
            for ov_len in range(max_check, min_overlap - 1, -1):
                suffix_key = c1[-ov_len:][:min_overlap]
                
                if suffix_key in prefix_map:
                    found_merge = False
                    for c2 in prefix_map[suffix_key]:
                        if c1 != c2 and c2 in contigs and c1.endswith(c2[:ov_len]):
                            contigs.remove(c1)
                            contigs.remove(c2)
                            contigs.add(c1 + c2[ov_len:])
                            merged_any = True
                            found_merge = True
                            break
                    if found_merge:
                        break
                        
        if not merged_any:
            break
            
    sorted_contigs = sorted(contigs, key=len, reverse=True)
    return sorted_contigs[0] if sorted_contigs else ""


def main():
    fasta_file = "genome.fasta"
    
    try:
        target_genome = read_fasta(fasta_file, max_bases=20000)
        print(f"Successfully loaded '{fasta_file}'. Working slice length: {len(target_genome):,} bp.")
    except FileNotFoundError:
        print(f"Error: '{fasta_file}' not found. Using test genome.")
        target_genome = "ATGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGA" * 100

    a, b = 100, 150
    reads = generate_overlapping_samples(target_genome, a=a, b=b, min_coverage=10)
    rebuilt_genome = assemble_genome(reads, min_overlap=10)
    
    print("\n" + "=" * 50)
    print("RESULTS")
    print("=" * 50)
    print(f"Original Target Length: {len(target_genome)} bp")
    print(f"Rebuilt Genome Length : {len(rebuilt_genome)} bp")
    
    
    if target_genome == rebuilt_genome:
        print("Genome Rebuilding: SUCCESSFUL (100% Sequence Identity)")
    else:
        print("Genome Rebuilding: PARTIAL (Contains repeat ambiguities)")


if __name__ == "__main__":
    main()