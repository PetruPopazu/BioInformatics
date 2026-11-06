import random

def reverse_complement(seq):
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return ''.join(complement[base] for base in reversed(seq))

def random_dna(length):
    return ''.join(random.choice('ACGT') for _ in range(length))

def build_transposon(left_ir, middle):
    return left_ir + middle + reverse_complement(left_ir)

def insert_transposon(sequence, te_sequence, position, tsd_length=5):
    tsd = sequence[position:position + tsd_length]
    new_sequence = sequence[:position] + tsd + te_sequence + tsd + sequence[position + tsd_length:]
    te_start = position + tsd_length
    te_end = te_start + len(te_sequence)
    return new_sequence, te_start, te_end

def detect_transposons(sequence, ir_length=8, tsd_length=5, min_te_length=25, max_te_length=90):
    detections = []
    n = len(sequence)

    for i in range(tsd_length, n - ir_length - tsd_length):
        left_tsd = sequence[i - tsd_length : i]
        left_ir = sequence[i : i + ir_length]
        expected_right_ir = reverse_complement(left_ir)

        search_from = i + min_te_length - ir_length
        search_to = min(i + max_te_length - ir_length, n - ir_length - tsd_length)

        for j in range(search_from, search_to + 1):
            if sequence[j : j + ir_length] == expected_right_ir:
                right_tsd = sequence[j + ir_length : j + ir_length + tsd_length]
                if left_tsd == right_tsd:
                    detections.append((i - tsd_length, j + ir_length + tsd_length))
                    break

    return list(dict.fromkeys(detections))

random.seed(42)
sequence = random_dna(350)
print("Initial size of DNA Sequence:",len(sequence),"\n")

te_definitions = [
    ("ATCGATCG", "GGGGCCCCAAAATTTT"),
    ("GCTAGCTA", "TTTTAAAACCCCGGGG"),
    ("TAGCTAGC", "AAAAGGGGCCCTTTTT"),
    ("CGATCGAT", "CCCCGGGGTTTTAAAA"),
]
te_sequences = [build_transposon(left, middle) for left, middle in te_definitions]

sequence, _, _ = insert_transposon(sequence, te_sequences[0], 50, 5)
sequence, s_nested, _ = insert_transposon(sequence, te_sequences[1], 62, 4)
sequence, s3, _ = insert_transposon(sequence, te_sequences[2], 180, 5)
sequence, _, _ = insert_transposon(sequence, te_sequences[3], s3 + len(te_sequences[2]) - 8, 5)
    
print("DNA Sequence with inserted transponsons:\n",sequence)
print()
print("Detected transposable elements:")
detected = detect_transposons(sequence)
for start, end in detected:
    print(f"start:{start} end:{end}")
