def load_fasta(filepath: str) -> str:
    seq_lines = []
    with open(filepath, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line.startswith(">"):
                seq_lines.append(line.upper())
    return "".join(seq_lines)
def detect_tandem_repeats(dna_seq: str, min_unit: int = 3, max_unit: int = 5,
                          min_total: int = 6, max_total: int = 10) -> list[dict]:
    dna_seq = dna_seq.upper()
    n = len(dna_seq)
    results = []

    for u in range(min_unit, max_unit + 1):
        i = 0
        while i <= n - min_total:
            unit = dna_seq[i : i + u]
            
            if len(set(unit)) == 1:
                i += 1
                continue

            copies = 1
            while i + (copies + 1) * u <= n and dna_seq[i + copies * u : i + (copies + 1) * u] == unit:
                copies += 1

            total_len = copies * u
            if copies >= 2 and min_total <= total_len <= max_total:
                if i < u or dna_seq[i - u : i] != unit:
                    results.append({
                        "start_1based": i + 1,
                        "end_1based": i + total_len,
                        "unit": unit,
                        "unit_len": u,
                        "copies": copies,
                        "total_len": total_len,
                        "sequence": dna_seq[i : i + total_len]
                    })
                i += total_len
            else:
                i += 1

    results.sort(key=lambda r: (r["start_1based"], r["unit_len"]))
    return results


if __name__ == "__main__":
    dna_seq = load_fasta("sequence.fasta")

    repeats = detect_tandem_repeats(dna_seq)

    print(f"Sequence Length: {len(dna_seq)} bp")
    print(f"Total Tandem Repeats Found (Unit: 3-5 bp | Total: 6-10 bp): {len(repeats)}\n")
    print(f"{'Pos (1-based)':<16} | {'Unit':<6} | {'Unit Len':<8} | {'Copies':<6} | {'Total Len':<9} | {'Repeated Sequence'}")
    print("-" * 75)
    
    for r in repeats:
        pos_str = f"{r['start_1based']}..{r['end_1based']}"
        print(f"{pos_str:<16} | {r['unit']:<6} | {r['unit_len']:<8} | {r['copies']:<6} | {r['total_len']:<9} | {r['sequence']}")