import re
import numpy as np
import matplotlib.pyplot as plt

ENZYMES = {
    "EcoRI":   ("GAATTC", 1),
    "BamHI":   ("GGATCC", 1),
    "HindIII": ("AAGCTT", 1),
    "TaqI":    ("TCGA",   1),
    "HaeIII":  ("GGCC",   2),
}

def load_fasta(filepath: str) -> str:
    seq_lines = []
    with open(filepath, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line.startswith(">"):
                seq_lines.append(line.upper())
    return "".join(seq_lines)

def digest_dna(dna_seq: str, motif: str, offset: int) -> tuple[list[int], list[int]]:

    dna_seq = dna_seq.upper()
    motif = motif.upper()
    
    cut_indices = [m.start() + offset for m in re.finditer(f"(?={motif})", dna_seq)]
    
    boundaries = [0] + cut_indices + [len(dna_seq)]
    fragments = [boundaries[i + 1] - boundaries[i] for i in range(len(boundaries) - 1)]
    
    cut_positions_1based = [idx + 1 for idx in cut_indices]
    return cut_positions_1based, fragments

def simulate_gel(digest_results: dict[str, list[int]], total_length: int):
    ladder = [10000, 8000, 6000, 5000, 4000, 3000, 2000, 1500, 1000, 750, 500, 250, 100]
    lanes = {"Ladder": ladder}
    lanes.update(digest_results)
    
    lane_names = list(lanes.keys())
    num_lanes = len(lane_names)

    fig, ax = plt.subplots(figsize=(2.2 * num_lanes, 7))
    
    ax.set_facecolor("#121417")
    fig.patch.set_facecolor("#1e2229")

    min_bp, max_bp = 50, max(12000, total_length * 1.2)
    log_min, log_max = np.log10(min_bp), np.log10(max_bp)

    def bp_to_y(bp: int) -> float:
        bp_clamped = max(min_bp, min(max_bp, bp))
        return (log_max - np.log10(bp_clamped)) / (log_max - log_min)

    for col_idx, (lane_name, fragments) in enumerate(lanes.items()):
        x_center = col_idx + 1
        
        ax.add_patch(plt.Rectangle((x_center - 0.3, -0.04), 0.6, 0.03, color="#2c323d", ec="#6c757d"))
        
        for bp in fragments:
            if bp <= 0:
                continue
            y_pos = bp_to_y(bp)
            
            alpha = min(0.95, max(0.35, np.log10(bp) / 4.2))
            lw = 3.2 if lane_name != "Ladder" else 2.2
            color = "#39ff14" if lane_name != "Ladder" else "#00d4ff"
            
            ax.hlines(y=y_pos, xmin=x_center - 0.28, xmax=x_center + 0.28,
                      colors=color, linewidth=lw, alpha=alpha)
            
            if lane_name == "Ladder":
                ax.text(x_center - 0.34, y_pos, f"{bp} bp", color="#adb5bd",
                        fontsize=8, va="center", ha="right")


    ax.set_xlim(0.2, num_lanes + 0.6)
    ax.set_ylim(1.05, -0.08)
    ax.set_xticks(range(1, num_lanes + 1))
    ax.set_xticklabels(lane_names, color="white", fontsize=11, fontweight="bold")
    ax.xaxis.tick_top()
    ax.set_yticks([])
    
    for spine in ax.spines.values():
        spine.set_edgecolor("#495057")

    plt.title("Simulated Agarose Gel Electrophoresis (Restriction Digest)",
              color="white", fontsize=13, pad=20, fontweight="bold")
    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    dna_sequence = load_fasta("sequence.fasta")
    

    print(f"Total DNA Sequence Length: {len(dna_sequence)} bp\n" + "=" * 55)
    
    gel_lanes = {}
    for enzyme_name, (motif, offset) in ENZYMES.items():
        positions, frag_lengths = digest_dna(dna_sequence, motif, offset)
        gel_lanes[enzyme_name] = frag_lengths
        
        print(f"Enzyme: {enzyme_name} (Site: {motif[:offset]}|{motif[offset:]})")
        print(f"  -> Number of cleavages : {len(positions)}")
        print(f"  -> Cleavage positions  : {positions if positions else 'None'}")
        print(f"  -> Fragment lengths    : {frag_lengths} bp\n")

    simulate_gel(gel_lanes, len(dna_sequence))