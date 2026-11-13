import random
import numpy as np
import matplotlib.pyplot as plt

def load_fasta_slice(filepath: str, max_len: int = 2500) -> str:
    lines = []
    with open(filepath, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line.startswith(">"):
                lines.append(line.upper())
    seq = "".join(lines)[:max_len]
    if not (1000 <= len(seq) <= 3000):
        raise ValueError(f"Sequence length ({len(seq)}) must be between 1000 and 3000 nt.")
    return seq

def sample_random_fragments(dna_seq: str, num_samples: int = 10, min_bp: int = 100, max_bp: int = 3000) -> list[str]:
    seq_len = len(dna_seq)
    upper_bound = min(max_bp, seq_len)
    samples = []
    
    for _ in range(num_samples):
        frag_len = random.randint(min_bp, upper_bound)
        start_idx = random.randint(0, seq_len - frag_len)
        samples.append(dna_seq[start_idx : start_idx + frag_len])
        
    return samples

def simulate_samples_gel(samples: list[str], seq_len: int):
    ladder = [3000, 2500, 2000, 1500, 1000, 750, 500, 250, 100]
    
    fig, ax = plt.subplots(figsize=(13, 7))
    ax.set_facecolor("#121417")
    fig.patch.set_facecolor("#1e2229")

    min_bp, max_bp = 80, max(3200, seq_len + 200)
    log_min, log_max = np.log10(min_bp), np.log10(max_bp)

    def bp_to_y(bp: int) -> float:
        bp_clamped = max(min_bp, min(max_bp, bp))
        return (log_max - np.log10(bp_clamped)) / (log_max - log_min)

    ax.add_patch(plt.Rectangle((1 - 0.3, -0.04), 0.6, 0.03, color="#2c323d", ec="#6c757d"))
    for bp in ladder:
        y = bp_to_y(bp)
        ax.hlines(y=y, xmin=1 - 0.28, xmax=1 + 0.28, colors="#00d4ff", linewidth=2.2)
        ax.text(1 - 0.34, y, f"{bp} bp", color="#adb5bd", fontsize=8, va="center", ha="right")

    for idx, frag in enumerate(samples, start=1):
        lane_x = idx + 1
        frag_len = len(frag)
        y = bp_to_y(frag_len)
        
        ax.add_patch(plt.Rectangle((lane_x - 0.3, -0.04), 0.6, 0.03, color="#2c323d", ec="#6c757d"))
        ax.hlines(y=y, xmin=lane_x - 0.28, xmax=lane_x + 0.28, colors="#39ff14", linewidth=3.5)
        ax.text(lane_x, y - 0.02, f"{frag_len}", color="#e9ecef", fontsize=7.5, ha="center")

    labels = ["Ladder"] + [f"S{i}" for i in range(1, len(samples) + 1)]
    ax.set_xlim(0.1, len(samples) + 1.7)
    ax.set_ylim(1.05, -0.08)
    ax.set_xticks(range(1, len(samples) + 2))
    ax.set_xticklabels(labels, color="white", fontsize=10, fontweight="bold")
    ax.xaxis.tick_top()
    ax.set_yticks([])

    for spine in ax.spines.values():
        spine.set_edgecolor("#495057")

    plt.title("Exercise 1: Simulated Agarose Gel Electrophoresis (10 Random DNA Samples)",
              color="white", fontsize=13, pad=20, fontweight="bold")
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":

    dna_seq = load_fasta_slice("sequence.fasta", max_len=2400)
    
    samples_array = sample_random_fragments(dna_seq, num_samples=10, min_bp=100, max_bp=3000)
    
    print(f"Parent NCBI Sequence Length: {len(dna_seq)} bp\n" + "-" * 50)
    for i, s in enumerate(samples_array, start=1):
        print(f"Sample {i:2d}: {len(s):4d} bp | Preview: {s[:25]}...")
        
    simulate_samples_gel(samples_array, len(dna_seq))