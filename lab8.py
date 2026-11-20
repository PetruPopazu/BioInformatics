#Implement a software app that alligns two DNA sequences by using the Needleman-Wunsch algorithm:
#s1="ACCGTGAAGCCAATAC"
#s2="AGCGTGCAGCCAATAC"

def needleman_wunsch(seq1, seq2, match=2, mismatch=-1, gap=-2):
    n = len(seq1)
    m = len(seq2)
    
    score = [[0] * (m + 1) for _ in range(n + 1)]
    
    for i in range(1, n + 1):
        score[i][0] = i * gap
    for j in range(1, m + 1):
        score[0][j] = j * gap
    
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            diag = score[i-1][j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch)
            up   = score[i-1][j] + gap
            left = score[i][j-1] + gap
            score[i][j] = max(diag, up, left)
    
    align1 = []
    align2 = []
    i, j = n, m
    
    while i > 0 or j > 0:
        if i > 0 and j > 0 and score[i][j] == score[i-1][j-1] + (match if seq1[i-1] == seq2[j-1] else mismatch):
            align1.append(seq1[i-1])
            align2.append(seq2[j-1])
            i -= 1
            j -= 1
        elif i > 0 and score[i][j] == score[i-1][j] + gap:
            align1.append(seq1[i-1])
            align2.append('-')
            i -= 1
        else:
            align1.append('-')
            align2.append(seq2[j-1])
            j -= 1
    
    return ''.join(reversed(align1)), ''.join(reversed(align2)), score[n][m]


s1="ACCGTGAAGCCAATAC"
s2="AGCGTGCAGCCAATAC"

aligned1, aligned2, final_score = needleman_wunsch(s1, s2)

print("Alignment:")
print(aligned1)
print(aligned2)
print("Score:", final_score)