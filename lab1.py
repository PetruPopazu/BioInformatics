import matplotlib.pyplot as plt

def frequencies(dna):
    freq2 = {}
    for i in range(len(dna) - 1):
        comb = dna[i:i+2]
        if comb in freq2:
            freq2[comb] += 1
        else:
            freq2[comb] = 1
            
    freq3 = {}
    for i in range(len(dna) - 2):
        comb = dna[i:i+3]
        if comb in freq3:
            freq3[comb] += 1
        else:
            freq3[comb] = 1
            
    return freq2, freq3

s = "adcdacdgagdacdgadacccccccccccccc"
cnta = 0
cntc = 0
cntd = 0
cntg = 0
   
for i in range(len(s)):
    if s[i] == 'a':
        cnta += 1
    if s[i] == 'd':
        cntd += 1
    if s[i] == 'c':
        cntc += 1
    if s[i] == 'g':
        cntg += 1
            
print(cnta)
print(cntd)
print(cntc)
print(cntg)

print(frequencies(s))
freq2,freq3=frequencies(s)
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
ax1.bar(freq2.keys(), freq2.values(), color='skyblue', edgecolor='black')
ax1.set_title('comb of 2 charac')
ax2.bar(freq3.keys(), freq3.values(), color='red', edgecolor='black')
ax2.set_title('comb of 3 charac')

plt.show()