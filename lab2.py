code = {
    "UUU":"Phe","UUC":"Phe","UUA":"Leu","UUG":"Leu","UCU":"Ser","UCC":"Ser","UCA":"Ser",
    "UCG":"Ser","UAU":"Tyr","UAC":"Tyr","UGU":"Cys","UGC":"Cys","CUU":"Leu","CUC":"Leu",
    "CUA":"Leu","CUG":"Leu","CCU":"Pro","CCA":"Pro","CCC":"Pro","CCG":"Pro","CAU":"His",
    "CAC":"His","CAA":"Gln","CAG":"Gln","CGU":"Arg","CGC":"Arg","CGA":"Arg","CGG":"Arg",
    "AUU":"Ile","AUC":"Ile","AUA":"Ile","AUG":"Met","ACU":"Thr","ACC":"Thr","ACA":"Thr",
    "ACG":"Thr","AAU":"Asn","AAC":"Asn","AAA":"Lys","AAG":"Lys","AGU":"Ser","AGC":"Ser",
    "AGA":"Arg","AGG":"Arg","GUU":"Val","GUC":"Val","GUA":"Val","GUG":"Val","GCU":"Ala",
    "GCC":"Ala","GCA":"Ala","GCG":"Ala","GAU":"Asp","GAC":"Asp","GAA":"Glu","GAG":"Glu",
    "GGU":"Gly","GGC":"Gly","GGA":"Gly","GGG":"Gly","UAA":"Stop","UAG":"Stop","UGA":"Stop"
}

def translation(secv):
    x=[]
    for i in range(0,len(secv)-2,3):
        a=secv[i:i+3]
        b=code[a]
        x.append(b)
        if b=="Stop":
            break
    return x

def translationv2(secv):
    x=[]
    started=False
    for i in range(0,len(secv)-2,3):
        a=secv[i:i+3]
        b=code[a]
        if b=="Stop":
            if started==False:
                started=True
            else:
                break
        elif started:
            x.append(b)
    return x
            
            

s="ACGAUGCAGUACUGAAUCGAUCGAUCGACGAUAGCAGCAUAGUCGAUAUC"
s1=translation(s)
for i in range(0,len(s1)):
    print(s1[i])
    
print("Exercise 2\n")
    
s2=translationv2(s)
for i in range(0,len(s2)):
    print(s2[i])
#trebuie translatat
