# Part III: 8-mer duplex stability of candidate 2-letter alphabets (Turner 2004 nearest-neighbour stacks, NNDB tables)
# and the equivalent Stage 4 constant T_m8 (temperature where k_off(8-mer) = k0 = 1e6/h). CPU only: python alphabet_tm.py
import itertools, math, statistics as st
R=1.98720e-3
# Turner 2004 WC stacks (NNDB tables), key = 5'XY3'/3'X'Y'5' top strand dinucleotide
dG={'AA':-0.9,'UU':-0.9,'AU':-1.1,'UA':-1.3,'CU':-2.1,'AG':-2.1,'CA':-2.1,'UG':-2.1,'GU':-2.2,'AC':-2.2,'GA':-2.4,'UC':-2.4,
    'CG':-2.4,'GG':-3.3,'CC':-3.3,'GC':-3.4}
dH={'AA':-6.8,'UU':-6.8,'AU':-9.4,'UA':-7.7,'CU':-10.5,'AG':-10.5,'CA':-10.4,'UG':-10.4,'GU':-11.4,'AC':-11.4,'GA':-12.4,'UC':-12.4,
    'CG':-10.6,'GG':-13.4,'CC':-13.4,'GC':-14.9}
INI=(4.09,3.61); AUEND=(0.45,3.72)
def duplex(s):
    g=INI[0]; h=INI[1]
    for i in range(len(s)-1):
        d=s[i:i+2]; g+=dG[d]; h+=dH[d]
    for e in (s[0],s[-1]):
        if e in 'AU': g+=AUEND[0]; h+=AUEND[1]
    S=(h-g)/310.15
    return g,h,S
def tm(h,S,CT):  # non-self-compl., equal strands, CT total
    return h/(S+R*math.log(CT/4))-273.15
def stack_only(s):  # Stage-4-like: (l-1) stacks, no initiation/end terms -> T where dG=0
    h=sum(dH[s[i:i+2]] for i in range(len(s)-1)); g=sum(dG[s[i:i+2]] for i in range(len(s)-1))
    S=(h-g)/310.15; return h/S-273.15, h/(len(s)-1)
for name,al in [('G/C','GC'),('A/U','AU'),('4-letter','ACGU')]:
    rows=[]
    for s in itertools.product(al,repeat=8):
        s=''.join(s); g,h,S=duplex(s); t0,hbp=stack_only(s)
        rows.append((g,tm(h,S,1e-6),tm(h,S,1e-3),t0,hbp*4.184))
    def q(i):
        v=sorted(r[i] for r in rows); return f"{st.mean(v):7.1f} [{v[len(v)//20]:6.1f},{v[-len(v)//20-1]:6.1f}]"
    print(f"{name:9s} dG37 {q(0)}  Tm@1uM {q(1)}  Tm@1mM {q(2)}  stack-only T(dG=0) {q(3)}  dH/pair kJ {q(4)}")
print("effective Stage-4 T_m8: T where real k_off(8-mer) = k0 = 1e6/h, k_on 1e6..1e7 /M/s")
for name,al in [('G/C','GC'),('A/U','AU'),('4-letter','ACGU')]:
    out=[]
    for kon in (3.6e9,3.6e10):
        v=[]
        for s in itertools.product(al,repeat=8):
            s=''.join(s); g,h,S=duplex(s)
            # solve h - T S = R T ln(1e6/kon)  -> T = h/(S + R ln(1e6/kon))
            v.append(h/(S+R*math.log(1e6/kon))-273.15)
        v.sort(); out.append(f"kon={kon:.1e}/M/h: mean {st.mean(v):6.1f} 5-95% [{v[len(v)//20]:6.1f},{v[-len(v)//20-1]:6.1f}]")
    print(name, ' | '.join(out))
