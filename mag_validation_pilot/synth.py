"""Minimal in-silico binning pilot: degrade a complete assembly into synthetic MAGs
by dropping whole contigs, optionally splicing in contigs from a donor genome."""
import os, random, sys

def read_fasta(p):
    recs=[]; h=None; buf=[]
    for line in open(p):
        if line.startswith('>'):
            if h is not None: recs.append((h,''.join(buf)))
            h=line[1:].strip(); buf=[]
        else: buf.append(line.strip())
    if h is not None: recs.append((h,''.join(buf)))
    return recs

def write_fasta(p, recs):
    with open(p,'w') as fh:
        for h,s in recs:
            fh.write('>'+h+'\n')
            for i in range(0,len(s),80): fh.write(s[i:i+80]+'\n')

def subset_to_fraction(recs, frac, rng):
    """Drop whole contigs until <= frac of the total bases remain (contig-order shuffled)."""
    total=sum(len(s) for _,s in recs); target=frac*total
    order=recs[:]; rng.shuffle(order)
    out=[]; acc=0
    for r in order:
        if acc+len(r[1])>target and out: continue
        out.append(r); acc+=len(r[1])
    return out, acc, total

HOME=os.path.expanduser("~")
src=f"{HOME}/Documents/sludge_genomes/1219030.9.fna"
don=f"{HOME}/Documents/sludge_genomes/1352.8344.fna"
S=read_fasta(src); D=read_fasta(don)
stot=sum(len(s) for _,s in S)
print(f"source: {len(S)} contigs, {stot/1e6:.2f} Mb")

os.makedirs("mags", exist_ok=True)
rows=[]
for comp in (1.0, 0.9, 0.7, 0.5):
    for contam in (0.0, 0.10):
        rng=random.Random(1729)
        sub,acc,_=subset_to_fraction(S, comp, rng)
        recs=[(f"src_{i}", s) for i,(h,s) in enumerate(sub)]
        cbases=0
        if contam>0:
            dsub,cbases,_=subset_to_fraction(D, contam*acc/sum(len(s) for _,s in D), rng)
            recs += [(f"contam_{i}", s) for i,(h,s) in enumerate(dsub)]
        name=f"mag_c{int(comp*100)}_x{int(contam*100)}"
        write_fasta(f"mags/{name}.fna", recs)
        rows.append((name,len(recs),acc,cbases,acc/stot,cbases/(acc+cbases) if acc+cbases else 0))
print(f"{'mag':22} {'contigs':>8} {'src_bp':>10} {'cont_bp':>9} {'compl':>7} {'contam':>7}")
for n,nc,a,c,cp,ct in rows: print(f"{n:22} {nc:8d} {a:10d} {c:9d} {cp:7.3f} {ct:7.3f}")
