#!/usr/bin/env python3
"""
EOA OEIS baseline experiment -- ready-to-run local script.
Compares OEIS prefix-match hit rates of the real EOA letter-name sequences
against random-matrix controls of matched size/density/degrees.

Run:  pip install numpy
      python3 baseline_oeis_lookup.py --mats_dir mats/EOA_incidence_matrices_csv --out results
Requires stripped.gz from https://oeis.org/stripped.gz (download by hand; the site
blocks datacenter IPs). Place it next to this script or pass --stripped.
"""
import numpy as np, gzip, os, json, argparse, collections, math

ap=argparse.ArgumentParser()
ap.add_argument('--mats_dir',default='.')
ap.add_argument('--stripped',default='stripped.gz')
ap.add_argument('--out',default='results')
ap.add_argument('--samples',type=int,default=100)
ap.add_argument('--windows',default='10,15,20,30,40')
a=ap.parse_args(); os.makedirs(a.out,exist_ok=True)
W=[int(w) for w in a.windows.split(',')]; WMAX=max(W)

def load(p):
    M=[]
    for line in open(p):
        line=line.strip().strip(',').replace('"','')
        if line: M.append([int(float(x)) for x in line.split(',') if x.strip()!=''])
    return np.array(M,dtype=np.int64)
mats={f.split('_matrix')[0].replace('EOA_',''):load(os.path.join(a.mats_dir,f))
      for f in sorted(os.listdir(a.mats_dir)) if f.endswith('_matrix.csv')}

def seqs(M,K):
    n=len(M); out=[]
    for i in range(n):
        v=np.zeros(n,dtype=np.int64); v[i]=1; s=[]
        for _ in range(K): s.append(int(v.sum())); v=v@M
        out.append(tuple(s))
    return out

# ---- index OEIS entry prefixes
idx={w:collections.defaultdict(list) for w in W}
n_in=0
with gzip.open(a.stripped,'rt',encoding='utf-8',errors='replace') as f:
    for line in f:
        if not line.startswith('A'): continue
        parts=line.split(); aid=parts[0]; terms=[]
        ok=True
        for t in parts[1:]:
            t=t.strip(',')
            try: terms.append(int(t))
            except: ok=False; break
            if len(terms)==WMAX: break
        if ok and len(terms)>=WMAX:
            for w in W: idx[w][tuple(terms[:w])].append(aid)
            n_in+=1
print(f'indexed {n_in} OEIS entries')

def wilson(k,n):
    if n==0: return (0,0)
    z=1.96; p=k/n
    h=1+z*z/n
    c=(p+z*z/(2*n))/h; r=z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/h
    return (round(c-r,4),round(c+r,4))

def hit_stats(groups_all,w):
    G=set(groups_all); m=sum(1 for g in G if g in idx[w])
    return len(G),m,m/max(len(G),1)

rows=[]; summary={}
for name,M in mats.items():
    n=len(M); dens=float((M>0).mean())
    rng=np.random.default_rng(abs(hash(name))%99991)
    # real
    rs=seqs(M,WMAX)
    for w in W:
        G,m,r=hit_stats([g[:w] for g in rs],w)
        rows.append(dict(script=name,kind='real',window=w,groups=G,matched=m,rate=round(r,4)))
    # nulls
    for kind in ('L1_density','L2_rowshuffle','L3_degreeswap'):
        for _ in range(a.samples):
            if kind=='L1_density':
                Mn=(rng.random((n,n))<dens).astype(np.int64)
            elif kind=='L2_rowshuffle':
                Mn=np.zeros((n,n),dtype=np.int64)
                for r0 in range(n):
                    row=M[r0][:].tolist(); rng.shuffle(row); Mn[r0]=row
            else:
                Mn=M.copy()
                for _ in range(2*n):
                    i,j,k,l=rng.integers(0,n,4)
                    if i!=k and Mn[i,j]>0 and Mn[k,l]>0 and Mn[i,l]==0 and Mn[k,j]==0:
                        Mn[i,l],Mn[k,j]=Mn[i,j],Mn[k,l]; Mn[i,j]=0; Mn[k,l]=0
            ss=seqs(Mn,WMAX)
            for w in W:
                G,m,r=hit_stats([g[:w] for g in ss],w)
                rows.append(dict(script=name,kind=kind,window=w,groups=G,matched=m,rate=round(r,4)))
allr=collections.defaultdict(list)
for r0 in rows: allr[(r0['kind'],r0['window'])].append(r0['rate'])
for k,v in sorted(allr.items()):
    lo,hi=wilson(sum(int(x*1000) for x in v),1000*len(v))
    summary[f'{k[0]}_w{k[1]}']={'mean_rate':round(float(np.mean(v)),4),'sd':round(float(np.std(v)),4),'ci':wilson(round(float(np.mean(v))*1000),1000),'n':len(v)}
json.dump({'summary':summary,'rows':rows},open(os.path.join(a.out,'baseline.json'),'w'),indent=1)
print(json.dumps(summary,indent=1))
print('NOTE: interpret real vs controls within each window; report both group-level rates.')
