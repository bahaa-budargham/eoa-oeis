import numpy as np, os, json

def load(p):
    M=[]
    for line in open(p):
        line=line.strip().strip(',').replace('"','')
        if line: M.append([int(float(x)) for x in line.split(',') if x.strip()!=''])
    return np.array(M,dtype=object)

# Automatically find the matrices whether they are in a folder or uploaded directly
mats_dir = 'EOA_incidence_matrices_csv'
if not os.path.isdir(mats_dir):
    mats_dir = '.'

mats={f.split('_matrix')[0].replace('EOA_',''):load(os.path.join(mats_dir, f)) 
      for f in sorted(os.listdir(mats_dir)) if f.endswith('_matrix.csv')}

K=90
def rowseq(M,i,K=K):
    n=len(M); v=np.zeros(n,dtype=object); v[i]=1; out=[]
    for _ in range(K): out.append(int(v.sum())); v=v@M
    return out

# --- target sequences (exact, K terms each) ---
Kt=92
KK=90
F=[1,1]
for _ in range(200): F.append(F[-1]+F[-2])
L=[2,1]
for _ in range(200): L.append(L[-1]+L[-2])
T={}
T['A000079']=[2**k for k in range(Kt)]                      # 1,2,4,8...
T['A001519']=[F[2*k-1] for k in range(1,Kt+1)]              # 1,2,5,13,34
T['A001906']=[F[2*k] for k in range(1,Kt+1)]                # 1,3,8,21,55
T['A002878']=[L[2*k+1] for k in range(Kt)]                  # 1,4,11,29,76 (odd-index Lucas)
T['A000225']=[2**k-1 for k in range(1,Kt+1)]                # 1,3,7,15,31
T['A000027']=[k+1 for k in range(Kt)]                       # 1,2,3,4
T['A005408']=[2*k+1 for k in range(Kt)]                     # 1,3,5,7
T['A000012']=[1]*Kt
pell=[0,1]
for _ in range(200): pell.append(2*pell[-1]+pell[-2])
T['A094723']=[pell[k+2]-2**k for k in range(Kt)]            # Pell(n+2)-2^n
T['A048739']=[]
acc=0
for k in range(Kt): acc+=pell[k]; T['A048739'].append(acc)  # partial sums of Pell (from 0)
T['A007070']=[0,1]
for k in range(2,Kt): T['A007070'].append(4*T['A007070'][-1]-2*T['A007070'][-2])
T['A000125']=[(k**3+5*k+6)//6 for k in range(Kt)]           # cake numbers 1,2,4,8,15,26
T['A000295']=[2**k-k-1 for k in range(Kt)]                  # 0,0,1,4,11,26

claims=[
 ('arabic','A094723','WOW1 alef'), ('arabic','A048739','WOW1 lam'),
 ('geez','A001906','WOW2'), ('geez','A001519','WOW2b'),
 ('thai','A001906','WOW2 thai'), ('thai','A001519','WOW2 thai-b'),
 ('thai','A002878','WOW2 Lucas ณ,ฟ'),
 ('hebrew','A000079','WOW3 hebrew'),
 ('japanese-latin','A000225','WOW3 mersenne'),
 ('hebrew','A000012','WOW4 he'),
 ('syriac','A000027','WOW4 syriac'),
 ('aramaic','A005408','WOW4 aramaic odd'),
 ('burmese','A007070','SW4 burmese 2+sqrt2'),
 ('armenian','A000125','SW1 cake'),
 ('syriac','A000295','SW1 eulerian/syriac'),
 ('aramaic','A000225','SW1 mersenne aramaic'),
 ('sanskrit','A000225','SW1 mersenne sanskrit'),
 ('kana','A000027','WOW4 kana shifted'),
 ('khmer','A000027','WOW4 khmer shifted')
]

res={}
for script,aid,tag in claims:
    if script not in mats:
        print(f"Warning: {script} matrix not found.")
        continue
    M=mats[script]; tgt=T[aid]; rows_hit=[]
    for i in range(len(M)):
        s=rowseq(M,i,K)
        if all(s[k]==tgt[k] for k in range(K)): rows_hit.append((i,'exact0'))
        elif aid in ('A001477','A000027') and all(s[k]==tgt[k]+1 for k in range(K)): rows_hit.append((i,'shift+1'))
        elif all(s[k]==tgt[k+1] for k in range(K-1)): rows_hit.append((i,'align1'))
        elif all(s[k+1]==tgt[k] for k in range(K-1)) and s[0]==1: rows_hit.append((i,'align-1'))
        elif all(s[k]==tgt[k+2] for k in range(K-2)): rows_hit.append((i,'align2'))
        elif all(s[k+2]==tgt[k] for k in range(K-2)) and s[0]==1: rows_hit.append((i,'align-2'))
    res[f'{script}:{aid} ({tag})']={'rows_matched':rows_hit,'n_rows':len(M)}

json.dump(res,open('verified_wow_output.json','w'),indent=1,default=str)
print("SUCCESS! verified_wow_output.json has been updated.")
for k,v in res.items():
    print(('✅' if v['rows_matched'] else '❌'),k,v['rows_matched'][:6])