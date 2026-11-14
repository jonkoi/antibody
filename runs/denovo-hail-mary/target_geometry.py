"""Check both HER2 fragments together against crystal geometry (one rigid fit)."""
from pathlib import Path
import sys,csv
import numpy as np
ROOT=Path(__file__).resolve().parents[2];OUT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'pipeline/denovo'))
from epitope import read_atoms,chain_seq
ref=read_atoms(str(ROOT/'inputs/12-her2/her2_fv_1n8z.cif'))
c=chain_seq(ref,'C');d=chain_seq(ref,'D')
def coords(a,ch):return np.array([v[4:7] for v in sorted(a,key=lambda v:v[1]) if v[0]==ch and v[3]=='CA'])
y=np.concatenate([coords(ref,'C'),coords(ref,'D')]);y-=y.mean(0)
rows=[]
for folder in ['unguided','confirmation']:
 for p in sorted((OUT/folder).glob('*.pdb')):
  a=read_atoms(str(p));seqs={ch:chain_seq(a,ch) for ch in {v[0] for v in a}}
  ch1=next(ch for ch,s in seqs.items() if s==c);ch2=next(ch for ch,s in seqs.items() if s==d)
  x=np.concatenate([coords(a,ch1),coords(a,ch2)]);x-=x.mean(0)
  u,s,v=np.linalg.svd(x.T@y);m=np.eye(3);m[-1,-1]=np.linalg.det(u@v);rot=u@m@v
  rmsd=float(np.sqrt(np.mean(np.sum((x@rot-y)**2,axis=1))))
  rows.append(dict(stage=folder,structure=p.name,target_rmsd_A=round(rmsd,3)))
with (OUT/'target_geometry.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
print([r for r in rows if 'mary_1.' in r['structure'] or 'mary_1_model' in r['structure'] or 'ctl_wt' in r['structure']])
