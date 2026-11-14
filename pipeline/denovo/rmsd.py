import sys, numpy as np
sys.path.insert(0,'/tmp/claude-10009/-home-duy-khoi-tran-esmartsystems-com-antibody/e98dac37-fc3f-4bee-9602-26339af2e4c1/scratchpad')
from epitope import read_atoms
def kabsch(P,Q):
    P=np.array(P,float); Q=np.array(Q,float)
    P=P-P.mean(0); Q=Q-Q.mean(0)
    U,S,Vt=np.linalg.svd(P.T@Q)
    d=np.sign(np.linalg.det(Vt.T@U.T))
    R=Vt.T@np.diag([1,1,d])@U.T
    return float(np.sqrt((((P@R.T)-Q)**2).sum(1).mean()))
ref=read_atoms('/home/duy.khoi.tran@esmartsystems.com/antibody/inputs/12-her2/her2_fv_1n8z.cif')
R={a[1]:a[4:7] for a in ref if a[3]=='CA' and a[0] in ('C','D')}
for f in sys.argv[1:]:
    des=read_atoms(f)
    D={a[1]+488:a[4:7] for a in des if a[0]=='A' and a[3]=='CA'}
    if not D: D={a[1]+488:a[4:7] for a in des if a[0]=='C' and a[3]=='CA'}
    c=sorted(set(R)&set(D))
    # epitope-local rmsd too
    epi=[557,558,560,561,569,570,571,572,573,579,591,593,596,598,600,602,603,604,605]
    ce=[x for x in c if x in epi]
    print(f.split('/')[-2]+'/'+f.split('/')[-1], 'n=%d'%len(c),
          'allRMSD=%.2f'%kabsch([D[x] for x in c],[R[x] for x in c]),
          'epiRMSD=%.2f'%kabsch([D[x] for x in ce],[R[x] for x in ce]) if len(ce)>3 else '')
