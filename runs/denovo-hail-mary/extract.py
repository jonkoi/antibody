"""Extract candidates from RFD4 outputs; never infer binding from design confidence."""
from pathlib import Path
import sys,csv,json
ROOT=Path(__file__).resolve().parents[2];OUT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'pipeline/denovo'))
from epitope import read_atoms,chain_seq
rows=[];audit=[]
for p in sorted((OUT/'design_complete').rglob('*_generated.cif')):
 atoms=read_atoms(str(p));seqs={ch:chain_seq(atoms,ch) for ch in {a[0] for a in atoms}}
 binder=[(ch,s) for ch,s in seqs.items() if len(s)>100 and 'EVQL' in s[:10]]
 assert len(binder)==1,(p,seqs)
 ch,s=binder[0];sid='dn_'+p.parent.name
 rows.append(dict(sequence_id=sid,kind='vhh',heavy_sequence=s,light_sequence='',design_description='One-shot RFD4 antibody-checkpoint HER2 VHH; all 3 CDRs redesigned, canonical-disulfide framework retained; four crystal hotspot residues conditioned'))
 audit.append(dict(sequence_id=sid,structure=str(p.relative_to(OUT)),chain=ch,length=len(s),cysteines=s.count('C'),sequences=seqs))
assert len(rows)==12,len(rows)
with (OUT/'candidates_12.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
(OUT/'sequence_audit.json').write_text(json.dumps(audit,indent=2)+'\n')
(OUT/'candidates.fa').write_text(''.join('>'+r['sequence_id']+'\n'+r['heavy_sequence']+'\n' for r in rows))
print('Extracted',len(rows),'designs, unique',len({r['heavy_sequence'] for r in rows}))
print('Cysteine counts:',[r['cysteines'] for r in audit])
