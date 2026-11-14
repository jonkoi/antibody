"""Complete canonical-disulfide VHH scaffold + crystal HER2, all CDRs redesigned."""
from pathlib import Path
import sys
ROOT=Path(__file__).resolve().parents[2];OUT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'pipeline/denovo'))
from epitope import read_atoms,chain_seq
v=read_atoms(str(ROOT/'runs/denovo-vhh/boltz2_unguided/cont_full_s1_r15.pdb'))
t=read_atoms(str(ROOT/'inputs/12-her2/her2_fv_1n8z.cif'))
rows=[a for a in v if a[0]=='H']+[a for a in t if a[0] in ['C','D']]
s=chain_seq(rows,'H');assert len(s)==125 and s[21]=='C' and s[95]=='C' and s[114]=='W',s
cols='group_PDB id type_symbol label_atom_id label_alt_id label_comp_id label_asym_id label_entity_id label_seq_id pdbx_PDB_ins_code Cartn_x Cartn_y Cartn_z occupancy B_iso_or_equiv auth_seq_id auth_comp_id auth_asym_id auth_atom_id pdbx_PDB_model_num'.split()
out=['data_vhh_her2','#','loop_','_entity.id','_entity.type','_entity.pdbx_description','1 polymer VHH','2 polymer HER2_C','3 polymer HER2_D','#','loop_','_struct_asym.id','_struct_asym.entity_id','H 1','C 2','D 3','#','loop_','_entity_poly_seq.entity_id','_entity_poly_seq.num','_entity_poly_seq.mon_id','_entity_poly_seq.hetero']
seen=set()
for ch,ri,rn,an,x,y,z in rows:
 if (ch,ri) not in seen:out.append(f"{ {'H':1,'C':2,'D':3}[ch]} {ri} {rn} n");seen.add((ch,ri))
out+=['#','loop_']+['_atom_site.'+c for c in cols]
for i,(ch,ri,rn,an,x,y,z) in enumerate(rows,1):
 out.append(f"ATOM {i} {an[0]} {an} . {rn} {ch} { {'H':1,'C':2,'D':3}[ch]} {ri} ? {x:.3f} {y:.3f} {z:.3f} 1.0 0.0 {ri} {rn} {ch} {an} 1")
out+=['#'];(OUT/'input.cif').write_text('\n'.join(out)+'\n')
spec='''vhh_her2_hail_mary:
  input: input.cif
  contig: "C489-580, /0, D591-607, /0, H1-22, 13-13, H36-49, 9-10, H60-96, 12-21, H115-125"
  conditions:
    C_CRDGRP:
      1:
        - select: "sel('H/*/[1-22,36-49,60-96,115-125]')"
      2:
        - select: "sel('C/*/*/[N,CA,C,O]')"
        - select: "sel('D/*/*/[N,CA,C,O]')"
    C_SEQ:
      true:
        - select: "sel('H/*/[1-22,36-49,60-96,115-125]')"
        - select: "sel('C')"
        - select: "sel('D')"
    C_HOT:
      true:
        - select: "sel('D/*/593/NZ')"
        - select: "sel('C/*/573/CD1')"
        - select: "sel('D/*/591/CG2')"
        - select: "sel('C/*/572/CD')"
    C_CATH:
      "2.60.40":
        - select: "sel('H')"
'''
(OUT/'design.yaml').write_text(spec)
print('Complete VHH',len(s),'aa; 2 cysteines retained; atoms',len(rows))
