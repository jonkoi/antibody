"""Combine the shipped VHH framework with crystal HER2; retain loop placeholders."""
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];OUT=Path(__file__).resolve().parent

def atomloop(path):
 lines=path.read_text().splitlines();cols=[l.strip().split('.',1)[1] for l in lines if l.startswith('_atom_site.')]
 rows=[dict(zip(cols,l.split())) for l in lines if l.startswith(('ATOM ','HETATM '))]
 return cols,rows
cols,s=atomloop(OUT/'pmvk_caplacizumab.cif')
_,t=atomloop(ROOT/'inputs/12-her2/her2_fv_1n8z.cif')
rows=[r for r in s if r['label_asym_id']=='H']
rows += [r for r in t if r['label_asym_id'] in ['C','D']]
# Retain canonical framework and all three original variable-length placeholders.
assert len({r['label_seq_id'] for r in rows if r['label_asym_id']=='H'})==86
out=['data_vhh_her2_hail_mary','#','loop_','_entity.id','_entity.type','_entity.pdbx_description','1 polymer VHH','2 polymer HER2_C','3 polymer HER2_D','#','loop_','_struct_asym.id','_struct_asym.entity_id','H 1','C 2','D 3','#','loop_']
out += ['_atom_site.'+c for c in cols]
for i,r in enumerate(rows,1):
 ch=r['label_asym_id'];r['label_entity_id']={'H':'1','C':'2','D':'3'}[ch];r['id']=str(i);r['atom_array_index']=str(i-1)
 if ch!='H':
  r.update(annotation_expsegmin_1_residue='-1',annotation_expsegmax_1_residue='-1',condition_sequence_1_residue=r['label_comp_id'],condition_hotspot_1_atom='False',**{'condition_cath-label_1_residue':'.','condition_CoordinateGroup_1_atom':'2'})
 out.append(' '.join(r.get(c,'.') for c in cols))
# Preserve component definitions, including expandable loop pseudo-residues.
src=(OUT/'pmvk_caplacizumab.cif').read_text();marker=src.index('loop_\n_chem_comp.id')
out += ['#',src[marker:]]
(OUT/'input.cif').write_text('\n'.join(out)+'\n')
(OUT/'design.yaml').write_text('''vhh_her2_hail_mary:
  input: input.cif
  contig: "C489-580, /0, D591-607, /0, H1-22, 13-13, H24-37, 9-10, H39-75, 12-21, H77-86"
  conditions:
    C_CRDGRP:
      1:
        - select: "sel('H') and not has_nan_coord()"
      2:
        - select: "sel('C/*/*/[N,CA,C,O]') and not has_nan_coord()"
        - select: "sel('D/*/*/[N,CA,C,O]') and not has_nan_coord()"
    C_SEQ:
      true:
        - select: "sel('H') and not has_nan_coord()"
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
''')
print('Wrote input with',len(rows),'atoms')
