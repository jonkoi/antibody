"""Preserve remote 24-design order; append unique validated backups; retain provenance."""
from pathlib import Path
import csv,json,re,subprocess
ROOT=Path(__file__).resolve().parents[2];OUT=Path(__file__).resolve().parent

def read(p):return list(csv.DictReader((ROOT/p).open()))
def write(p,rows,fields=None):
 with (ROOT/p).open('w') as f:
  w=csv.DictWriter(f,fieldnames=fields or list(rows[0]),lineterminator='\n',extrasaction='ignore');w.writeheader();w.writerows(rows)
def key(r):return (r.get('kind',r.get('format')),r['heavy_sequence'],r.get('light_sequence') or '')
remote=read('runs/reconciliation/team15_01.remote.csv');current=read('runs/all-ranked/valid_candidates_ranked_with_notes.csv');old={key(r):r for r in read('runs/reconciliation/local_ranked_before_sync.csv')};lookup={key(r):r for r in current}
assert len(remote)==24 and len(current)==546
selected={key(r) for r in remote};ordered=[lookup[key(r)] for r in remote]+[r for r in current if key(r) not in selected]
assert len({key(r) for r in ordered})==546
# Raw predictor output, not laboratory measurements.
ddg={r['variant_id']:r for r in read('runs/stab-ddg/coworker-result/9490750b-9332-5d42-b9be-05b4011c311f.csv')}
gen={r['variant_id']:r for r in read('runs/stab-ddg/designs-result/f0e3df20-70c2-5f9e-a25e-c1c46e59c494.csv')}
info={'WGHDGPYWMD':('cw_1',4.07),'WGGKGPYWMD':('cw_2',6.74),'WGGMGPYWMD':('cw_3',7.38),'WGGRGPYWMD':('cw_4',3.48),'WGPDGPYWMD':('cw_5',3.20)}
submission=[];notes=[];cross=[]
for i,r in enumerate(ordered,1):
 rr=remote[i-1] if i<=24 else {};prior=old.get(key(r),{});sid=rr.get('sequence_id') or f'r{i:03d}_'+r['sequence_id']
 window=r['heavy_sequence'][98:108] if r['format']=='paired' else ''
 datum={};distance='';parent=''
 if window in info:
  ident,distance=info[window];datum=ddg[ident]
  if r['framework_edits']:parent='CDRH3-only parent; framework combination not scored';distance=''
 elif window=='SLGSADYAFD':datum=gen['her2_cdrh3_61'];distance=2.51
 rawreason=rr.get('design_description','')
 reason=rawreason or r['reason_picked']
 for a,b in [('tier 1: measured','tier 1: computationally characterised'),('measured binding ddG','predicted binding ddG'),('measured ddG','predicted ddG'),('Binding ddG','Predicted binding ddG'),('statistically indistinguishable from wild-type','model estimate near zero; experimental equivalence untested'),('so a glycan cannot reach either interface','protein-only distance does not establish glycan safety'),('benefit is glycosylation-dependent and unconfirmed','glycosylation dependence and transfer remain unresolved'),('HER2 contacts kept','crystal-contact residue identities retained')]:reason=reason.replace(a,b)
 desc=f'Rank {i}; {r["track"]}. '+reason+' Latest independent contact evidence: '+r['epitope_evidence']+' '+r['risk_notes']
 submission.append(dict(sequence_id=sid,kind=r['format'],heavy_sequence=r['heavy_sequence'],light_sequence=r['light_sequence'],design_description=desc))
 note=dict(r);note.update(submission_rank=i,submission_id=sid,colleague_rank=i if i<=24 else '',local_previous_rank=prior.get('rank',''),local_previous_number=prior.get('candidate_number',''),selection_background=reason,predicted_binding_ddg_kcal_mol=datum.get('binding_ddg_kcal_mol',''),ddg_model_mc_standard_error=datum.get('mc_standard_error_kcal_mol',''),ddg_construct_caveat=parent,reported_template_O103_NZ593_A=distance,experimental_measurements_on_design='none',historical_probability_threshold_nM=0.146 if r.get('model_p_below_0_146_nM') else '',actual_wt_kd_median_nM=0.1345)
 notes.append(note)
 if i<=24:cross.append({k:note[k] for k in ['submission_rank','submission_id','sequence_id','local_previous_rank','local_previous_number','selection_background','predicted_binding_ddg_kcal_mol','ddg_model_mc_standard_error','ddg_construct_caveat','reported_template_O103_NZ593_A','epitope_evidence','risk_notes','heavy_sequence','light_sequence']})
assert [key(r) for r in submission[:24]]==[key(r) for r in remote]
accepted=json.loads((ROOT/'runs/all-ranked/accepted.json').read_text());assert {key(r) for r in submission}=={key(r) for r in accepted}
assert len({r['sequence_id'] for r in submission})==546
write('team15_01.csv',submission);write('SUBMISSION.csv',submission)
write('runs/reconciliation/ranked_with_evidence.csv',notes);write('runs/reconciliation/priority24_crosswalk.csv',cross)
write('runs/reconciliation/de_novo_ranked_with_evidence.csv',[r for r in notes if r['track']=='de_novo'])
summary=dict(remote_commit='ab8ade5',original_priority_count=24,matching_previous_candidates=10,new_shortlisted_candidates=14,union_count=653,screen_pass=546,screen_rejected=107,optimisation=369,de_novo=177,first24_sequences_and_order_preserved=True,all_exported_genotypes_match_fresh_screen=True,ranking_policy='Original 24 in unchanged sequence/order; remaining validated candidates appended in existing heuristic order; track-specific ranks retained.',no_design_has_new_experimental_measurements=True)
(OUT/'sync_summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
