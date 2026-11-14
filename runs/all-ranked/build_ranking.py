#!/usr/bin/env python3
"""Consolidate actual generated/shortlisted HER2 designs; exclude enumeration pools."""
import csv,json,hashlib,sys,re,collections
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];OUT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/'pipeline'))
from build_submission_table import WT_HEAVY,WT_LIGHT,graft,COLUMNS
sys.path.insert(0,str(ROOT/'pipeline/denovo'))
from epitope import read_atoms,chain_seq

def table(path):
 p=ROOT/path
 return list(csv.DictReader(p.open())) if p.exists() else []
def js(path):
 p=ROOT/path
 return json.loads(p.read_text()) if p.exists() else []

pool={};controls={};aliases={}
def key(kind,h,l=''):return (kind,h,l)
def add(r,source,background):
 k=r.get('kind','paired');h=r.get('heavy_sequence','').strip();l=r.get('light_sequence','').strip()
 if not h:return
 q=key(k,h,l)
 control=k=='paired' and len(h)==120 and h[98:108]=='WGGDGFYAMD'
 dest=controls if control else pool
 if q not in dest:
  digest=hashlib.sha256(('|'.join(q)).encode()).hexdigest()[:12]
  native=k=='paired' and len(h)==120 and l==WT_LIGHT and h==graft(h[98:108])
  fw=k=='paired' and len(h)==120 and l==WT_LIGHT and h[84]=='N' and h[86]=='T' and h[:84]+WT_HEAVY[84]+h[85:86]+WT_HEAVY[86]+h[87:]==graft(h[98:108])
  sid=('opt_'+h[98:108]) if native else ('opt_fw_'+h[98:108] if fw else 'vhh_'+digest)
  dest[q]=dict(sequence_id=sid,kind=k,heavy_sequence=h,light_sequence=l,design_description=r.get('design_description','HER2 antibody design'),native=native,framework_variant=fw,sources=set(),source_ids=set(),background=set())
 d=dest[q];d['sources'].add(source);d['source_ids'].add(r.get('sequence_id',''));d['background'].add(background)
 if r.get('sequence_id'):aliases[(source,r['sequence_id'])]=q
 return q

sources=[
 ('runs/reconciliation/team15_01.remote.csv','team15-priority','Latest team priority list from remote ab8ade5; source descriptions contain historical terminology; presentation audit distinguishes predictions from measurements.'),
 ('runs/final-shortlist/submission_table_18.csv','final18','Final 18-model optimisation shortlist; includes later epitope failures retained for completeness.'),
 ('runs/framework-hail-mary/validation_panel.csv','framework','Published S85N/R87T stabilisation hypothesis combined with shortlisted CDRH3s; matched originals included.'),
 ('runs/denovo-hail-mary/candidates_12.csv','rfd4-vhh','One-shot antibody-specific RFD4 VHH campaign; all three CDRs redesigned with four explicit HER2 hotspot atoms.'),
 ('runs/denovo-vhh/w1_ss_abguard.csv','boltzgen-disulfide','Earlier BoltzGen VHH campaign with restored canonical disulfide; included as exploratory/long-shot pool.'),
 ('runs/team15-claude-n100-seed1/submission_all.csv','rfd4-seed1','Original 100-design HER2-finetuned RFD4 production batch.')]
for path,source,bg in sources:
 for r in table(path):
  if source=='team15-priority':
   r=dict(r);text=r['design_description']
   text=text.replace('tier 1: measured','tier 1: computationally characterised').replace('measured binding ddG','predicted binding ddG').replace('measured ddG','predicted ddG').replace('Binding ddG','Predicted binding ddG')
   text=text.replace('statistically indistinguishable from wild-type','model estimate close to zero; no experimental equivalence established')
   text=text.replace('so a glycan cannot reach either interface','protein-only distance does not establish glycan safety')
   text=text.replace('benefit is glycosylation-dependent and unconfirmed','glycosylation dependence and transfer are unresolved')
   text=text.replace('HER2 contacts kept','crystal-contact residue identities retained')
   r['design_description']=text
  add(r,source,bg)
for r in table('data/shortlist_optimisation.csv'):
 c=r['cdrh3'];add(dict(sequence_id='prior_'+c,kind='paired',heavy_sequence=graft(c),light_sequence=WT_LIGHT,design_description='Earlier winner-prior/expression shortlist '+c),'winner-prior55','Earlier shortlist based on measured SPR winner positions, FACS model and chemical-liability reasoning; later ML shortlist supersedes its top ranking.')
for r in table('data/af2ig_scores.csv'):
 if r['source']=='ml_ranked':
  c=r['cdrh3'];add(dict(sequence_id=r['design_id'],kind='paired',heavy_sequence=graft(c),light_sequence=WT_LIGHT,design_description='Evaluated ML top-100 CDRH3 '+c),'ml-top100','One of the 100 ML-ranked CDRH3s actually rescored by AF2-IG; retained beyond the final 18.')
for base,source,bg in [('runs/all-ranked/source/seed7/inference_samples','rfd4-seed7','Second 100-design HER2 RFD4 batch, seed 7; F104/Y105 sequence conserved in the input.'),('runs/her2-validation-seed1/designs','pilot4','Initial four-design pipeline validation batch.')]:
 if source=='pilot4' and (OUT/'source/pilot_candidates.csv').exists():
  for r in table('runs/all-ranked/source/pilot_candidates.csv'):add(r,source,bg)
  continue
 if source=='rfd4-seed7' and (OUT/'source/seed7_candidates.csv').exists():
  for r in table('runs/all-ranked/source/seed7_candidates.csv'):add(r,source,bg)
  continue
 for p in sorted((ROOT/base).rglob('*_generated.cif')):
  a=read_atoms(str(p));add(dict(sequence_id=p.parent.name,kind='paired',heavy_sequence=chain_seq(a,'B'),light_sequence=chain_seq(a,'A'),design_description=bg),source,bg)
for r in js('runs/denovo-vhh/rows_w1.json'):
 add(dict(sequence_id=r['_job']+'_r'+str(r['final_rank']),kind='vhh',heavy_sequence=r['designed_sequence'],light_sequence='',design_description='BoltzGen original-sequence HER2 VHH'),'boltzgen-original','Original BoltzGen VHH sequence; canonical framework cysteines were often absent; independent epitope targeting failed for the tested earlier campaign members.')
for j in js('runs/all-ranked/source/previous_jobs_summary.json'):
 for r in js(j['path']).get('designs',[]):
  add(dict(sequence_id=j['name']+'_r'+str(r['final_rank']),kind='vhh',heavy_sequence=r.get('designed_sequence',''),light_sequence='',design_description='Previously generated BoltzGen HER2 VHH '+j['name']),'boltzgen-jobs','Previously completed HER2 BoltzGen job '+j['name']+'; no new generation; independent contact evidence is absent unless stated.')
assert len(pool)==653,(len(pool),'Source inventory changed; review count')
assert len({d['sequence_id'] for d in pool.values()})==len(pool)

# Exact full-sequence evidence, never transferred across framework/format variants.
contacts={}
def contact(q,good,n,evidence,partial=False):
 if q in pool:contacts[q]=dict(good=good,n=n,evidence=evidence,partial=partial)
for r in table('runs/final-shortlist/epitope_check_18.csv'):
 if len(r['cdrh3'])==10:
  contact(key('paired',graft(r['cdrh3']),WT_LIGHT),int(r['verdict']=='PASS'),1,'Single unguided epitope check: '+r['verdict']+'; recall '+r['epi_recall']+'; core '+r['core4_n']+'/4.',r['verdict']=='CAUTION')
for r in table('runs/framework-hail-mary/epitope_summary.csv'):
 q=aliases.get(('framework',r['target']))
 contact(q,int(r['good']),int(r['n']),r['good']+'/'+r['n']+' strict epitope-contact samples in matched framework experiment; minimum recall '+r['min_recall']+'.')
# Actual prior independent VHH structures: these were already found off-site.
for p in (ROOT/'runs/denovo-vhh/boltz2_unguided').glob('*.pdb'):
 if p.name.startswith('ctl_'):continue
 a=read_atoms(str(p));seqs={ch:chain_seq(a,ch) for ch in {v[0] for v in a}}
 for ch,s in seqs.items():
  if s.startswith('EVQL'):contact(key('vhh',s,''),0,1,'Earlier independent unguided co-fold: off-epitope (0/4 core hotspots).')
for r in table('runs/denovo-hail-mary/epitope_summary.csv'):
 q=aliases.get(('rfd4-vhh',r['target']))
 if q:contact(q,int(r['vhh_overlap']),int(r['n']),r['vhh_overlap']+'/'+r['n']+' VHH-overlap samples (>=50% epitope recall + all 4 core residues); strict Fab-footprint passes '+r['strict']+'.')
for r in table('runs/denovo-hail-mary/confirmation_epitope_summary.csv'):
 q=aliases.get(('rfd4-vhh',r['target']))
 if q:
  prev=contacts[q];g=prev['good']+int(r['vhh_overlap']);n=prev['n']+int(r['n'])
  contact(q,g,n,f'{g}/{n} VHH-overlap samples across seeds 73/1073 (>=50% recall + all 4 core); 0/6 strict Fab-footprint passes.')

needed={d['heavy_sequence'][98:108] for d in pool.values() if d['native'] or d['framework_variant']}
ml={r['cdrh3']:r for r in table('data/legal_3mutants_scored_ml.csv') if r['cdrh3'] in needed}
prior={r['cdrh3']:r for r in table('data/legal_3mutants_ranked.csv') if r['cdrh3'] in needed}
af={}
for r in table('data/af2ig_scores.csv'):
 if r['cdrh3'] in needed and r['source'] in ['ml_ranked','production','ranked'] and r['status']=='ok':af[r['cdrh3']]=r['iptm']
tm={}
for name in ['runs/denovo-vhh/w1_tm.csv','runs/denovo-vhh/w1_ss_tm.csv','runs/denovo-hail-mary/nanomelt.csv']:
 for r in table(name):tm[r['Sequence']]=r['NanoMelt Tm (C)']
# Fast incumbent sequence safety score for generated paired candidates outside ML pool.
from rank_candidates import load_facs,fit,score
weights=fit(load_facs(ROOT/'data/her2-trastuzumab/tz_her2_facs_3classes_524346.csv.zip'),10)
accepted={r['sequence_id'] for r in js('runs/all-ranked/accepted.json')}
nov={r['sequence_id']:r for r in js('runs/all-ranked/novelty.json')}
liab={r['sequence_id']:r for r in js('runs/all-ranked/liabilities.json')}
final='--final' in sys.argv
if final:
 assert (OUT/'accepted.json').exists() and (OUT/'novelty.json').exists()
 assert accepted <= {d['sequence_id'] for d in pool.values()}

team_rows=table('runs/reconciliation/team15_01.remote.csv')
team_by_key={key(r['kind'],r['heavy_sequence'],r.get('light_sequence','')):(i,r) for i,r in enumerate(team_rows,1)}
previous={key(r['format'],r['heavy_sequence'],r.get('light_sequence','')):r for r in table('runs/reconciliation/local_ranked_before_sync.csv')}
ranked=[]
for q,d in pool.items():
 h=d['heavy_sequence'];l=d['light_sequence'];sid=d['sequence_id'];c=h[98:108] if d['native'] or d['framework_variant'] else ''
 m=ml.get(c,{}) if d['native'] else {};parent=ml.get(c,{}) if d['framework_variant'] else {}
 e=contacts.get(q,{});nr=nov.get(sid,{});track={'novel':'optimisation','de_novo':'de_novo'}.get(nr.get('novelty_class'),'pending_screen')
 status=('pass' if sid in accepted else 'rejected') if final else 'pending'
 risks=[]
 if h.count('C')==0:risks.append('Canonical framework disulfide absent; expression/stability risk.')
 if h.count('C')%2 or l.count('C')%2:risks.append('Odd cysteine count; inspect screen rejection.')
 sequons=[str(x.start()+1)+':'+x.group() for x in re.finditer(r'N[^P][ST]',h)]
 if sequons:risks.append('VH potential N-glycosylation sequon(s) '+','.join(sequons)+'; glycans not modelled.')
 if d['framework_variant']:risks.append('Framework thermal benefit on this combined construct unmeasured; glycosylation dependence unresolved.')
 if not d['native']:risks.append('No measured affinity; CDR-only affinity model not transferred to this construct.')
 if e and e['good']==0:risks.append('Predicted epitope concern: '+e['evidence'])
 if not e:risks.append('No independent epitope-contact validation for this exact sequence.')
 if tm.get(h):risks.append('NanoMelt apparent Tm prediction, not measured or a proven improvement over trastuzumab; CC BY-NC-SA 4.0.')
 ll=liab.get(sid,{})
 if ll.get('errors'):risks.append('Live screen errors: '+str(ll['errors']))
 p=float(m.get('p_beat_wt') or 0);safety=score(c,weights) if len(c)==10 else -999
 if e and e['good']==e['n'] and e['n']>=3:
  group=12 if d['framework_variant'] else 10;why='Prioritised for reproducible independent epitope contacts; '+('literature-based framework stability hypothesis.' if d['framework_variant'] else 'CDR affinity-model support where applicable.')
 elif e and e['good']>0 and e['n']>=3:
  group=15 if d['kind']=='vhh' else 30;why='Exploratory: partial reproducibility of independent epitope overlap; retained as a higher-risk option.'
 elif e and e['good']>0:
  group=20;why='Shortlisted optimisation design with one favourable independent epitope-contact prediction.'
 elif e:
  group=60 if e.get('partial') else 80;why='Long shot retained at user request despite predicted epitope caution/off-site docking.'
 elif d['native'] and m:
  group=40;why='Native-framework candidate supported by the existing sequence affinity model; independent epitope validation absent unless separately reported.'
 elif d['kind']=='paired':
  group=45;why='Generated or earlier-shortlisted optimisation backup; retained for sequence diversity and completeness.'
 else:
  group=55 if h.count('C')==0 else 50;why='Untested VHH long shot retained for diversity/completeness; NanoMelt support shown only when available.'
 if status=='rejected':group=95;why='Retained in audit only: rejected by live abguard; excluded from submission export.'
 row=dict(rank='',track_rank='',candidate_number='',sequence_id=sid,track=track,format=d['kind'],screen_status=status,reason_picked=why,background_notes=' | '.join(sorted(d['background'])),risk_notes=' | '.join(risks),cdrh3=nr.get('cdrh3') or c,heavy_sequence=h,light_sequence=l,heavy_length=len(h),light_length=len(l),framework_edits='S85N R87T (one-based VH sequence)' if d['framework_variant'] else '',epitope_evidence=e.get('evidence','untested'),epitope_good_samples=e.get('good',''),epitope_samples=e.get('n',''),nanomelt_tm_C=tm.get(h,''),model_p_below_0_146_nM=m.get('p_beat_wt',''),model_pred_kd_nM=m.get('pred_kd_nM',''),parent_cdr_model_p_below_0_146_nM=parent.get('p_beat_wt',''),facs_safety_logodds=round(safety,4) if len(c)==10 else '',af2ig_iptm=af.get(c,'') if d['native'] else '',sources=';'.join(sorted(d['sources'])),original_ids=';'.join(sorted(d['source_ids'])))
 priority,teamrow=team_by_key.get(q,('',{}))
 old=previous.get(q,{})
 row.update(model_probability_threshold_nM=0.146 if (m or parent) else '',actual_wt_median_nM=0.1345,team15_priority=priority,team15_sequence_id=teamrow.get('sequence_id',''),previous_rank=old.get('rank',''),previous_candidate_number=old.get('candidate_number',''),team15_original_description=teamrow.get('design_description',''))
 if priority:row['reason_picked']='Team priority #'+str(priority)+'; '+row['reason_picked']
 row['_sort']=((0 if priority and status!='rejected' else 1),int(priority or 999),group,0 if m else 1,-p,-safety if len(c)==10 else -float(tm.get(h) or 0),sid)
 ranked.append(row)
ranked.sort(key=lambda r:r['_sort']);counts=collections.Counter()
for i,r in enumerate(ranked,1):
 r['rank']=i;r['candidate_number']=f'AB{i:03d}'
 if r['screen_status']=='pass':counts[r['track']]+=1;r['track_rank']=counts[r['track']]
 r.pop('_sort')

def write(path,rows,fields=None):
 if not rows:return
 with (OUT/path).open('w') as f:
  w=csv.DictWriter(f,fieldnames=fields or list(rows[0]),extrasaction='ignore');w.writeheader();w.writerows(rows)
write('all_candidates_ranked.csv',ranked)
screenrows=[{k:d[k] for k in COLUMNS} for d in sorted(pool.values(),key=lambda d:d['sequence_id'])]
write('screen_input.csv',screenrows,COLUMNS)
write('controls_not_for_submission.csv',[{k:d[k] for k in COLUMNS} for d in controls.values()],COLUMNS)
if final:
 valid=[r for r in ranked if r['screen_status']=='pass']
 lookup={d['sequence_id']:d for d in pool.values()}
 submission=[{k:lookup[r['sequence_id']][k] for k in COLUMNS} for r in valid]
 assert len(submission)==len(accepted)
 assert {r['sequence_id'] for r in submission}==accepted
 write('submission_all_valid_ranked.csv',submission,COLUMNS)
 write('valid_candidates_ranked_with_notes.csv',valid)
 write('antibody_designs_ranked.csv',valid,['candidate_number','rank','track_rank','sequence_id','track','format','cdrh3','framework_edits','reason_picked','background_notes','risk_notes','epitope_evidence','nanomelt_tm_C','model_p_below_0_146_nM','model_pred_kd_nM','heavy_sequence','light_sequence','original_ids','sources','team15_priority','team15_sequence_id','previous_rank','previous_candidate_number','team15_original_description'])
 for track in ['optimisation','de_novo']:
  rr=[r for r in valid if r['track']==track]
  write(track+'_ranked_with_notes.csv',rr)
  write(track+'_submission_ranked.csv',[{k:lookup[r['sequence_id']][k] for k in COLUMNS} for r in rr],COLUMNS)
 report=dict(unique_candidates=len(pool),controls=len(controls),formats=dict(collections.Counter(d['kind'] for d in pool.values())),screen_pass=len(accepted) if final else None,screen_reject=len(pool)-len(accepted) if final else None,valid_tracks=dict(counts),source_counts=dict(collections.Counter(s for d in pool.values() for s in d['sources'])),enumeration_pool_included=False,ranking_is_heuristic=True)
else:report=dict(unique_candidates=len(pool),controls=len(controls),screen_status='pending',enumeration_pool_included=False)
(OUT/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2));print('TOP',[(r['rank'],r['sequence_id'],r['track'],r['screen_status']) for r in ranked[:10]])
