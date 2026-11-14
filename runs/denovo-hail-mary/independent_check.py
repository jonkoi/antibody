"""Unguided co-fold of screen-passing de novo VHHs, with positive/failed-route controls."""
import csv,json,os,sys,urllib.request,io,zipfile,statistics
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];OUT=Path(__file__).resolve().parent
for l in (ROOT/'.env').read_text().splitlines():
 if '=' in l and not l.startswith('#'):
  k,v=l.split('=',1);os.environ[k]=v.strip().strip('\"\'')
sys.path.insert(0,str(ROOT/'pipeline/denovo'))
from jf import get,post,BASE,KEY
from epitope import read_atoms,chain_seq,contacts,TRAST_EPI,CORE4

def save(name,value):(OUT/name).write_text(json.dumps(value,indent=2)+'\n')

def target_sequences():
 a=read_atoms(str(ROOT/'runs/final-shortlist/boltz2_unguided/ctl_wt.pdb'))
 return chain_seq(a,'C'),chain_seq(a,'D')

act=sys.argv[1]
confirm=len(sys.argv)>2 and sys.argv[2]=='confirm'
RUN='confirmation' if confirm else 'fold'
PREFIX='confirmation_' if confirm else ''
if act=='submit':
 if (OUT/(RUN+'_job.json')).exists():raise SystemExit('Refusing duplicate submission')
 novelty={r['sequence_id']:r for r in json.loads((OUT/'novelty.json').read_text())}
 accepted=json.loads((OUT/'accepted.json').read_text())
 chosen=[r for r in accepted if novelty[r['sequence_id']]['novelty_class']=='de_novo']
 if confirm:chosen=[r for r in chosen if r['sequence_id']=='dn_vhh_her2_hail_mary_1']
 assert chosen,'No eligible de novo candidates'
 c,d=target_sequences();targets=[]
 for r in chosen:
  targets.append(dict(name=r['sequence_id'],content=f">A\n{r['heavy_sequence']}\n>C\n{c}\n>D\n{d}\n"))
 sys.path.insert(0,str(ROOT/'pipeline'));from build_submission_table import WT_HEAVY,WT_LIGHT
 targets.append(dict(name='ctl_wt',content=f">A\n{WT_LIGHT}\n>B\n{WT_HEAVY}\n>C\n{c}\n>D\n{d}\n"))
 old=read_atoms(str(ROOT/'runs/denovo-vhh/boltz2_unguided/cont_full_s1_r15.pdb'))
 if not confirm:targets.append(dict(name='ctl_old_offsite',content=f">A\n{chain_seq(old,'H')}\n>C\n{c}\n>D\n{d}\n"))
 payload=dict(model='boltz2',name='team15-denovo-confirm' if confirm else 'team15-one-shot-denovo-check',targets=targets,params=dict(use_msa_server=True,diffusion_samples=3,seed=1073 if confirm else 73,output_format='pdb'))
 save((RUN+'_request.json'),payload);res=post('/v1/predictions',payload,timeout=300);save((RUN+'_job.json'),res);print(json.dumps(res)[:700])
elif act=='poll':
 job=json.loads((OUT/(RUN+'_job.json')).read_text());r=get('/v1/jobs/'+job['id']+'?wait=45',timeout=65);save((RUN+'_status.json'),r);print(json.dumps({k:r.get(k) for k in ['id','status','done','total','error']},indent=2))
elif act=='collect':
 job=json.loads((OUT/(RUN+'_job.json')).read_text());jid=job['id']
 res=get('/v1/jobs/'+jid+'/results',timeout=180);save((RUN+'_results.json'),res)
 req=urllib.request.Request(BASE+'/v1/jobs/'+jid+'/archive',headers={'Authorization':'Bearer '+KEY})
 raw=urllib.request.urlopen(req,timeout=120).read();folder=OUT/('confirmation' if confirm else 'unguided');folder.mkdir(exist_ok=True)
 with zipfile.ZipFile(io.BytesIO(raw)) as z:
  for name in z.namelist():
   if name.endswith('.pdb'):(folder/Path(name).name).write_bytes(z.read(name))
 save((RUN+'_usage.json'),get('/v1/jobs/'+jid+'/usage'));print('Structures:',len(list(folder.glob('*.pdb'))))
elif act=='analyse':
 c,d=target_sequences();rows=[]
 for path in sorted((OUT/('confirmation' if confirm else 'unguided')).glob('*.pdb')):
  a=read_atoms(str(path));seqs={ch:chain_seq(a,ch) for ch in {v[0] for v in a}}
  name=path.stem.split('_model_')[0]
  request=json.loads((OUT/(RUN+'_request.json')).read_text())
  content=next(t['content'] for t in request['targets'] if t['name']==name)
  expected=[part.split('\n',1)[1].replace('\n','') for part in content.split('>')[1:]]
  assert sorted(seqs.values())==sorted(expected),(path,seqs)
  cch=next(ch for ch,s in seqs.items() if s==c);dch=next(ch for ch,s in seqs.items() if s==d)
  binders=[ch for ch in seqs if ch not in [cch,dch]];hit=set()
  for b in binders:
   for (ch,ri,rn),dist in contacts(a,[cch,dch],b).items():hit.add(ri+(488 if ch==cch else 590))
  epi=hit&set(TRAST_EPI);core=hit&set(CORE4);rec=len(epi)/19
  rows.append(dict(target=path.stem.split('_model_')[0],structure=path.name,recall=round(rec,4),core4=len(core),strict=int(rec>=0.85 and len(core)==4),vhh_overlap=int(rec>=0.5 and len(core)==4),hits=' '.join(map(str,sorted(epi)))))
 with (OUT/(PREFIX+'epitope_samples.csv')).open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
 summary=[]
 for target in sorted({r['target'] for r in rows}):
  rr=[r for r in rows if r['target']==target]
  summary.append(dict(target=target,n=len(rr),strict=sum(r['strict'] for r in rr),vhh_overlap=sum(r['vhh_overlap'] for r in rr),median_recall=statistics.median(r['recall'] for r in rr),min_recall=min(r['recall'] for r in rr),min_core4=min(r['core4'] for r in rr)))
 with (OUT/(PREFIX+'epitope_summary.csv')).open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(summary[0]));w.writeheader();w.writerows(summary)
 print(json.dumps(summary,indent=2))
