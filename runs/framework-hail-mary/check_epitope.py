#!/usr/bin/env python3
"""Small, unguided, matched Boltz-2 check; no affinity/Tm inference."""
import csv,json,os,sys,urllib.request
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
OUT=Path(__file__).resolve().parent
for line in (ROOT/'.env').read_text().splitlines():
    if '=' in line and not line.startswith('#'):
        k,v=line.split('=',1);os.environ[k]=v.strip().strip('\"\'')
sys.path.insert(0,str(ROOT/'pipeline/denovo'))
from jf import get,post
from epitope import read_atoms,chain_seq

def save(name,data):
    (OUT/name).write_text(json.dumps(data,indent=2)+'\n')

if sys.argv[1]=='submit':
    if (OUT/'job.json').exists(): raise SystemExit('Job already recorded; refusing duplicate submission')
    atoms=read_atoms(str(ROOT/'runs/final-shortlist/boltz2_unguided/ctl_wt.pdb'))
    c,d=chain_seq(atoms,'C'),chain_seq(atoms,'D')
    assert len(c)==92 and len(d)==17,(len(c),len(d))
    rows=list(csv.DictReader((OUT/'validation_panel.csv').open()))
    sys.path.insert(0,str(ROOT/'pipeline'))
    from build_submission_table import WT_HEAVY,WT_LIGHT
    rows.append(dict(sequence_id='ctl_wt',heavy_sequence=WT_HEAVY,light_sequence=WT_LIGHT))
    targets=[dict(name=r['sequence_id'],content=f">A\n{r['light_sequence']}\n>B\n{r['heavy_sequence']}\n>C\n{c}\n>D\n{d}\n") for r in rows]
    payload=dict(model='boltz2',name='team15-framework-hail-mary',params=dict(use_msa_server=True,diffusion_samples=3,output_format='pdb',seed=42),targets=targets)
    save('request.json',payload)
    response=post('/v1/predictions',payload,timeout=300)
    save('job.json',response);print(json.dumps(response)[:1500])
elif sys.argv[1]=='poll':
    job=json.loads((OUT/'job.json').read_text());jid=job['id']
    response=get('/v1/jobs/'+jid+'?wait=45',timeout=60)
    save('status.json',response);print(json.dumps(response)[:1800])
elif sys.argv[1]=='results':
    job=json.loads((OUT/'job.json').read_text())
    response=get('/v1/jobs/'+job['id']+'/results',timeout=180)
    save('results.json',response);print(json.dumps(response)[:2000])

elif sys.argv[1]=='download':
    import io,zipfile
    job=json.loads((OUT/'job.json').read_text())
    url=os.environ.get('JAPANFOLD_BASE_URL','https://api.japanfold.aiand.com').rstrip('/')+'/v1/jobs/'+job['id']+'/archive'
    req=urllib.request.Request(url,headers={'Authorization':'Bearer '+os.environ['JAPANFOLD_API_KEY']})
    raw=urllib.request.urlopen(req,timeout=120).read()
    folder=OUT/'structures';folder.mkdir(exist_ok=True)
    count=0
    with zipfile.ZipFile(io.BytesIO(raw)) as z:
        for name in z.namelist():
            if name.endswith('.pdb'):
                (folder/Path(name).name).write_bytes(z.read(name));count+=1
    save('usage.json',get('/v1/jobs/'+job['id']+'/usage'))
    print('Downloaded structures:',count)
elif sys.argv[1]=='analyse':
    import statistics
    from epitope import TRAST_EPI,CORE4,contacts
    rows=[]
    panel={r['sequence_id']:r for r in csv.DictReader((OUT/'validation_panel.csv').open())}
    for path in sorted((OUT/'structures').glob('*.pdb')):
        target=path.stem.split('_model_')[0]
        atoms=read_atoms(str(path))
        if target in panel:
            assert chain_seq(atoms,'B')==panel[target]['heavy_sequence'],path
            assert chain_seq(atoms,'A')==panel[target]['light_sequence'],path
        hit=set()
        for binder in ['A','B']:
            for (ch,ri,rn),dist in contacts(atoms,['C','D'],binder).items():
                hit.add(ri+(488 if ch=='C' else 590))
        recall=len(hit & set(TRAST_EPI))/len(TRAST_EPI)
        core=len(hit & set(CORE4))
        rows.append(dict(target=target,structure=path.name,epitope_recall=round(recall,4),core4=core,pass_epitope=int(recall>=0.85 and core==4),epitope_hits=' '.join(map(str,sorted(hit & set(TRAST_EPI))))))
    assert len(rows)==24,len(rows)
    with (OUT/'epitope_samples.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    summary=[]
    for target in sorted({r['target'] for r in rows}):
        rr=[r for r in rows if r['target']==target]
        summary.append(dict(target=target,n=len(rr),good=sum(r['pass_epitope'] for r in rr),min_recall=min(r['epitope_recall'] for r in rr),median_recall=statistics.median(r['epitope_recall'] for r in rr),min_core4=min(r['core4'] for r in rr)))
    with (OUT/'epitope_summary.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(summary[0]));w.writeheader();w.writerows(summary)
    print(json.dumps(summary,indent=2))
