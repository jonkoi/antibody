"""Read results of previously submitted HER2 design jobs; starts no jobs."""
from pathlib import Path
import os,sys,json,concurrent.futures
ROOT=Path(__file__).resolve().parents[2];OUT=Path(__file__).resolve().parent/'source'
for l in (ROOT/'.env').read_text().splitlines():
 if '=' in l and not l.startswith('#'):
  k,v=l.split('=',1);os.environ[k]=v.strip().strip('\"\'')
sys.path.insert(0,str(ROOT/'pipeline/denovo'));from jf import get
jobs={}
for p in (ROOT/'runs/denovo-vhh').glob('bg_jobs*.json'):
 for j in json.loads(p.read_text()):
  if not j['name'].startswith('ctl-lyz'):jobs[j['id']]=j

def fetch(j):
 path=OUT/('previous_job_'+j['id']+'.json')
 if path.exists():r=json.loads(path.read_text())
 else:
  r=get('/v1/jobs/'+j['id']+'/results',timeout=120);path.write_text(json.dumps(r,indent=2)+'\n')
 designs=r.get('designs') or []
 return dict(name=j['name'],id=j['id'],designs=len(designs),ready=r.get('ready'),path=str(path.relative_to(ROOT)))
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:summary=list(pool.map(fetch,jobs.values()))
(OUT/'previous_jobs_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
