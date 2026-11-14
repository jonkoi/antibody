"""Compute presentation facts from raw local data; no new models or experiments."""
import csv,zipfile,io,json,statistics,math,collections
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2];OUT=Path(__file__).resolve().parent
wt='WGGDGFYAMD';counts=collections.defaultdict(collections.Counter)
with zipfile.ZipFile(ROOT/'data/her2-trastuzumab/tz_her2_facs_3classes_524346.csv.zip') as z:
 with z.open(z.namelist()[0]) as f:
  for r in csv.DictReader(io.TextIOWrapper(f)):
   counts[int(r['edit_distance'])][r['binding_class']]+=1
facs=[dict(substitutions=d,n=sum(c.values()),high=c['high'],high_rate=c['high']/sum(c.values())) for d,c in sorted(counts.items())]
with (OUT/'facs_full_table_by_distance.csv').open('w') as f:
 w=csv.DictWriter(f,fieldnames=list(facs[0]),lineterminator='\n');w.writeheader();w.writerows(facs)
spr=list(csv.DictReader((ROOT/'data/her2-trastuzumab/trastuzumab_her2_spr_all_768.csv').open()))
q=[r for r in spr if r.get('kd_nM') and len(r['cdrh3'])==10]
values=[float(r['kd_nM']) for r in q if r['cdrh3']==wt]
logsd=statistics.stdev(math.log10(x) for x in values)
priority=list(csv.DictReader((OUT/'team15_01.remote.csv').open()))
kept=[sum(r['heavy_sequence'][98+i]==wt[i] for i in [0,4,5,6]) for r in priority]
seqs=[r['heavy_sequence'][98:108] for r in priority[8:]]
pairmin=min(sum(x!=y for x,y in zip(a,b)) for i,a in enumerate(seqs) for b in seqs[i+1:])
facts=dict(facs_rows=sum(r['n'] for r in facs),facs_distance3=next(r for r in facs if r['substitutions']==3),spr_quantified_all_lengths=sum(bool(r.get('kd_nM')) for r in spr),spr_quantified_10mers=len(q),wt_kd_values_nM=values,wt_kd_median_nM=statistics.median(values),wt_kd_geometric_mean_nM=statistics.geometric_mean(values),wt_kd_range_nM=[min(values),max(values)],wt_replicate_log10_sd=logsd,rough_95pct_single_measurement_ratio=10**(1.96*math.sqrt(2)*logsd),team_priorities=len(priority),tier2_all4_contact_identities=sum(k==4 for k in kept[8:]),tier2_min_pairwise_hamming=pairmin)
(OUT/'verified_facts.json').write_text(json.dumps(facts,indent=2)+'\n');print(json.dumps(facts,indent=2))
