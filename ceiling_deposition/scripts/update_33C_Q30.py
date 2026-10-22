"""Idempotent full update. Never replace historical runs, reports or parameters."""
from pathlib import Path
import re,json,time,sys,subprocess,concurrent.futures,os
import pandas as pd
from run_33C_Q30 import ROOT,S,run_flow
from run_suite import run_one
import process_ceiling as pc
import clearance_integral as ci
PARAM=json.loads((S/'parameters.json').read_text())
REF=json.loads((ROOT/'config/emission_scenario.json').read_text());REF.update(label_zh='中式住家一次有機氣膠排放外推，非純油脂實測',rate_per_burner_mg_min=1.8,total_rate_mg_min=3.6,total_rate_ug_s=60,interpretation='1.8mg/min/head is an explicitly selected residential Chinese primary-organic-aerosol example, extrapolated to two simultaneous burners. Mass-median2.5um and assumed GSD1.5 retained; not pure grease or PM2.5 measurement.',rate_sensitivity_factors=[1.1/1.8,1,2.2/1.8,7.2/1.8])
(S/'emission_scenario.json').write_text(json.dumps(REF,ensure_ascii=False,indent=2))
pc.REF=REF;ci.REF=REF
pc.OUT=S/'results';pc.OUT.mkdir(exist_ok=True);ci.OUT=pc.OUT

def finished_air(c):
 for q in Path('/proc').glob('[0-9]*/cmdline'):
  try:a=q.read_bytes().split(b'\0')
  except OSError:continue
  if a and a[0]==str(ROOT/'bin/kitchenCeilingAirFoam').encode() and str(c).encode() in a:return False
 return True

def replay(c,end=240,alpha=1,dt=.1,stride=2,suffix=''):
 name=c.name+('_stick'+str(alpha) if alpha!=1 else '')+suffix
 spec=dict(source=str(c),name=name,end=end,dt=dt,stride=stride,ceilingSticking=alpha,impactionMode='parallel')
 run_one(spec)
 out=pc.OUT/name
 if not (out/'summary.json').exists():pc.process(ROOT/'runs'/name)
 if end==240:ci.merge(ci.complete(name))
 return name

def task(c):
 while not finished_air(c):time.sleep(5)
 a=run_flow(c)
 n=replay(c)
 if 'both_closed' in c.name:
  for alpha in [.1,.5]:replay(c,alpha=alpha)
 return a,n

def main():
 cases=[ROOT/'flow_sensitivity'/f'{l}_{d}_T33_Q30' for d in ['both_closed','both_open','left_open','right_open'] for l in ['corner','center']]
 with concurrent.futures.ThreadPoolExecutor(2) as pool:
  futures=[pool.submit(task,c) for c in cases];records=[]
  for f in concurrent.futures.as_completed(futures):
   a,n=f.result();records.append(a);pd.DataFrame(records).to_csv(S/'flow_validation.csv',index=False)
   (S/'PROGRESS.json').write_text(json.dumps(dict(status='RUNNING',completed=len(records),expected=8,cases=records),indent=2))
 # A full-history aerosol timestep and frame-spacing test on the same air field.
 c=ROOT/'flow_sensitivity/corner_both_closed_T33_Q30'
 replay(c,dt=.05,suffix='_dt005');replay(c,stride=4,suffix='_frames4')
 # Low-cost mixed mesh/startup diagnostic: fine retains its initial flux projection. Not an isolated grid convergence test.
 fine=ROOT/'flow_sensitivity/corner_both_closed_T33_Q30_fine'
 while not finished_air(fine):time.sleep(5)
 if (fine/'20/U').exists():
  replay(c,end=20,suffix='_startup20');replay(fine,end=20,suffix='_startup20')
 subprocess.run([sys.executable,str(ROOT/'scripts/report_33C_Q30.py')],check=True)
 (S/'PROGRESS.json').write_text(json.dumps(dict(status='COMPLETE_PHYSICAL_VALIDATION_LIMITED',completed=8),indent=2))
if __name__=='__main__':main()
