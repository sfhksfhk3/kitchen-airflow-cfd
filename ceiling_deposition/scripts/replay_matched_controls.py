import json,time
from pathlib import Path
from run_suite import ROOT,run_one
for layout in ['corner','center']:
 source=ROOT/'flow_sensitivity'/f'{layout}_both_closed_control'
 while not (source/'140/U').exists():time.sleep(5)
 starts=[float(p.name) for p in source.iterdir() if p.is_dir() and p.name.replace('.','',1).isdigit()];start=min(starts,key=lambda t:abs(t-100))
 print(json.dumps(run_one(dict(source=str(source),name=f'{layout}_fanwindow_control_parallel',start=start,end=140,dt=.1,useSavedInitialAerosols=True,impactionMode='parallel'))),flush=True)
