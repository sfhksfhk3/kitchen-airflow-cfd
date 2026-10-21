from pathlib import Path
import json,concurrent.futures
from run_suite import ROOT,run_one
base=Path('/workspace/kitchen_openfoam/cases');specs=[]
for layout in ['corner','center']:
 source=base/f'{layout}_both_closed';ts=sorted(float(p.name) for p in source.iterdir() if p.is_dir() and p.name.replace('.','',1).isdigit());start=min(ts,key=lambda x:abs(x-100))
 specs.append(dict(source=str(source),name=f'{layout}_fanwindow_baseline',start=start,end=140,dt=.1,useSavedInitialAerosols=True))
 specs.append(dict(source=str(ROOT/'flow_sensitivity'/f'{layout}_both_closed_Qtarget120'),name=f'{layout}_fanwindow_boost',start=start,end=140,dt=.1,useSavedInitialAerosols=True))
specs.append(dict(source=str(ROOT/'flow_sensitivity/corner_both_closed_heatOff100'),name='corner_fanwindow_heatOff',start=100,end=140,dt=.1,useSavedInitialAerosols=True))
(ROOT/'config/interventions.json').write_text(json.dumps(specs,indent=2))
# Delay until the coupled-air intervention is actually complete.
import time
while not (ROOT/'flow_sensitivity/corner_both_closed_heatOff100/140/U').exists():time.sleep(5)
with concurrent.futures.ThreadPoolExecutor(max_workers=1) as p:
 for result in p.map(run_one,specs):print(json.dumps(result),flush=True)
