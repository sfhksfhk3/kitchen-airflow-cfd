"""Run only the five new short coupled-flow cases, never the old main cases."""
from pathlib import Path
import concurrent.futures,subprocess
from prepare_flow_sensitivity import ROOT,prepare,prepare_heat_off
from run_matched_flow_controls import run as control

def run(case):
 case=Path(case)
 if (case/'140/U').exists():return str(case)
 with (ROOT/'verification'/f'{case.name}_coupled.log').open('w') as f:subprocess.run(['docker','exec','kitchen-openfoam-2406','bash','-lc',f'/workspace/kitchen_openfoam/bin/kitchenAerosolFoam -case {case}'],stdout=f,stderr=subprocess.STDOUT,check=True)
 return str(case)
if __name__=='__main__':
 cases=[prepare(l) for l in ['corner','center']]+[prepare_heat_off()]
 with concurrent.futures.ThreadPoolExecutor(max_workers=2) as p:
  for c in p.map(run,cases):print(c,flush=True)
 with concurrent.futures.ThreadPoolExecutor(max_workers=2) as p:
  for c in p.map(control,['corner','center']):print(c,flush=True)
