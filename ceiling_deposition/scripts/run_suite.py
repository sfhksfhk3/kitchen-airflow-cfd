from pathlib import Path
import json,subprocess,time,concurrent.futures,sys,os
import pandas as pd
from prepare_replay import prepare,ROOT

def run_one(spec):
 name=spec['name'];case=ROOT/'runs'/name;out=ROOT/'verification'/f'{name}.log'
 if (case/'RUN_COMPLETE.json').exists():
  saved=json.loads((case/'replay_manifest.json').read_text())
  for k,v in spec.items():
   if k in saved['parameters'] and saved['parameters'][k]!=v:raise ValueError(f'{name}: changed {k}; choose a new case name')
  return json.loads((case/'RUN_COMPLETE.json').read_text())
 if not case.exists():prepare(**spec)
 t=time.monotonic()
 # A stopped scheduler does not invalidate a completed solver. Adopt its output.
 def active():
  for q in Path('/proc').glob('[0-9]*/cmdline'):
   try:args=q.read_bytes().split(b'\0')
   except OSError:continue
   if args and args[0]==str(ROOT/'bin/ceilingReplayFoam').encode() and str(case).encode() in args:return True
  return False
 while active():time.sleep(5)
 already=(out.exists() and 'End replay' in out.read_text()[-1000:])
 if not already:
  with out.open('w') as f:
   r=subprocess.run(['docker','exec','kitchen-openfoam-2406','bash','-lc',f'{ROOT}/bin/ceilingReplayFoam -case {case}'],stdout=f,stderr=subprocess.STDOUT)
 if not already and r.returncode:raise RuntimeError(f'{name} failed; see {out}')
 d=pd.read_csv(case/'diagnostics/ceiling_budget.csv');assert len(d)>0
 total=d.initial_kg+d.released_kg
 err=float((d.mass_error_kg.abs()/total.clip(lower=1e-30)).max());assert err<1e-7,(name,err)
 assert d.min_Y.min()>-1e-18
 last=d.groupby('diameter_um').tail(1);assert abs(last.time_s.min()-spec.get('end',240))<1e-6
 result=dict(name=name,wall_seconds=time.monotonic()-t,max_relative_mass_error=err,end_s=float(last.time_s.min()),status='PASS')
 (case/'RUN_COMPLETE.json').write_text(json.dumps(result,indent=2));print(json.dumps(result),flush=True);return result

def suite():
 base=Path(os.environ.get('KITCHEN_FLOW_ROOT','/workspace/kitchen_openfoam'))/'cases';a=[]
 for layout in ['corner','center']:
  for doors in ['both_open','both_closed','left_open','right_open']:
   a.append(dict(source=str(base/f'{layout}_{doors}'),name=f'{layout}_{doors}'))
   a.append(dict(source=str(base/f'{layout}_{doors}'),name=f'{layout}_{doors}_parallel',impactionMode='parallel'))
 for layout in ['corner','center']:
  src=str(base/f'{layout}_both_closed')
  for alpha in [.1,.5]:a.append(dict(source=src,name=f'{layout}_closed_stick{alpha:g}',ceilingSticking=alpha,impactionMode='parallel'))
  a.append(dict(source=src,name=f'{layout}_closed_dt005',dt=.05,impactionMode='parallel'))
 c=str(base/'corner_both_closed')
 for name,opt in [('frames20',dict(stride=20)),('ScT05',dict(SchmidtT=.5)),('ceiling_cold5',dict(wallTemperatureOffsetK=-5)),('ceiling_warm5',dict(wallTemperatureOffsetK=5)),('no_thermo',dict(thermophoresisFactor=0)),('no_inertia',dict(inertiaFactor=0,impactionFactor=0)),('coeff05',dict(coefficientUpdateSeconds=.5)),('quad144',dict(wallQuadrature=144))]:
  a.append(dict(source=c,name='corner_window_'+name,start=110,end=140,impactionMode='parallel',**opt))
 for layout in ['corner','center']:
  a.append(dict(source=str(base/f'{layout}_both_closed'),name=f'{layout}_gridwindow_coarse',start=110,end=140,dt=.1,impactionMode='parallel'))
  a.append(dict(source=f'/workspace/kitchen_openfoam/sensitivity/{layout}_both_closed_grid',name=f'{layout}_gridwindow_fine',start=110,end=140,dt=.1,impactionMode='parallel'))
 return a
if __name__=='__main__':
 pilot=ROOT/'runs/pilot_final/diagnostics/ceiling_budget.csv';p=pd.read_csv(pilot);assert (p.mass_error_kg.abs()/p.released_kg).max()<1e-7
 (ROOT/'verification/pilot_validation.json').write_text(json.dumps({'passed':True,'duration_s':20,'dt_s':.2,'relative_mass_error':float((p.mass_error_kg.abs()/p.released_kg).max()),'min_Y':float(p.min_Y.min())},indent=2))
 specs=suite();(ROOT/'config/suite.json').write_text(json.dumps(specs,indent=2))
 with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
  futures=[pool.submit(run_one,s) for s in specs]
  results=[]
  for f in concurrent.futures.as_completed(futures):
   try:results.append(f.result())
   except Exception as e:print('FAIL',str(e),flush=True);raise
   (ROOT/'verification/suite_progress.json').write_text(json.dumps(results,indent=2))
