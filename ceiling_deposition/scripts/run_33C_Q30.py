"""Run the specified-actual-flow scenario, preserving previous experiments."""
from pathlib import Path
import json,subprocess,re,concurrent.futures,time
import pandas as pd
import numpy as np
from foam_mesh import mesh_geometry
ROOT=Path(__file__).resolve().parents[1]
S=ROOT/'scenarios/ceiling33_Q30';S.mkdir(exist_ok=True,parents=True)
def run_flow(c,end=240):
 p=c/'system/controlDict';t=p.read_text();t=t.replace('startFrom startTime','startFrom latestTime');t=re.sub(r'endTime\s+[^;]+;',f'endTime {end};',t);p.write_text(t)
 log=ROOT/'verification'/f'{c.name}_air.log'
 if not (c/f'{end:g}'/'U').exists():
  with log.open('w') as f:subprocess.run(['docker','exec','kitchen-openfoam-2406','bash','-lc',f'{ROOT}/bin/kitchenCeilingAirFoam -case {c}'],stdout=f,stderr=subprocess.STDOUT,check=True)
 ds=[pd.read_csv(p) for p in (c/'diagnostics').glob('aerosol_*.csv')];d=pd.concat(ds).drop_duplicates(['time_s','diameter_um']).sort_values('time_s');d=d[d.diameter_um==.1];late=d[d.time_s>max(10,end-20)]
 # Compare actual static pressure to ambient hydrostatic pressure at each
 # hood face height. 101325-p alone falsely includes ~19 Pa hydrostatic head.
 g=mesh_geometry(c);h=next(x for x in g['patches'] if x['name']=='hood');ids=np.arange(h['start'],h['start']+h['count']);area=g['area'][ids];ambient=101325*np.exp(-9.81*g['centres'][ids,2]/(287.05*298.15));pressure=[]
 for f in c.iterdir():
  if not f.is_dir() or not re.fullmatch(r'\d+(?:\.\d+)?',f.name) or not (f/'p').exists():continue
  if float(f.name)<end-20 or float(f.name)>end+1e-6:continue
  block=re.search(r'\bhood\s*\{([^{}]*)\}',(f/'p').read_text())[1]
  uniform=re.search(r'value\s+uniform\s+([^;]+);',block)
  if uniform:v=np.full(len(ids),float(uniform[1]))
  else:
   m=re.search(r'value\s+nonuniform\s+List<scalar>\s+(\d+)\s*\((.*?)\)\s*;',block,re.S);assert m,f;v=np.fromstring(m[2],sep=' ');assert len(v)==len(ids),(f,len(v),len(ids))
  pressure.append(dict(time_s=float(f.name),suction_vs_same_height_ambient_Pa=float(np.average(ambient-v,weights=area))))
 assert pressure,c
 (S/'pressure_history').mkdir(exist_ok=True)
 pd.DataFrame(pressure).sort_values('time_s').to_csv(S/'pressure_history'/(c.name+'.csv'),index=False)
 suction=float(np.mean([x['suction_vs_same_height_ambient_Pa'] for x in pressure]))
 result=dict(case=c.name,end_s=float(d.time_s.max()),Q_late_m3_min=float(late.Qhood_m3_s.mean()*60),hood_static_suction_late_Pa=suction,maximum_101325_minus_hood_pressure_Pa=float(d.hood_pressure_drop_Pa.max()),exceeds_fan_max_1300Pa=bool(suction>1300),air_mass_residual_abs_p95_kg_s=float(d.air_mass_residual_kg_s.abs().quantile(.95)),energy_residual_abs_p95_W=float(d.energy_balance_W.abs().quantile(.95)),ceiling_fixed_C=33)
 assert abs(result['Q_late_m3_min']-30)<.03,result
 assert result['end_s']>=end-1e-6
 (c/'FLOW_COMPLETE.json').write_text(json.dumps(result,indent=2));print(json.dumps(result),flush=True);return result
if __name__=='__main__':
 cases=[ROOT/'flow_sensitivity'/f'{l}_{d}_T33_Q30' for l in ['corner','center'] for d in ['both_closed','both_open','left_open','right_open']]
 with concurrent.futures.ThreadPoolExecutor(2) as pool:
  rs=list(pool.map(run_flow,cases))
 pd.DataFrame(rs).to_csv(S/'flow_validation.csv',index=False)
