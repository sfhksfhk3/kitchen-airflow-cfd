"""Matched 80–100s low-flow intervention: mathematical actual Q24 vs Q30.
A zero-initial-aerosol 20s source pulse, NOT an entire cooking event comparison.
"""
from pathlib import Path
import re,shutil,json,subprocess
import pandas as pd
from update_33C_Q30 import ROOT,S,pc
from run_33C_Q30 import run_flow
from run_suite import run_one
src=ROOT/'flow_sensitivity/corner_both_closed_T33_Q30';c=ROOT/'flow_sensitivity/corner_closed_T33_Q24_window'
if not c.exists():
 c.mkdir()
 for d in ['constant','system','80']:shutil.copytree(src/d,c/d)
 for f in ['mesh_reference.npz','boundary_update.json']:shutil.copy2(src/f,c/f)
 p=c/'80/U';p.write_text(re.sub(r'volumetricFlowRate\s+constant\s+[^;]+;','volumetricFlowRate constant 0.4;',p.read_text()))
 p=c/'constant/aerosolProperties';p.write_text(re.sub(r'hoodTargetQ\s+[^;]+;','hoodTargetQ 0.4;',p.read_text()))
 p=c/'system/controlDict';t=p.read_text();t=re.sub(r'endTime\s+[^;]+;','endTime 100;',t);p.write_text(t)
 (c/'intervention.json').write_text(json.dumps(dict(source=str(src),initial_air_checkpoint_s=80,end_s=100,actual_target_m3_min=24,control_target_m3_min=30,initial_aerosol='zero in both replay pulses; earlier emitted particles excluded',note='Step down from common Q30 air state. Not a full-meal or distinct steady-Q24 initial state.'),indent=2))
log=ROOT/'verification/fan_window_33C_air.log'
if not (c/'100/U').exists():
 with log.open('w') as f:subprocess.run(['docker','exec','kitchen-openfoam-2406','bash','-lc',f'{ROOT}/bin/kitchenCeilingAirFoam -case {c}'],stdout=f,stderr=subprocess.STDOUT,check=True)
rows=[]
for source,q in [(c,24),(src,30)]:
 n=f'corner_closed_T33_Q{q}_pulse80_100';spec=dict(source=str(source),name=n,start=80,end=100,dt=.1,stride=2,sourceStopTime=100,impactionMode='parallel')
 run_one(spec)
 if not (pc.OUT/n/'summary.json').exists():pc.process(ROOT/'runs'/n)
 s=json.loads((pc.OUT/n/'summary.json').read_text());rows.append(dict(case=n,Q_late_m3_min=s['Q_late_mean_m3_s']*60,M_unit_mg_per_mg=s['M_total_unit_mg_per_mg'],M_reference_mg=s['M_total_reference_mg'],input_reference_mg=s['input_reference_mg'],mass_error=s['max_relative_mass_error'],hood_fraction=s['hood_fraction'],other_wall_fraction=s['other_wall_fraction'],airborne_fraction=s['airborne_fraction'],ceiling_fraction=s['ceiling_fraction']))
a,b=rows
value=100*(a['M_unit_mg_per_mg']-b['M_unit_mg_per_mg'])/a['M_unit_mg_per_mg'] if a['M_unit_mg_per_mg'] else None
for r in rows:r.update(model_deposition_reduction_percent_Q24_to_Q30=value,scope='80–100s source pulse; initial oil mass zero; common Q30 air initial state at80s, lower-Q branch steps down to24; 33C ceiling. Not full-event or grid-validated improvement.')
pd.DataFrame(rows).to_csv(S/'fan_window_comparison.csv',index=False);print(json.dumps(rows,ensure_ascii=False),flush=True)
