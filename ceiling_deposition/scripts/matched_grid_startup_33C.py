"""Same rest initialization and time-step controls on existing coarse/fine grids.
A 0–10s startup diagnostic, not full-event grid convergence.
"""
from pathlib import Path
import shutil,re,subprocess,json
from update_33C_Q30 import ROOT,S,pc
from run_suite import run_one
base=ROOT/'flow_sensitivity/corner_both_closed_T33_Q30_fine';c=ROOT/'flow_sensitivity/corner_T33_Q30_fine_matched10'
if not c.exists():
 c.mkdir()
 for d in ['constant','system','0']:shutil.copytree(base/d,c/d)
 for f in ['mesh_reference.npz','boundary_update.json']:shutil.copy2(base/f,c/f)
 p=c/'constant/aerosolProperties';p.write_text(re.sub(r'projectInitialMassFlux\s+[^;]+;','projectInitialMassFlux false;',p.read_text()))
 p=c/'system/controlDict';s=p.read_text();s=s.replace('startFrom latestTime','startFrom startTime');s=re.sub(r'endTime\s+[^;]+;','endTime 10;',s);s=re.sub(r'deltaT\s+[^;]+;','deltaT 0.002;',s);p.write_text(s)
log=ROOT/'verification/matched_grid_startup33_air.log'
if not (c/'10/U').exists():
 with log.open('w') as f:subprocess.run(['docker','exec','kitchen-openfoam-2406','bash','-lc',f'{ROOT}/bin/kitchenCeilingAirFoam -case {c}'],stdout=f,stderr=subprocess.STDOUT,check=True)
ss=[]
for source,tag in [(ROOT/'flow_sensitivity/corner_both_closed_T33_Q30','coarse'),(c,'fine')]:
 n=f'corner_T33_Q30_matched10_{tag}';spec=dict(source=str(source),name=n,start=0,end=10,dt=.1,stride=2,sourceStopTime=10,impactionMode='parallel')
 run_one(spec)
 if not (pc.OUT/n/'summary.json').exists():pc.process(ROOT/'runs'/n)
 ss.append(json.loads((pc.OUT/n/'summary.json').read_text()))
a,b=ss
(S/'matched_grid_startup10.json').write_text(json.dumps(dict(coarse=a,fine=b,ratio=b['M_total_unit_mg_per_mg']/a['M_total_unit_mg_per_mg'],scope='0–10s, both start from rest and zero oil, no initial flux projection, common maxCo12/maxDt0.04/dtInitial0.002; actual adaptive air time steps differ. Not a full-event mesh/GCI validation.'),indent=2))
