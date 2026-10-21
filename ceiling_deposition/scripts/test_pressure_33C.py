from pathlib import Path
import subprocess,shutil,re,json,time
import numpy as np
R=Path(__file__).resolve().parents[1];src=R/'flow_sensitivity/corner_both_closed_T33_Q30';out=R/'scenarios/ceiling33_Q30/linear_solver_tests';out.mkdir(exist_ok=True)
rows=[]
for name,smoother,direct in [('baseline','GaussSeidel',False),('dic','DICGaussSeidel',True),('symmetric','symGaussSeidel',True)]:
 c=R/'flow_sensitivity'/('linear_T33_'+name)
 if not c.exists():
  c.mkdir()
  for d in ['constant','system','40']:shutil.copytree(src/d,c/d)
  p=c/'system/controlDict';s=p.read_text();s=re.sub(r'endTime\s+[^;]+;','endTime 40.5;',s);p.write_text(s)
  p=c/'system/fvSolution';s=p.read_text().replace('smoother GaussSeidel;',f'smoother {smoother};')
  if direct:s=s.replace('cacheAgglomeration true;','cacheAgglomeration true; directSolveCoarsest yes;')
  p.write_text(s)
 log=R/'verification'/('linear_T33_'+name+'.log');t=time.monotonic()
 if not (c/'40.5/U').exists():
  with log.open('w') as f:subprocess.run(['docker','exec','kitchen-openfoam-2406','bash','-lc',f'{R}/bin/kitchenCeilingAirFoam -case {c}'],stdout=f,stderr=subprocess.STDOUT,check=True)
 cpu=float(re.findall(r'ExecutionTime = ([\d.eE+-]+)',log.read_text())[-1]);rows.append(dict(case=name,cpu_seconds=cpu,wall_seconds=time.monotonic()-t));print(rows[-1],flush=True)
(out/'timing.json').write_text(json.dumps(rows,indent=2))
