"""Independent checks of exported areas, units, mass budgets and VTK arrays.

These checks verify output integrity, not physical mesh convergence.
Run after report_33C_Q30.py; no CFD work is repeated.
"""
from pathlib import Path
import json,re,xml.etree.ElementTree as ET
import numpy as np
import pandas as pd
import vtk
from vtk.util.numpy_support import vtk_to_numpy

ROOT=Path(__file__).resolve().parents[1]
S=ROOT/'scenarios/ceiling33_Q30'

def main():
 cases=[f'{l}_{d}_T33_Q30' for d in ['both_open','both_closed','left_open','right_open'] for l in ['corner','center']]
 cases += [f'{l}_both_closed_T33_Q30_stick{a}' for l in ['corner','center'] for a in [.1,.5]]
 records=[]
 for c in cases:
  out=S/'results'/c
  f=pd.read_csv(out/'ceiling_complete_event_faces.csv')
  s=json.loads((out/'complete_event_summary.json').read_text())
  assert np.isfinite(f.select_dtypes('number').to_numpy()).all(),c
  assert (f.area_m2>0).all() and np.isclose(f.area_m2.sum(),4.245,atol=1e-12),c
  assert (f.M_unit_mg_m2_per_mg>=0).all(),c
  total=float(f.M_unit_mg_m2_per_mg@f.area_m2)
  assert np.isclose(total,s['M_total_unit_mg_per_mg'],rtol=1e-10,atol=1e-300),c
  assert np.allclose(f.M_reference_mg_m2,7.2*f.M_unit_mg_m2_per_mg,rtol=1e-10,atol=1e-300),c
  assert np.isclose(s['M_mean_unit_mg_m2_per_mg'],total/f.area_m2.sum(),rtol=1e-10,atol=1e-300),c
  fate=sum(s[k] for k in ['ceiling_fraction','other_wall_fraction','hood_fraction','left_fraction','right_fraction','airborne_fraction'])
  assert abs(fate-1)<1e-7,(c,fate)
  reader=vtk.vtkXMLPolyDataReader();reader.SetFileName(str(out/'ceiling_complete_event.vtp'));reader.Update();poly=reader.GetOutput()
  assert poly.GetNumberOfCells()==len(f),c
  a=poly.GetCellData().GetArray('M_unit_mg_m2_per_mg')
  assert a is not None,c
  assert np.allclose(vtk_to_numpy(a),f.M_unit_mg_m2_per_mg,rtol=1e-10,atol=1e-300),c
  pvd=next(out.glob('*.pvd'));times=[float(x.attrib['timestep']) for x in ET.parse(pvd).findall('.//DataSet')]
  assert all(t<=240.000001 for t in times) and all(a<b for a,b in zip(times,times[1:])),(c,times)
  source=Path(json.loads((ROOT/'runs'/c/'replay_manifest.json').read_text())['source'])
  t=(source/'240/T').read_text();block=re.search(r'\bceiling\s*\{([^{}]*)\}',t)[1]
  assert re.search(r'type\s+fixedValue\s*;',block) and re.search(r'value\s+uniform\s+306\.15\s*;',block),c
  records.append(dict(case=c,faces=len(f),area_m2=float(f.area_m2.sum()),total_unit_mg_per_mg=total,budget_fraction_sum=fate,last_physical_VTK_time_s=max(times),status='PASS'))
 long=pd.read_csv(S/'results/long_term_scenarios.csv')
 base=long[(long.material=='gypsum_matte')&(long.emitted_mg_per_event==1)&(long.clean_efficiency==0)]
 for row in base.itertuples():
  s=json.loads((S/'results'/row.case/'complete_event_summary.json').read_text())
  assert np.isclose(row.retained_total_mg,s['M_total_unit_mg_per_mg']*row.days*row.events_per_day,rtol=1e-10,atol=1e-300)
 # Shared effective attachment factors cannot produce a fabricated material ranking.
 keys=[c for c in long.columns if c not in ['material','retained_total_mg']]
 assert long.groupby(keys,dropna=False).retained_total_mg.nunique().max()==1
 result=dict(status='PASS',scope='Export integrity only; physical deposition is not validated by these checks.',cases=records,long_term_rows=len(long))
 (S/'EXPORT_VALIDATION.json').write_text(json.dumps(result,ensure_ascii=False,indent=2))
 print(json.dumps(result,ensure_ascii=False,indent=2))

if __name__=='__main__':main()
