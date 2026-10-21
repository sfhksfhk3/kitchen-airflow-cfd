"""Generate quantitative validation and comparison tables from completed runs."""
from pathlib import Path
import json,numpy as np,pandas as pd,sys
from process_ceiling import ROOT,OUT

def get(name,full=False):
 p=OUT/name/('complete_event_summary.json' if full else 'summary.json')
 return json.loads(p.read_text()) if p.exists() else None

def main():
 rows=[];paired=[];fan=[]
 for layout in ['corner','center']:
  for doors in ['both_open','both_closed','left_open','right_open']:
   n=f'{layout}_{doors}_parallel';obs=get(n);full=get(n,True)
   if not obs:continue
   row={**obs,**(full or {}),'observed_to_240_mg':obs['M_total_reference_mg'],'observed_ceiling_percent':obs['ceiling_percent']};rows.append(row)
  ref=get(f'{layout}_both_closed_parallel');alt=get(f'{layout}_closed_dt005')
  if ref and alt:paired.append(dict(test=f'{layout}_time_step_0.1_to_0.05',scope='full 0–240 s aerosol replay',baseline=ref['M_total_reference_mg'],variant=alt['M_total_reference_mg'],relative_change_percent=100*(alt['M_total_reference_mg']/ref['M_total_reference_mg']-1)))
  coarse=get(f'{layout}_gridwindow_coarse');fine=get(f'{layout}_gridwindow_fine')
  if coarse and fine:paired.append(dict(test=f'{layout}_late_mesh_refinement',scope='110–140 s short source pulse, existing coarse/fine coupled flows; not full-history grid independence',baseline=coarse['M_total_unit_mg_per_mg'],variant=fine['M_total_unit_mg_per_mg'],relative_change_percent=100*(fine['M_total_unit_mg_per_mg']/coarse['M_total_unit_mg_per_mg']-1)))
  control=get(f'{layout}_fanwindow_control_parallel');boost=get(f'{layout}_fanwindow_boost_parallel')
  if control and boost:fan.append(dict(layout=layout,scope=f"{control['start_s']}–140 s matched air restart and old aerosol initial state; 120 s source off; not a whole-meal flow sweep",Q_baseline_m3_s=control['Q_late_mean_m3_s'],Q_boost_m3_s=boost['Q_late_mean_m3_s'],Q_increase_percent=100*(boost['Q_late_mean_m3_s']/control['Q_late_mean_m3_s']-1),deposition_baseline_reference_mg=control['M_total_reference_mg'],deposition_boost_reference_mg=boost['M_total_reference_mg'],deposition_reduction_percent=100*(1-boost['M_total_reference_mg']/control['M_total_reference_mg']),difference_mg=boost['M_total_reference_mg']-control['M_total_reference_mg']))
 ref=get('corner_gridwindow_coarse')
 if ref:
  for p in OUT.glob('corner_window_*/summary.json'):
   alt=json.loads(p.read_text());paired.append(dict(test=alt['case'],scope='110–140 s matched-window aerosol/closure sensitivity',baseline=ref['M_total_unit_mg_per_mg'],variant=alt['M_total_unit_mg_per_mg'],relative_change_percent=100*(alt['M_total_unit_mg_per_mg']/ref['M_total_unit_mg_per_mg']-1)))
 h0=get('corner_fanwindow_control_parallel');h1=get('corner_fanwindow_heatOff_parallel')
 if h0 and h1:paired.append(dict(test='heat_to_air_off_from100',scope='100–140 s matched coupled-air intervention; warm initial state and emissions until120 retained',baseline=h0['M_total_reference_mg'],variant=h1['M_total_reference_mg'],relative_change_percent=100*(h1['M_total_reference_mg']/h0['M_total_reference_mg']-1)))
 pd.DataFrame(rows).to_csv(OUT/'primary_complete_event_summary.csv',index=False)
 pd.DataFrame(paired).to_csv(OUT/'numerical_and_physical_sensitivity.csv',index=False)
 pd.DataFrame(fan).to_csv(OUT/'fan_intervention_comparison.csv',index=False)
 runs=[]
 for p in (ROOT/'runs').glob('*/RUN_COMPLETE.json'):runs.append(json.loads(p.read_text()))
 tails=[]
 for p in (ROOT/'verification').glob('*_clearance_*.json'):tails.append(json.loads(p.read_text()))
 # Read VTP in the same VTK library ParaView uses, not merely as XML.
 import vtk
 count=0
 for p in OUT.glob('*/ceiling_complete_event.vtp'):
  r=vtk.vtkXMLPolyDataReader();r.SetFileName(str(p));r.Update();v=r.GetOutput();assert v.GetNumberOfCells()>0
  for field in ['M_unit_mg_m2_per_mg','M_reference_mg_m2']:
   ar=v.GetCellData().GetArray(field);assert ar and ar.GetNumberOfTuples()==v.GetNumberOfCells()
   assert ar.GetRange()[0]>=0
  count+=1
 validation={'completed_transient_runs':len(runs),'completed_clearance_integrals':len(tails),'primary_cases_with_full_event':sum(get(f'{l}_{d}_parallel',True) is not None for l in ['corner','center'] for d in ['both_open','both_closed','left_open','right_open']),'worst_transient_relative_mass_error':max([x['max_relative_mass_error'] for x in runs],default=0),'worst_clearance_relative_mass_error':max([x['max_relative_mass_error'] for x in tails],default=0),'VTK_complete_event_files_read':count,'ceiling_grid_sensitivity':'FAILED_ROBUSTNESS_COARSE_FINE_LATE_BRANCHES_DIFFER_ABOUT_27X','full_history_grid_independence':'NOT_DEMONSTRATED','strict_mesh_checks':'FAILED_2_CHECKS_IN_REUSED_COARSE_MESH','experimental_ceiling_deposition_calibration':'NOT_AVAILABLE','condensable_vapour_model':'NOT_INCLUDED','native_replay_and_closed_forms':'PASS','physical_absolute_prediction_confidence':'LOW; conditional model response, not a measured grease exposure or paint-aging prediction'}
 markers=[]
 for H in ['1e+09','1e+12']:
  fp=ROOT/'verification'/f'corner_both_closed_parallel_clearance_{H}.json'
  if fp.exists():
   mm=json.loads(fp.read_text());bb=pd.read_csv(Path(mm['tail_case'])/'diagnostics/ceiling_budget.csv').sort_values('diameter_um');markers.append(float(np.array(json.loads((ROOT/'config/emission_scenario.json').read_text())['mass_weights'])@bb.ceiling_kg.to_numpy()))
 if len(markers)==2:validation['clearance_regularization_relative_change']=markers[0]/markers[1]-1
 (ROOT/'verification/final_validation.json').write_text(json.dumps(validation,indent=2));print(json.dumps(validation,indent=2))
 if len(rows)==8 and count>=12:
  assert all(abs(sum(r[k] for k in ['ceiling_fraction','other_wall_fraction','hood_fraction','left_fraction','right_fraction','airborne_fraction'])-1)<1e-7 for r in rows)
 if len(paired):print(pd.DataFrame(paired)[['test','relative_change_percent']].to_string(index=False))
 if fan:print(pd.DataFrame(fan).to_string(index=False))
if __name__=='__main__':main()
