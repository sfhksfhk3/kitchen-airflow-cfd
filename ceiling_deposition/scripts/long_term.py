"""Configurable repeated-event model. No yellowing, thickness or aging claims."""
from pathlib import Path
import json,numpy as np,pandas as pd
from process_ceiling import ROOT,OUT,size_weights,DS,REF

def cycle_units(days,events,interval,eff):
 stock=0.;peak=0.;before=0.;removed=0.
 for day in range(1,days+1):
  stock+=events;before=stock;peak=max(peak,before)
  if interval and day%interval==0:
   clean=eff*stock;removed+=clean;stock-=clean
 return stock,peak,before,removed

def main():
 cfg=json.loads((ROOT/'config/long_term.json').read_text());rows=[];sizes=[]
 for p in sorted(OUT.glob('*/summary.json')):
  s=json.loads(p.read_text());case=s['case']
  complete=p.parent/'complete_event_summary.json'
  if complete.exists():s.update(json.loads(complete.read_text()))
  ismain=case.endswith('_parallel') and 'fanwindow' not in case
  if ismain:alpha=1.
  elif '_closed_stick' in case:alpha=float(case.split('_closed_stick')[1])
  else:continue
  f=pd.read_csv(p.parent/('ceiling_complete_event_faces.csv' if complete.exists() else 'ceiling_faces.csv'));area=f.area_m2.to_numpy();m=f.M_unit_mg_m2_per_mg.to_numpy()
  for gsd in [1.2,1.5,2.0]:
   msize=f[[f'M_{d:g}um_unit_mg_m2_per_mg' for d in DS]].to_numpy()@size_weights(gsd)
   sizes.append(dict(case=case,gsd=gsd,ceiling_mg_per_emitted_mg=float(msize@area),mean_mg_m2_per_emitted_mg=float(msize@area/area.sum()),max_mg_m2_per_emitted_mg=float(msize.max())))
  for material in cfg['materials']:
   for events in cfg['daily_events']:
    for days in cfg['report_days']:
     for interval,eff in [(0,0)]+[(cfg['cleaning_interval_days'],e) for e in cfg['cleaning_efficiency_sensitivity']]:
      retained,peak,before,removed=cycle_units(days,events,interval,eff)
      for dose in cfg['emitted_mass_per_event_mg']:
       rows.append(dict(case=case,material=material['id'],material_zh=material['label_zh'],sticking_probability=alpha,events_per_day=events,days=days,emitted_mg_per_event=dose,clean_every_days=interval,clean_efficiency=eff,retained_total_mg=s['M_total_unit_mg_per_mg']*dose*retained,retained_mean_mg_m2=s['M_mean_unit_mg_m2_per_mg']*dose*retained,retained_p95_mg_m2=s['M_p95_unit_mg_m2_per_mg']*dose*retained,retained_max_mg_m2=s['M_max_unit_mg_m2_per_mg']*dose*retained,peak_before_clean_mean_mg_m2=s['M_mean_unit_mg_m2_per_mg']*dose*peak,removed_total_mg=s['M_total_unit_mg_per_mg']*dose*removed,note='same alpha has identical CFD deposition across materials; cleaning efficiency assumed, not material measured'))
 pd.DataFrame(rows).to_csv(OUT/'long_term_scenarios.csv',index=False)
 pd.DataFrame(sizes).to_csv(OUT/'size_distribution_sensitivity.csv',index=False)
 # Expose a machine-readable unit-response table before absolute mass scenarios.
 data=pd.DataFrame(rows)
 if len(data):data[data.emitted_mg_per_event==1].to_csv(OUT/'long_term_per_1mg_emitted.csv',index=False)
 print('long-term scenarios',len(rows),'size distribution cases',len(sizes))
if __name__=='__main__':main()
