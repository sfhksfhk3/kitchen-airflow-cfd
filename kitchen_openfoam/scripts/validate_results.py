"""Quantitative conservation, endpoint, export and late-stage sensitivity checks."""
from pathlib import Path
import json,numpy as np,csv
from foam_io import history,times
from emission_scenario import weights,DIAMETERS
ROOT=Path(__file__).resolve().parents[1]

def avg(t,v,a,b):
    inside=(t>a)&(t<b);tt=np.r_[a,t[inside],b];vv=np.r_[np.interp(a,t,v),v[inside],np.interp(b,t,v)]
    return float(np.trapezoid(vv,tt)/(b-a))

def stats(p,a=110,b=120):
    rows=history(p);rs=[]
    for d in DIAMETERS:
        r=[x for x in rows if abs(x['diameter_um']-d)<1e-8];t=np.array([x['time_s'] for x in r]);v=lambda name:np.array([x[name] for x in r])
        rs.append({'diameter_um':float(d),'stock_at_120_over_R_s':float(np.interp(120,t,v('stock_kg'))/1e-9),'walkway_mean_c_over_R_s_m3':avg(t,v('walkway_c_over_R_s_m3'),a,b),'max_particle_mass_error_relative':float(np.max(abs(v('mass_error_kg'))/np.maximum(v('released_kg'),1e-25))),'min_Y':float(min(v('min_Y'))),'Q_mean_m3_s':avg(t,v('Qhood_m3_s'),a,b),'Umax_mean_m_s':avg(t,v('Umax_m_s'),a,b),'Tmax_mean_C':avg(t,v('Tmax_K'),a,b)-273.15,'energy_rms_window_W':float(np.sqrt(avg(t,v('energy_balance_W')**2,a,b)))})
    return rs

def main():
    checks=[]
    for p in sorted((ROOT/'cases').iterdir()):
        if not (p/'diagnostics').exists():continue
        rows=history(p)
        if not rows:continue
        ts=times(p);s=stats(p);error=max(x['max_particle_mass_error_relative'] for x in s)
        base=[x for x in rows if x['diameter_um']==.1];t=np.array([x['time_s'] for x in base]);energy=np.array([x['energy_balance_W'] for x in base]);active=(t>10)&(t<120)
        out={'case':p.name,'end_time_history_s':max(t),'end_time_saved_s':max(ts),'particle_mass_error_max_relative':error,'minimum_Y':min(x['min_Y'] for x in s),'air_mass_balance_max_kg_s_after10':max(abs(x['air_mass_residual_kg_s']) for x in base if x['time_s']>10),'energy_rms_after10_cooking_W':float(np.sqrt(avg(t,energy**2,10,120))),'energy_residual_integral_fraction_input':float(np.trapezoid(energy,t)/(4000*120)),'source_release_end_kg_per_bin':base[-1]['released_kg'],'passed_particle_mass_balance':error<1e-6,'complete_history':max(t)>=239.99,'complete_saved_fields':max(ts)>=239.99}
        export=ROOT/'deliverables/data'/p.name/'field_export_validation.json'
        if export.exists():out['max_relative_export_mass_error']=max(abs(v) for x in json.loads(export.read_text()) for v in x['relative_field_mass_vs_history'])
        checks.append(out)
    sensitivities=[]
    for p in sorted((ROOT/'sensitivity').glob('*')):
        meta=p/'sensitivity_manifest.json'
        if not meta.exists() or not (p/'diagnostics').exists() or max(times(p))<139.99:continue
        m=json.loads(meta.read_text());src=Path(m['source_case']);orig=stats(src);new=stats(p)
        for a,b in zip(orig,new):
            sensitivities.append({'test':p.name,'kind':m['kind'],'diameter_um':a['diameter_um'],'parent_stock_R_120_s':a['stock_at_120_over_R_s'],'test_stock_R_120_s':b['stock_at_120_over_R_s'],'stock_relative_change':b['stock_at_120_over_R_s']/a['stock_at_120_over_R_s']-1,'parent_walkway_mean_R':a['walkway_mean_c_over_R_s_m3'],'test_walkway_mean_R':b['walkway_mean_c_over_R_s_m3'],'walkway_relative_change':b['walkway_mean_c_over_R_s_m3']/max(a['walkway_mean_c_over_R_s_m3'],1e-25)-1,'parent_energy_RMS_110_120_W':a['energy_rms_window_W'],'test_energy_RMS_110_120_W':b['energy_rms_window_W'],'test_max_particle_mass_error':b['max_particle_mass_error_relative'],'limitation':m['transfer'].get('limitation','matched-state late-stage timestep sensitivity')})
    out={'cases':checks,'sensitivity':sensitivities,'interpretation':'Conservation and numerical sensitivity, not experimental validation; large local changes must limit placement recommendations.'}
    (ROOT/'verification/final_validation.json').write_text(json.dumps(out,indent=2,default=lambda o:o.item() if hasattr(o,'item') else str(o)))
    if sensitivities:
        p=ROOT/'deliverables/numerical_sensitivity.csv'
        with p.open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(sensitivities[0]));w.writeheader();w.writerows(sensitivities)
    print(json.dumps({'completed':sum(x['complete_saved_fields'] for x in checks),'conservation_passed':sum(x['passed_particle_mass_balance'] for x in checks),'sensitivity_rows':len(sensitivities)},indent=2))

if __name__=='__main__':main()
