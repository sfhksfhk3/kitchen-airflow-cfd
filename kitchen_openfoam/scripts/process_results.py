"""Export independent size responses, mass scenario, raw cell fields and checks.

The exporter never launches or modifies a CFD solve.
"""
from pathlib import Path
import argparse,json,csv,numpy as np,meshio
from foam_io import history,fields,load_mesh,times
from emission_scenario import config,weights,DIAMETERS
from validate_results import avg
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'deliverables';OUT.mkdir(exist_ok=True)
LAYOUTS=['corner','center'];DOORS=['both_open','both_closed','left_open','right_open']

def csv_write(path,rows):
    if not rows:return
    with Path(path).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)

def stable_crossing(t,y,frac=.05):
    start=np.interp(120,t,y);mask=t>=120;t=t[mask]-120;y=y[mask]
    good=np.maximum.accumulate(y[::-1])[::-1]<=frac*start
    k=np.flatnonzero(good)
    return float(t[k[0]]) if len(k) else None

def process_case(p,export_frames=True):
    rows=history(p);all_t=sorted(set(r['time_s'] for r in rows));names=list(rows[0])
    a=np.array([[r[k] for k in names] for r in rows]);index={k:i for i,k in enumerate(names)}
    out=OUT/'data'/p.name;out.mkdir(parents=True,exist_ok=True)
    csv_write(out/'history_all_sizes.csv',rows)
    bins=[];arrays=[]
    for d in DIAMETERS:
        b=a[np.isclose(a[:,index['diameter_um']],d)];arrays.append(b)
        t=b[:,index['time_s']];source=b[:,index['released_kg']];j=np.argmin(abs(t-120));end=b[-1]
        last=(t>=90)&(t<=120)
        trapezoid=np.trapezoid
        bins.append({'case':p.name,'layout':p.name.split('_')[0],'doors':p.name.split('_',1)[1],'diameter_um':float(d),'end_time_s':float(t[-1]),'stock_120_over_R_s':float(np.interp(120,t,b[:,index['stock_kg']])/1e-9),'walkway_90_120_c_over_R_s_m3':avg(t,b[:,index['walkway_c_over_R_s_m3']],90,120),'walkway_120_c_over_R_s_m3':float(np.interp(120,t,b[:,index['walkway_c_over_R_s_m3']])),'hood_fraction_at_120':float(np.interp(120,t,b[:,index['hood_out_kg']])/(120e-9)),'hood_fraction_at_end':float(end[index['hood_out_kg']]/source[-1]),'door_escape_fraction_end':float((end[index['left_out_kg']]+end[index['right_out_kg']])/source[-1]),'deposited_fraction_end':float(end[index['deposited_kg']]/source[-1]),'remaining_fraction_end':float(end[index['stock_kg']]/source[-1]),'retention_integral_observed_s':float(trapezoid(np.r_[0,b[:,index['stock_kg']]],np.r_[0,t])/source[-1]),'stock_t95_after_stop_s':stable_crossing(t,b[:,index['stock_kg']]),'walkway_t95_after_stop_s':stable_crossing(t,b[:,index['walkway_c_over_R_s_m3']]),'max_relative_particle_mass_error':float(np.max(abs(b[:,index['mass_error_kg']])/np.maximum(source,1e-20))),'minimum_mass_fraction':float(b[:,index['min_Y']].min()),'Q_90_120_mean_m3_s':avg(t,b[:,index['Qhood_m3_s']],90,120),'Tair_max_90_120_mean_C':avg(t,b[:,index['Tmax_K']],90,120)-273.15})
    # Timestamps are identical across independent bin solves at each CFD step.
    arr=np.array(arrays);t=arr[0,:,index['time_s']];w=np.array(config()['mass_weights']);rate=config()['total_rate_ug_s'];source=arr[0,:,index['released_kg']]
    def mixcol(name):return np.einsum('b,bt->t',w,arr[:,:,index[name]])
    mc={name:mixcol(name) for name in ['stock_kg','hood_out_kg','left_out_kg','right_out_kg','deposited_kg','walkway_c_over_R_s_m3']}
    mix=[]
    for j,tt in enumerate(t):
        mix.append({'time_s':tt,'airborne_oil_mg':mc['stock_kg'][j]*rate*1e6,'released_oil_mg':source[j]*rate*1e6,'hood_out_oil_mg':mc['hood_out_kg'][j]*rate*1e6,'doors_out_oil_mg':(mc['left_out_kg'][j]+mc['right_out_kg'][j])*rate*1e6,'deposited_oil_mg':mc['deposited_kg'][j]*rate*1e6,'walkway_oil_ug_m3':mc['walkway_c_over_R_s_m3'][j]*rate,'Qhood_m3_s':arr[0,j,index['Qhood_m3_s']],'Tair_max_C':arr[0,j,index['Tmax_K']]-273.15,'Umax_m_s':arr[0,j,index['Umax_m_s']],'energy_balance_W':arr[0,j,index['energy_balance_W']],'hood_pressure_drop_Pa':arr[0,j,index['hood_pressure_drop_Pa']]})
    csv_write(out/'mixture_history_reference_scenario.csv',mix)
    stock=mixcol('stock_kg');walk=mixcol('walkway_c_over_R_s_m3');k=np.argmin(abs(t-120));last=(t>=90)&(t<=120)
    summary={'case':p.name,'layout':p.name.split('_')[0],'doors':p.name.split('_',1)[1],'end_time_s':float(t[-1]),'source_total_mg_min':config()['total_rate_mg_min'],'airborne_oil_at_120_mg':float(np.interp(120,t,stock)*rate*1e6),'walkway_90_120_oil_ug_m3':avg(t,walk,90,120)*rate,'walkway_120_oil_ug_m3':float(np.interp(120,t,walk)*rate),'hood_fraction_at_end':float(mixcol('hood_out_kg')[-1]/source[-1]),'door_escape_fraction_end':float((mixcol('left_out_kg')[-1]+mixcol('right_out_kg')[-1])/source[-1]),'deposited_fraction_end':float(mixcol('deposited_kg')[-1]/source[-1]),'remaining_fraction_end':float(stock[-1]/source[-1]),'retention_integral_observed_s':float(np.trapezoid(np.r_[0,stock],np.r_[0,t])/source[-1]),'airborne_t95_after_stop_s':stable_crossing(t,stock),'walkway_t95_after_stop_s':stable_crossing(t,walk),'Q_90_120_mean_m3_s':avg(t,arr[0,:,index['Qhood_m3_s']],90,120),'max_particle_mass_relative_error':max(b['max_relative_particle_mass_error'] for b in bins),'energy_balance_rms_cooking_after10_W':float(np.sqrt(np.mean(arr[0,(t>10)&(t<=120),index['energy_balance_W']]**2)))}
    (out/'summary.json').write_text(json.dumps(summary,indent=2))
    raw_checks=[]
    if export_frames:
        m=load_mesh(p);n=len(m['cells']);native_times=times(p)
        target_times=[0,10,30,60,100,120,130,150,180,240]
        chosen=[];requested=[]
        for target in target_times:
            pool=[v for v in native_times if 90<=v<=120+1e-7] if target==120 else native_times
            if not pool or target>native_times[-1]+5:continue
            tt=min(pool,key=lambda v:abs(v-target))
            if tt not in chosen:chosen.append(tt);requested.append(target)
        cooking_time=min([v for v in native_times if 90<=v<=120+1e-7],key=lambda v:abs(v-120))
        cache={};V=m['volumes']
        for tt in chosen:
            f=fields(p,tt,m);conc=f['C_over_R'];mc=conc@w*rate
            key=str(tt);cache[key]=np.c_[f['T']-273.15,f['speed'],f['p_rgh']-101325,conc,f['U']]
            if tt==cooking_time or abs(tt-240)<1e-7:
                cdata={'velocity_m_s':f['U'],'temperature_C':f['T']-273.15,'pressure_reduced_relative_Pa':f['p_rgh']-101325,'eddy_viscosity_m2_s':f['nut'],'oil_scenario_ug_m3':mc,**{f'C_over_R_{d:g}um_s_m3':conc[:,j] for j,d in enumerate(DIAMETERS)}}
                blocks=[];data={k:[] for k in cdata}
                for nv,typ in [(8,'hexahedron'),(6,'wedge')]:
                    ids=m['nverts']==nv
                    if not ids.any():continue
                    blocks.append((typ,m['cells'][ids,:nv]))
                    for k,val in cdata.items():data[k].append(val[ids])
                label=f'{tt:.8f}'.rstrip('0').rstrip('.').replace('.','p')
                meshio.write(out/f'fields_{label}s.vtu',meshio.Mesh(m['points'],blocks,cell_data=data))
                full=np.c_[m['centres'],V,f['U'],f['T']-273.15,f['p_rgh']-101325,f['nut'],conc,mc]
                np.savetxt(out/f'cells_{label}s.csv',full,delimiter=',',header='x_m,y_m,z_m,volume_m3,Ux_m_s,Uy_m_s,Uz_m_s,Tair_C,p_rgh_relative_Pa,nut_m2_s,C_R_0p1um_s_m3,C_R_0p3um_s_m3,C_R_1um_s_m3,C_R_3um_s_m3,C_R_10um_s_m3,oil_scenario_ug_m3',comments='')
                row_index=np.argmin(abs(t-tt));inferred=np.sum(conc*V[:,None],axis=0)*1e-9
                raw=arr[:,row_index,index['stock_kg']]
                raw_checks.append({'time_s':tt,'relative_field_mass_vs_history':((inferred-raw)/np.maximum(raw,1e-25)).tolist()})
        np.savez_compressed(out/'visualization_frames.npz',centres=m['centres'],volumes=V,times_s=chosen,target_times_s=requested,values=np.array([cache[str(tt)] for tt in chosen]))
        (out/'visualization_frames_metadata.json').write_text(json.dumps({'arrays':{'centres':'[cell,xyz] metres','volumes':'[cell] m3','times_s':'actual saved physical times','target_times_s':'requested viewer ticks; not actual save times','values':'[time,cell,column]'},'value_columns':['Tair_C','speed_m_s','p_rgh_relative_Pa','C_R_0p1um_s_m3','C_R_0p3um_s_m3','C_R_1um_s_m3','C_R_3um_s_m3','C_R_10um_s_m3','Ux_m_s','Uy_m_s','Uz_m_s'],'reference_emission_kg_s_per_independent_size':1e-9},indent=2))
        (out/'field_export_validation.json').write_text(json.dumps(raw_checks,indent=2))
    sensitivity=[]
    for gsd in config()['width_sensitivity_gsd']:
        ww=weights(gsd=gsd)
        sw=ww@arr[:,:,index['walkway_c_over_R_s_m3']]
        ss=ww@arr[:,:,index['stock_kg']]
        for factor in config()['rate_sensitivity_factors']:
            sensitivity.append({'case':p.name,'mass_median_um':2.5,'assumed_mass_GSD':gsd,'emission_rate_factor':factor,'total_emission_mg_min':config()['total_rate_mg_min']*factor,'walkway_90_120_oil_ug_m3':avg(t,sw,90,120)*rate*factor,'airborne_oil_at_120_mg':float(np.interp(120,t,ss)*rate*factor*1e6),'retention_integral_observed_s':float(np.trapezoid(np.r_[0,ss],np.r_[0,t])/source[-1]),'deposited_fraction_end':float(ww@arr[:,-1,index['deposited_kg']]/source[-1]),'method':'linear recombination of independently solved monodisperse responses; same thermal flow'})
    csv_write(out/'emission_distribution_sensitivity.csv',sensitivity)
    return bins,summary

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--partial',action='store_true');ap.add_argument('--skip-fields',action='store_true');a=ap.parse_args();allbins=[];summaries=[]
    for layout in LAYOUTS:
        for doors in DOORS:
            p=ROOT/'cases'/f'{layout}_{doors}'
            if not (p/'diagnostics').exists():
                if a.partial:continue
                raise RuntimeError(f'No computation for {p.name}')
            if max(times(p))<240-1e-7 and not a.partial:raise RuntimeError(f'Case incomplete: {p.name}')
            if a.partial and max(times(p))<120-1e-7:continue
            bins,s=process_case(p,not a.skip_fields);allbins.extend(bins);summaries.append(s);print(p.name,s['end_time_s'],flush=True)
    csv_write(OUT/'independent_size_summary.csv',allbins);csv_write(OUT/'mixture_summary_reference_scenario.csv',summaries)
    (OUT/'summary.json').write_text(json.dumps({'emission_scenario':config(),'cases':summaries,'size_cases':allbins},ensure_ascii=False,indent=2))
