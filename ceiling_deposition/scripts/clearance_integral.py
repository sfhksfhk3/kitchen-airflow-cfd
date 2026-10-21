"""Integrate remaining deposition on frozen terminal air flow.

For autonomous linear transport M'=-A M, the integrated residence satisfies
A F=M0. A large-pseudo-step backward Euler computes (A+I/H)^-1 M0; H is
an algebraic regularizer, NOT a cooking duration. Account for its remaining
inventory explicitly. Never label this field as a physical future timestamp.
"""
from pathlib import Path
import subprocess,json,time,shutil,re,numpy as np,pandas as pd
from prepare_replay import prepare,ROOT,NAMES
from process_ceiling import OUT,WEIGHTS,REF,weighted_quantile,write_vtp
from foam_mesh import mesh_geometry

def complete(case_name,H=1e12):
 base=ROOT/'runs'/case_name
 while not (base/'RUN_COMPLETE.json').exists():time.sleep(5)
 marker=ROOT/'verification'/f'{case_name}_clearance_{H:g}.json'
 if marker.exists():return json.loads(marker.read_text())
 manifest=json.loads((base/'replay_manifest.json').read_text());source=Path(manifest['source']);t0=manifest['end_s'];end=t0+H
 src=ROOT/'tail_sources'/case_name;src.mkdir(parents=True,exist_ok=True)
 if not (src/f'{t0:g}').exists():
  shutil.copytree(base/'constant',src/'constant');shutil.copy2(base/'mesh_reference.npz',src/'mesh_reference.npz');fr=src/f'{t0:g}';fr.mkdir()
  original=min((p for p in source.iterdir() if p.is_dir() and re.fullmatch(r'\d+(?:\.\d+)?',p.name)),key=lambda p:abs(float(p.name)-t0))
  finished=min((p for p in base.iterdir() if p.is_dir() and re.fullmatch(r'\d+(?:\.\d+)?',p.name)),key=lambda p:abs(float(p.name)-t0))
  for field in ['U','T','p','k','nut','phi']:(fr/field).symlink_to(original/field)
  for field in NAMES:(fr/field).symlink_to(finished/field)
  (fr/'sourceShape').symlink_to(base/f"{manifest['start_s']:.12g}"/'sourceShape')
 if not (src/'system').exists():shutil.copytree(base/'system',src/'system')
 alias=src/f'{end:.0f}'
 if not alias.exists():alias.symlink_to(src/f'{t0:g}',target_is_directory=True)
 name=case_name+'_clearance_'+f'{H:g}'.replace('+','');case=ROOT/'runs'/name
 if case.exists() and not (case/'replay_manifest.json').exists():shutil.rmtree(case)
 if not case.exists():
  opts=manifest['parameters'].copy();opts.update(useSavedInitialAerosols=True,sourceStopTime=120,coefficientUpdateSeconds=H*2)
  prepare(src,name,start=t0,end=end,dt=H,**opts)
  p=case/'system/fvSolution';p.write_text(p.read_text().replace('1e-20','1e-30').replace('relTol 1e-10','relTol 1e-12'))
 p=case/'system/controlDict';text=p.read_text();text=re.sub(r'writeControl\s+[^;]+;', 'writeControl timeStep;',text);text=re.sub(r'writeInterval\s+[^;]+;', 'writeInterval 1;',text);p.write_text(text)
 log=ROOT/'verification'/f'{name}.log'
 if not log.exists() or 'End replay' not in log.read_text()[-1000:]:
  with log.open('w') as f:subprocess.run(['docker','exec','kitchen-openfoam-2406','bash','-lc',f'{ROOT}/bin/ceilingTailFoam -case {case}'],stdout=f,stderr=subprocess.STDOUT,check=True)
 d=pd.read_csv(case/'diagnostics/ceiling_budget.csv');err=float((d.mass_error_kg.abs()/d.initial_kg.clip(lower=1e-300)).max());res=float((d.stock_kg/d.initial_kg.clip(lower=1e-300)).max())
 assert err<1e-7,(name,err)
 result=dict(case=case_name,algebraic_H_s=H,method='frozen-terminal-flow linear clearance resolvent; H is not a physical observation time',remaining_fraction_of_tail_input_max=res,max_relative_mass_error=err,tail_case=str(case))
 marker.write_text(json.dumps(result,indent=2));print(json.dumps(result),flush=True);return result

def merge(result):
 name=result['case'];out=OUT/name;out.mkdir(exist_ok=True);base=ROOT/'runs'/name;tail=Path(result['tail_case']);m=json.loads((base/'replay_manifest.json').read_text());end=m['end_s']
 def read(case,t=None):
  pp=list((case/'diagnostics').glob('ceilingFaces_*.csv'));p=max(pp,key=lambda p:float(p.stem.split('_',1)[1]));return pd.read_csv(p)
 a=read(base).set_index(['face_id','diameter_um']);b=read(tail).set_index(['face_id','diameter_um']);assert a.index.equals(b.index)
 raw=a.reset_index();raw['M_observed_kg_m2']=raw.M_kg_m2;raw['M_tail_kg_m2']=b.M_kg_m2.to_numpy();raw['M_kg_m2']+=raw.M_tail_kg_m2
 raw.to_csv(out/'ceiling_complete_event_by_size.csv',index=False)
 budget=pd.read_csv(base/'diagnostics/ceiling_budget.csv').groupby('diameter_um').tail(1).sort_values('diameter_um');tailb=pd.read_csv(tail/'diagnostics/ceiling_budget.csv').sort_values('diameter_um');E=float(WEIGHTS@budget.released_kg.to_numpy());scale=REF['total_rate_ug_s']*1e-9/m['parameters']['emissionKgPerSecPerBin']
 geom=mesh_geometry(base);ids=geom['ceiling_ids'];area=geom['area'][ids];xyz=geom['centres'][ids];matrix=raw.pivot(index='face_id',columns='diameter_um',values='M_kg_m2').loc[ids].to_numpy();obs=raw.pivot(index='face_id',columns='diameter_um',values='M_observed_kg_m2').loc[ids].to_numpy();unit=(matrix@WEIGHTS)/E;ref=(matrix@WEIGHTS)*scale*1e6
 f=pd.DataFrame({'original_face_id':ids,'x_m':xyz[:,0],'y_m':xyz[:,1],'z_m':xyz[:,2],'area_m2':area,'M_unit_mg_m2_per_mg':unit,'M_reference_mg_m2':ref,'M_observed_reference_mg_m2':(obs@WEIGHTS)*scale*1e6,'M_tail_reference_mg_m2':((matrix-obs)@WEIGHTS)*scale*1e6})
 for i,ds in enumerate([.1,.3,1,3,10]):f[f'M_{ds:g}um_unit_mg_m2_per_mg']=matrix[:,i]/budget.released_kg.iloc[i]
 f.to_csv(out/'ceiling_complete_event_faces.csv',index=False)
 write_vtp(out/'ceiling_complete_event.vtp',geom,ids,{'original_face_id':ids,'area_m2':area,'M_unit_mg_m2_per_mg':unit,'M_reference_mg_m2':ref,'M_observed_reference_mg_m2':f.M_observed_reference_mg_m2.to_numpy(),'M_tail_reference_mg_m2':f.M_tail_reference_mg_m2.to_numpy()})
 j=int(np.argmax(unit));s=dict(case=name,scope='0–240 s transient deposition + frozen-240 s flow integrated residual clearance',M_total_unit_mg_per_mg=float(unit@area),M_total_reference_mg=float(ref@area),M_mean_unit_mg_m2_per_mg=float(unit@area/area.sum()),M_p95_unit_mg_m2_per_mg=weighted_quantile(unit,area),M_max_unit_mg_m2_per_mg=float(unit[j]),M_mean_reference_mg_m2=float(ref@area/area.sum()),M_p95_reference_mg_m2=weighted_quantile(ref,area),M_max_reference_mg_m2=float(ref[j]),max_x_m=float(xyz[j,0]),max_y_m=float(xyz[j,1]),ceiling_percent=float(unit@area*100),tail_reference_mg=float(f.M_tail_reference_mg_m2@area),tail_deposition_share=float(f.M_tail_reference_mg_m2@area/(ref@area)) if ref@area else 0,clearance_residual_fraction=result['remaining_fraction_of_tail_input_max'])
 totalinput=WEIGHTS@(budget.released_kg.to_numpy()+budget.initial_kg.to_numpy())
 for key,col in [('ceiling','ceiling_kg'),('other_wall','other_wall_kg'),('hood','hood_out_kg'),('left','left_out_kg'),('right','right_out_kg')]:s[key+'_fraction']=float(WEIGHTS@(budget[col].to_numpy()+tailb[col].to_numpy())/totalinput)
 s['airborne_fraction']=float(WEIGHTS@tailb.stock_kg.to_numpy()/totalinput)
 (out/'complete_event_summary.json').write_text(json.dumps(s,indent=2))
 return s
if __name__=='__main__':
 cases=[f'{l}_{d}_parallel' for l in ['corner','center'] for d in ['both_open','both_closed','left_open','right_open']]+[f'{l}_closed_stick{a}' for l in ['corner','center'] for a in ['0.1','0.5']]
 for name in cases:merge(complete(name))
 a=complete('corner_both_closed_parallel',1e9)
 print('resolvent_regularization_sensitivity',json.dumps(a),flush=True)
