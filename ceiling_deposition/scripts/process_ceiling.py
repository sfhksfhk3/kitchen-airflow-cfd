"""Per-face deposition, conservative budgets, normalized and reference results."""
from pathlib import Path
import json,sys,numpy as np,pandas as pd,xml.etree.ElementTree as ET
from scipy.special import ndtr
from foam_mesh import mesh_geometry
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'results';OUT.mkdir(exist_ok=True)
DS=np.array([.1,.3,1,3,10]);REF=json.loads((ROOT/'config/emission_scenario.json').read_text());WEIGHTS=np.array(REF['mass_weights'])

def size_weights(gsd=1.5):
 edges=np.r_[0,np.sqrt(DS[:-1]*DS[1:]),np.inf]
 with np.errstate(divide='ignore'):p=ndtr(np.log(edges/2.5)/np.log(gsd))
 return np.diff(p)
def weighted_quantile(x,weights,q=.95):
 order=np.argsort(x);x=np.array(x)[order];a=np.array(weights)[order]
 return float(x[np.searchsorted(np.cumsum(a),q*a.sum(),side='left')])
def write_vtp(path,mesh,ids,arrays):
 used=np.unique(np.concatenate([mesh['faces'][i] for i in ids]));index={int(i):k for k,i in enumerate(used)}
 polys=[[index[int(i)] for i in mesh['faces'][j]] for j in ids];connect=[i for p in polys for i in p];off=np.cumsum([len(p) for p in polys])
 root=ET.Element('VTKFile',type='PolyData',version='0.1',byte_order='LittleEndian');pdv=ET.SubElement(root,'PolyData');piece=ET.SubElement(pdv,'Piece',NumberOfPoints=str(len(used)),NumberOfPolys=str(len(ids)),NumberOfVerts='0',NumberOfLines='0',NumberOfStrips='0')
 def array(parent,values,name=None,typ='Float64',nc=None):
  kw={'type':typ,'format':'ascii'}
  if name:kw['Name']=name
  if nc:kw['NumberOfComponents']=str(nc)
  ET.SubElement(parent,'DataArray',**kw).text=' '.join(str(x) if typ=='Int64' else f'{x:.14g}' for x in np.asarray(values).ravel())
 points=ET.SubElement(piece,'Points');array(points,mesh['points'][used],nc=3)
 pol=ET.SubElement(piece,'Polys');array(pol,connect,'connectivity','Int64');array(pol,off,'offsets','Int64')
 cd=ET.SubElement(piece,'CellData',Scalars='M_unit_mg_m2_per_mg')
 for n,values in arrays.items():array(cd,values,n,'Int64' if n=='original_face_id' else 'Float64')
 ET.ElementTree(root).write(path,encoding='utf-8',xml_declaration=True)

def process(case):
 case=Path(case);name=case.name;manifest=json.loads((case/'replay_manifest.json').read_text());out=OUT/name;out.mkdir(exist_ok=True)
 d=pd.read_csv(case/'diagnostics/ceiling_budget.csv');last=d.groupby('diameter_um').tail(1).set_index('diameter_um').loc[DS]
 source=last.released_kg.to_numpy();initial=last.initial_kg.to_numpy();den=source+initial;total_in=float(WEIGHTS@den)
 assert (den>0).all();t_end=float(last.time_s.min());refscale=REF['total_rate_ug_s']*1e-9/manifest['parameters']['emissionKgPerSecPerBin'] # original independent rate 1e-9 kg/s
 geometry=mesh_geometry(case);ids=geometry['ceiling_ids'];area=geometry['area'][ids];x=geometry['centres'][ids]
 snaps=sorted(case.joinpath('diagnostics').glob('ceilingFaces_*.csv'),key=lambda p:float(p.stem.split('_',1)[1]));records=[];timeline=[];previous=np.zeros(len(ids));cookingM=None
 for p in snaps:
  tt=float(p.stem.split('_',1)[1]);raw=pd.read_csv(p);matM=raw.pivot(index='face_id',columns='diameter_um',values='M_kg_m2').loc[ids,DS].to_numpy();matJ=raw.pivot(index='face_id',columns='diameter_um',values='J_kg_m2_s').loc[ids,DS].to_numpy()
  # denominator is total input of the full observation, fixed across the PVD series.
  unitM=(matM@WEIGHTS)/total_in;unitJ=(matJ@WEIGHTS)/total_in;refM=(matM@WEIGHTS)*refscale*1e6;refJ=(matJ@WEIGHTS)*refscale
  arrays={'original_face_id':ids,'area_m2':area,'M_unit_mg_m2_per_mg':unitM,'M_reference_mg_m2':refM,'J_reference_kg_m2_s':refJ,'J_unit_mg_m2_s_per_mg':unitJ}
  for i,ds in enumerate(DS):arrays[f'M_{ds:g}um_unit_mg_m2_per_mg']=matM[:,i]/den[i];arrays[f'J_{ds:g}um_kg_m2_s_reference']=matJ[:,i]*refscale
  file=f'ceiling_{tt:g}.vtp';write_vtp(out/file,geometry,ids,arrays);timeline.append((tt,file))
  records.append(dict(time_s=tt,total_ceiling_unit_mg_per_mg=float(unitM@area),total_ceiling_reference_mg=float(refM@area)))
  assert np.min(unitM-previous)>-1e-18;previous=unitM
  if abs(tt-120)<1e-5:cookingM=matM.copy()
 # Verify per-face integrals equal the independently accumulated step budgets.
 integral=area@matM;expected=last.ceiling_kg.to_numpy();face_err=float(np.max(abs(integral-expected)/np.maximum(abs(expected),1e-300)))
 assert np.allclose(integral,expected,rtol=2e-10,atol=1e-25),(case,face_err)
 pvd=ET.Element('VTKFile',type='Collection',version='0.1',byte_order='LittleEndian');col=ET.SubElement(pvd,'Collection')
 for tt,file in timeline:ET.SubElement(col,'DataSet',timestep=str(tt),group='',part='0',file=file)
 ET.ElementTree(pvd).write(out/'ceiling_deposition.pvd',encoding='utf-8',xml_declaration=True)
 face=pd.DataFrame({'original_face_id':ids,'x_m':x[:,0],'y_m':x[:,1],'z_m':x[:,2],'area_m2':area,'M_unit_mg_m2_per_mg':unitM,'M_reference_mg_m2':refM,'J_reference_kg_m2_s':refJ,'J_mean_reference_kg_m2_s':refM*1e-6/(t_end-manifest['start_s'])})
 for i,ds in enumerate(DS):face[f'M_{ds:g}um_unit_mg_m2_per_mg']=matM[:,i]/den[i];face[f'J_{ds:g}um_reference_kg_m2_s']=matJ[:,i]*refscale
 face.to_csv(out/'ceiling_faces.csv',index=False);pd.DataFrame(records).to_csv(out/'ceiling_history.csv',index=False)
 bins=[]
 for i,ds in enumerate(DS):
  row=last.loc[ds];r=dict(diameter_um=float(ds),input_kg=float(den[i]),emitted_kg=float(source[i]),initial_kg=float(initial[i]),ceiling_fraction=float(expected[i]/den[i]),other_wall_fraction=float(row.other_wall_kg/den[i]),hood_fraction=float(row.hood_out_kg/den[i]),left_fraction=float(row.left_out_kg/den[i]),right_fraction=float(row.right_out_kg/den[i]),airborne_fraction=float(row.stock_kg/den[i]),mass_error_fraction=float(row.mass_error_kg/den[i]));bins.append(r)
 pd.DataFrame(bins).to_csv(out/'size_mass_balance.csv',index=False)
 # Mixture budget must weight raw kg, not ratios if different initial size stocks.
 total_in=float(WEIGHTS@den);mix={}
 for key,col in [('ceiling','ceiling_kg'),('other_wall','other_wall_kg'),('hood','hood_out_kg'),('left','left_out_kg'),('right','right_out_kg'),('airborne','stock_kg')]:mix[key+'_fraction']=float(WEIGHTS@last[col].to_numpy()/total_in)
 j=int(np.argmax(unitM));hot=weighted_quantile(unitM,area,.95);hot_ids=np.flatnonzero(unitM>=hot)
 cookingTotal=float((area@cookingM)@WEIGHTS*refscale*1e6) if cookingM is not None else None
 q=d[d.time_s>=min(120,manifest['end_s']-5)].Qhood_m3_s
 summary=dict(case=name,start_s=manifest['start_s'],end_s=t_end,ceiling_area_m2=float(area.sum()),ceiling_face_count=len(ids),ceiling_fraction=mix['ceiling_fraction'],ceiling_percent=100*mix['ceiling_fraction'],M_total_unit_mg_per_mg=float(unitM@area),M_total_reference_mg=float(refM@area),M_mean_unit_mg_m2_per_mg=float(unitM@area/area.sum()),M_p95_unit_mg_m2_per_mg=weighted_quantile(unitM,area),M_max_unit_mg_m2_per_mg=float(unitM[j]),M_mean_reference_mg_m2=float(refM@area/area.sum()),M_p95_reference_mg_m2=weighted_quantile(refM,area),M_max_reference_mg_m2=float(refM.max()),max_face_id=int(ids[j]),max_x_m=float(x[j,0]),max_y_m=float(x[j,1]),high_5pct_area_centroid_xy_m=np.average(x[hot_ids,:2],weights=area[hot_ids],axis=0).tolist(),deposited_by_source_stop_reference_mg=cookingTotal,post_stop_deposition_reference_mg=float(refM@area)-cookingTotal if cookingTotal is not None else None,max_relative_mass_error=float((d.mass_error_kg.abs()/(d.initial_kg+d.released_kg).clip(lower=1e-30)).max()),per_face_budget_relative_error=face_err,air_continuity_L1_kg_s_p95=float(d.air_continuity_L1_kg_s.quantile(.95)),max_cell_Stokes=float(d.max_cell_Stokes.max()),ceiling_yplus_max=float(d.ceiling_yplus_max.max()),ceiling_delta_max_m=float(d.ceiling_delta_max_m.max()),Q_late_mean_m3_s=float(q.mean()),input_reference_mg=total_in*refscale*1e6,**{k:v for k,v in mix.items() if k!='ceiling_fraction'})
 (out/'summary.json').write_text(json.dumps(summary,indent=2));return summary

def main():
 cases=[p for p in sorted((ROOT/'runs').iterdir()) if (p/'RUN_COMPLETE.json').exists()];results=[]
 for p in cases:
  if p.name.startswith('pilot'):continue
  cached=OUT/p.name/'summary.json'
  if cached.exists() and cached.stat().st_mtime >= (p/'RUN_COMPLETE.json').stat().st_mtime:
   results.append(json.loads(cached.read_text()));continue
  results.append(process(p));print(p.name,'processed',flush=True)
 pd.DataFrame(results).to_csv(OUT/'case_summary.csv',index=False)
 (OUT/'all_summaries.json').write_text(json.dumps(results,ensure_ascii=False,indent=2))
if __name__=='__main__':main()
