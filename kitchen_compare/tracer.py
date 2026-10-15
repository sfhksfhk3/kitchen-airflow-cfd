"""Conservative transient P1 FEM gas-like tracer for layout/door comparisons.
Danckwerts clean-air inflow, natural advective outflow, impermeable walls.
M-lumped backward Euler + symmetric algebraic edge diffusion; no deposition.
"""
from pathlib import Path
import argparse,json,csv,time
import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import splu,spsolve
from skfem import MeshTet,Basis,FacetBasis,ElementTetP1,BilinearForm,LinearForm,asm
from skfem.helpers import dot,grad
from conservative_flux import reconstruct
ROOT=Path(__file__).resolve().parent
@BilinearForm
def adv_diff(u,v,w):return w.diff*dot(grad(u),grad(v))-u*dot(w.vel,grad(v))
@BilinearForm
def outward_adv(u,v,w):return w.vn*u*v
@LinearForm
def load(v,w):return w.source*v
@LinearForm
def weights(v,w):return v
@LinearForm
def outflux_weight(v,w):return w.vn*v

def graph_diffusion(A):
 off=A.tocsr().copy();off.setdiag(0);off.eliminate_zeros()
 pos=off.maximum(0);edge=pos.maximum(pos.T)
 return (A-edge+diags(np.asarray(edge.sum(axis=1)).ravel())).tocsc()

def classify(mesh,c,case):
 f=mesh.boundary_facets();p=mesh.p[:,mesh.facets[:,f]].mean(axis=1);d=c['doors'];h=c['hood'];L=c['room']['length']
 ishood=np.isclose(p[2],h['z'],atol=1e-8)&(p[0]>h['x0']+1e-8)&(p[0]<h['x1']-1e-8)&(p[1]>h['y0']+1e-8)&(p[1]<h['y1']-1e-8)
 left=np.isclose(p[0],0,atol=1e-8)&(p[1]>d['y0'])&(p[1]<d['y1']);right=np.isclose(p[0],L,atol=1e-8)&(p[1]>d['y0'])&(p[1]<d['y1'])
 if case=='both_closed':left&=p[2]<d['homogenized_inlet_height_m'];right&=p[2]<d['homogenized_inlet_height_m']
 if case=='left_open':right[:]=False
 if case=='right_open':left[:]=False
 return {k:f[v] for k,v in [('hood',ishood),('left',left),('right',right)]}

def simulate(layout,case,field_suffix='',Sc=None,dt=None,suffix=''):
 folder=ROOT/layout;r=folder/'results';c=json.loads((folder/'config.json').read_text());tr=c['tracer'];Sc=tr['Sc_t'] if Sc is None else Sc;dt=tr['dt_s'] if dt is None else dt
 name=case+field_suffix;metric=json.loads((r/(name+'_metrics.json')).read_text())
 if not metric['converged']:raise ValueError('Tracer requires converged airflow')
 d=dict(np.load(r/(name+'.npz')));m=MeshTet(d['points'].T,d['cells'].T);bs=Basis(m,ElementTetP1(),intorder=5);n=bs.N;U=d['velocity'].T
 boundary=classify(m,c,case)
 urt,projected_flux,face_area,projection_meta=reconstruct(m,U,np.r_[boundary['left'],boundary['right']],boundary['hood'])
 vel=np.broadcast_to(urt.T[:,:,None],(3,m.nelements,bs.dx.shape[1])).copy();D=tr['molecular_diffusivity_m2_s']+metric['nu_t_m2_s']/Sc
 A=asm(adv_diff,bs,vel=vel,diff=D);fluxweights={}
 for key,fac in boundary.items():
  if len(fac)==0:fluxweights[key]=np.zeros(n);continue
  fb=FacetBasis(m,ElementTetP1(),facets=fac,intorder=5)
  area_vectors=np.cross(m.p[:,m.facets[1,fac]].T-m.p[:,m.facets[0,fac]].T,m.p[:,m.facets[2,fac]].T-m.p[:,m.facets[0,fac]].T)/2
  sign=np.sign(np.sum(area_vectors*fb.normals[:,:,0].T,axis=1))
  vn=np.broadcast_to(np.maximum(sign*projected_flux[fac]/face_area[fac],0)[:,None],fb.dx.shape).copy()
  A+=asm(outward_adv,fb,vn=vn);fluxweights[key]=asm(outflux_weight,fb,vn=vn)
 A=graph_diffusion(A);M=asm(weights,bs)
 x=bs.global_coordinates();s=np.zeros_like(x[0]);source_share=tr['total_release_reference_units_s']/len(c['burners'])
 for bx,by,bz in c['burners']:
  q=np.exp(-((x[0]-bx)**2+(x[1]-by)**2)/(2*.085**2)-((x[2]-bz-.07)**2)/(2*.055**2))
  q*=((x[2]>bz)&(x[2]<bz+.25)&(x[1]<c['counter']['depth']))
  q*=source_share/np.sum(q*bs.dx);s+=q
 S=asm(load,bs,source=s);R=float(S.sum());assert abs(R-1)<1e-10 and S.min()>-1e-13
 # Column sum equals true boundary outflux, since added edge diffusion is conservative.
 outweights=sum(fluxweights.values());column_error=float(np.max(np.abs(np.asarray(A.sum(axis=0)).ravel()-outweights)))
 if column_error>1e-9:raise RuntimeError('Conservative tracer column-sum check failed')
 Aoff=A.copy();Aoff.setdiag(0);assert Aoff.data.max(initial=0)<1e-10
 css=spsolve(A,S);assert css.min()>-1e-10;css=np.maximum(css,0)
 coords=m.p.T;k=c['counter'];node_occ=(coords[:,1]>k['depth'])&(coords[:,2]>k['height'])&(coords[:,2]<k['height']+.9)
 occq=(x[1]>k['depth'])&(x[2]>k['height'])&(x[2]<k['height']+.9);occweights=asm(load,bs,source=occq.astype(float));occvol=float(occweights.sum())
 # Above-counter apron outside hood projection, more directly related to spill.
 spillq=(x[1]>.40)&(x[1]<.90)&(x[2]>.90)&(x[2]<1.60);spillweights=asm(load,bs,source=spillq.astype(float));spillvol=float(spillweights.sum())
 cook=float(tr['cooking_duration_s']);end=cook+float(tr['clearance_duration_s']);steps=int(round(end/dt))
 Astep=A+diags(M/dt);lu=splu(Astep);C=np.zeros(n);cum={'hood':0.,'left':0.,'right':0.};inventory_stop=None;t90=None;t95=None;t99=None;history=[];frames=[];frame_times=[]
 wanted={0.,5.,10.,20.,30.,60.,120.,121.,122.,124.,126.,128.,130.,135.,140.,150.,180.,240.,300.,420.,600.,1020.}
 for i in range(steps+1):
  t=round(i*dt,8)
  if i:
   on=t<=cook+1e-8;C=lu.solve(M*C/dt+(S if on else 0))
   if C.min()<-1e-10:raise RuntimeError('Tracer negative concentration')
   for key in cum:cum[key]+=dt*float(fluxweights[key]@C)
  inventory=float(M@C);released=R*min(t,cook);balance=inventory+sum(cum.values())-released
  if abs(t-cook)<dt/2:inventory_stop=inventory
  if t>cook and inventory_stop is not None:
   frac=inventory/inventory_stop
   if t90 is None and frac<=.1:t90=t-cook
   if t95 is None and frac<=.05:t95=t-cook
   if t99 is None and frac<=.01:t99=t-cook
  if t in wanted:frames.append(C.copy());frame_times.append(t)
  history.append([t,inventory,*cum.values(),released,balance,float(occweights@C/occvol),float(spillweights@C/spillvol),float(C.max())])
 H=np.array(history);atcook=H[np.argmin(abs(H[:,0]-cook))];steady_inv=float(M@css)
 final_balance=float(np.max(np.abs(H[:,6])));ss_out={k:float(w@css) for k,w in fluxweights.items()}
 summary={'velocity_transport':projection_meta,'layout':layout,'case':case,'field_dataset':name,'label':c['case_labels'][case],'source_reference_units_s':R,'concentration_unit':'s/m³ per unit total emission rate; multiply by actual emission units/s for actual units/m³','inventory_unit':'equivalent source-seconds; divide by release_duration_s for retained fraction','cooking_duration_s':cook,'clearance_heat_condition':'Hold steady 4000W heat and fan; stop only tracer release','Sc_t':Sc,'effective_diffusivity_m2_s':D,'dt_s':dt,'steady_equivalent_residence_time_s':steady_inv/R,'retained_source_seconds_at_120s':float(atcook[1]/R),'retained_fraction_at_120s':float(atcook[1]/(R*cook)),'cumulative_hood_removed_fraction_at_120s':float(atcook[2]/(R*cook)),'cumulative_door_escape_fraction_at_120s':float((atcook[3]+atcook[4])/(R*cook)),'occupied_zone_normalized_concentration_at_120s':float(atcook[7]),'spill_zone_normalized_concentration_at_120s':float(atcook[8]),'steady_occupied_zone_normalized_concentration':float(occweights@css/occvol),'steady_spill_zone_normalized_concentration':float(spillweights@css/spillvol),'clearance_t90_s':t90,'clearance_t95_s':t95,'clearance_t99_s':t99,'mass_balance_max_reference_units':final_balance,'mass_balance_max_relative_to_emitted':final_balance/(R*cook),'steady_outflux_reference_units_s':ss_out,'matrix_conservation_column_error':column_error,'min_concentration':float(min(css.min(),np.min(frames))),'note':'Steady total hood fraction may be forced to unity if it is the only tracer outflow; use finite-time retention, occupied-zone concentration and decay instead of calling this real capture efficiency.'}
 outname=name+suffix
 np.savez_compressed(r/(outname+'_tracer.npz'),points=d['points'],cells=d['cells'],velocity=d['velocity'],temperature=d['temperature'],pressure=d['pressure'],concentration_steady=css,concentration_frames=np.array(frames),frame_times_s=np.array(frame_times),mass_lumped=M,velocity_RT0_cells=urt,projected_face_flux_m3_s=projected_flux)
 np.savetxt(r/(outname+'_tracer_history.csv'),H,delimiter=',',header='time_s,retained_source_seconds,cumulative_hood_removed_source_seconds,cumulative_left_door_out_source_seconds,cumulative_right_door_out_source_seconds,total_released_source_seconds,mass_balance_error,occupied_zone_normalized_c_s_m3,spill_zone_normalized_c_s_m3,max_normalized_c_s_m3',comments='')
 np.savetxt(r/(outname+'_tracer_nodes.csv'),np.c_[d['points'],css,frames[frame_times.index(120.)]],delimiter=',',header='x_m,y_m,z_m,normalized_c_steady_s_m3,normalized_c_120s_s_m3',comments='')
 import meshio
 meshio.write(r/(outname+'_combined.vtu'),meshio.Mesh(d['points'],[('tetra',d['cells'])],point_data={'velocity_FEM_m_s':d['velocity'],'temperature_C':d['temperature'],'pressure_room_reference_Pa':d['pressure'],'normalized_tracer_steady_s_m3':css,'normalized_tracer_120s_s_m3':frames[frame_times.index(120.)]},cell_data={'velocity_RT0_transport_m_s':[urt]}))
 for column,tag in [(7,'occupied_zone'),(8,'spill_zone')]:
  for fraction,label in [(.1,'t90'),(.05,'t95')]:
   rows=H[(H[:,0]>cook)&(H[:,column]<=atcook[column]*fraction)]
   summary['clearance_'+tag+'_'+label+'_s']=float(rows[0,0]-cook) if len(rows) else None
 (r/(outname+'_tracer_metrics.json')).write_text(json.dumps(summary,ensure_ascii=False,indent=2));print(json.dumps(summary,ensure_ascii=False),flush=True)
 return summary
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--layout',choices=['corner','center'],required=True);ap.add_argument('--case',default='all');ap.add_argument('--field-suffix',default='');ap.add_argument('--Sc',type=float);ap.add_argument('--dt',type=float);ap.add_argument('--suffix',default='');a=ap.parse_args()
 c=json.loads((ROOT/a.layout/'config.json').read_text())
 for case in c['cases'] if a.case=='all' else [a.case]:simulate(a.layout,case,a.field_suffix,a.Sc,a.dt,a.suffix)
