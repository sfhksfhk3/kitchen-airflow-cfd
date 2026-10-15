"""Independently verify exported meshes/fields, fluxes and file consistency."""
from pathlib import Path
import json,numpy as np,meshio
ROOT=Path(__file__).resolve().parent
config=json.loads((ROOT/'config.json').read_text())
checks={}
for case in config['cases']+(['both_open_refined'] if (ROOT/'results/both_open_refined.npz').exists() else []):
 r=ROOT/'results';d=np.load(r/(case+'.npz'));metric=json.loads((r/(case+'_metrics.json')).read_text())
 p,t,U,T=d['points'],d['cells'],d['velocity'],d['temperature']
 E=p[t[:,1:]]-p[t[:,0],None,:];vol=np.abs(np.linalg.det(E))/6
 allfaces=np.concatenate([t[:,ix] for ix in [[0,1,2],[0,1,3],[0,2,3],[1,2,3]]])
 owners=np.tile(np.arange(len(t)),4)
 _,indices,counts=np.unique(np.sort(allfaces,axis=1),axis=0,return_index=True,return_counts=True)
 bf=allfaces[indices[counts==1]];owner=owners[indices[counts==1]]
 vp=p[bf];cf=vp.mean(axis=1);normal=np.cross(vp[:,1]-vp[:,0],vp[:,2]-vp[:,0])/2
 outward=cf-p[t[owner]].mean(axis=1)
 normal[np.sum(normal*outward,axis=1)<0]*=-1
 faceflux=np.sum(U[bf].mean(axis=1)*normal,axis=1)
 independent_flux=float(faceflux.sum());relative_flux=abs(independent_flux)/metric['fan_flow_m3_s']
 hood=np.isclose(cf[:,2],config['hood']['z'],atol=1e-10)&(cf[:,0]>config['hood']['x0'])&(cf[:,0]<config['hood']['x1'])&(cf[:,1]<config['hood']['y1'])
 hoodflux=float(faceflux[hood].sum());hoodarea=float(np.linalg.norm(normal[hood],axis=1).sum())
 G=np.linalg.solve(E,U[t[:,1:]]-U[t[:,0],None,:]);div=np.trace(G,axis1=1,axis2=2)
 divratio=float(np.sqrt(np.sum(div**2*vol)/np.sum(G**2*vol[:,None,None])))
 csv=np.loadtxt(r/(case+'_nodes.csv'),delimiter=',',skiprows=1)
 vtu=meshio.read(r/(case+'.vtu'))
 checks[case]={'independent_fluid_volume_m3':float(vol.sum()),'minimum_tetra_volume_m3':float(vol.min()),'independent_net_boundary_flux_m3_s':independent_flux,'independent_relative_flux_error':relative_flux,'independent_hood_flow_m3_s':hoodflux,'independent_hood_area_m2':hoodarea,'independent_divergence_gradient_ratio':divratio,'CSV_matches_NPZ':bool(np.allclose(csv,np.c_[p,U,d['pressure'],T],rtol=0,atol=1e-12)),'VTU_mesh_and_velocity_match':bool(np.array_equal(vtu.cells_dict['tetra'],t) and np.allclose(vtu.points,p) and np.allclose(vtu.point_data['velocity_m_s'],U)),'direct_heat_balance_relative_error':abs(metric['discrete_heat_balance_residual_W'])/metric['heat_input_W'],'minimum_temperature_C':float(T.min()),'converged':metric['converged']}
 c=checks[case]
 assert np.isfinite(U).all() and np.isfinite(T).all() and vol.min()>0
 assert c['CSV_matches_NPZ'] and c['VTU_mesh_and_velocity_match']
 assert abs(vol.sum()-metric['fluid_volume_m3'])<1e-9
 assert abs(divratio-metric['divergence_gradient_norm_ratio'])<1e-9
 assert relative_flux<1e-9 and abs(hoodflux-metric['fan_flow_m3_s'])<1e-9
 assert c['direct_heat_balance_relative_error']<1e-4 and T.min()>=config['physics']['ambient_C']-1e-8
 assert metric['converged'] and metric['nonlinear_relative_residual']<1e-3 and metric['thermal_relative_residual']<1e-3
(ROOT/'results/validation.json').write_text(json.dumps({'method':'Independent tetrahedron volumes and boundary-face flux integration, gradients reconstructed from exported nodal values, CSV/VTU/NPZ consistency; no solver assembly reused.','cases':checks},ensure_ascii=False,indent=2))
print(json.dumps(checks,ensure_ascii=False,indent=2))
