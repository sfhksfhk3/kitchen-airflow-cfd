"""Independent verification of exported CFD / tracer data and geometry."""
from pathlib import Path
import json,numpy as np,meshio
from skfem import MeshTet
ROOT=Path(__file__).resolve().parent
LAYOUTS=['corner','center'];CASES=['both_open','both_closed','left_open','right_open'];checks={}
for layout in LAYOUTS:
 cfg=json.loads((ROOT/layout/'config.json').read_text());r=ROOT/layout/'results'
 for case in CASES:
  for suffix in ['', '_refined']:
   name=case+suffix;flow=json.loads((r/(name+'_metrics.json')).read_text());smoke=json.loads((r/(name+'_tracer_metrics.json')).read_text());d=dict(np.load(r/(name+'_tracer.npz')))
   p,t=d['points'],d['cells'];U=d['velocity'];UR=d['velocity_RT0_cells'];E=p[t[:,1:]]-p[t[:,0],None,:];volume=np.abs(np.linalg.det(E))/6
   m=MeshTet(p.T,t.T);f=m.facets.T;fv=p[f];normal=np.cross(fv[:,1]-fv[:,0],fv[:,2]-fv[:,0])/2;cf=fv.mean(axis=1);owner=m.f2t[0];outward=cf-p[t[owner]].mean(axis=1);sgn=np.sign(np.sum(normal*outward,axis=1));normal_out=normal*sgn[:,None]
   boundary=m.f2t[1]<0;hood=boundary&np.isclose(cf[:,2],cfg['hood']['z'],atol=1e-8)&(cf[:,0]>cfg['hood']['x0'])&(cf[:,0]<cfg['hood']['x1'])&(cf[:,1]>0)&(cf[:,1]<cfg['hood']['y1'])
   flux=np.sum(UR[owner]*normal_out,axis=1);normaljump=np.sum((UR[m.f2t[0,~boundary]]-UR[m.f2t[1,~boundary]])*normal[~boundary],axis=1)
   c=d['concentration_steady'];mass_weights=d['mass_lumped'];times=d['frame_times_s'];C120=d['concentration_frames'][list(times).index(120.)]
   H=np.loadtxt(r/(name+'_tracer_history.csv'),delimiter=',',skiprows=1)
   net_flux=float(flux[boundary].sum());hood_flux=float(flux[hood].sum());tracer_out=float(np.sum(np.maximum(flux[boundary],0)*c[f[boundary]].mean(axis=1)))
   csvdata=np.loadtxt(r/(name+'_tracer_nodes.csv'),delimiter=',',skiprows=1)
   nodal_consistency=bool(np.allclose(csvdata,np.c_[p,c,C120],rtol=0,atol=1e-12))
   vtu_ok=True
   combined=r/(name+'_combined.vtu')
   if combined.exists():
    vtu=meshio.read(combined);vtu_ok=bool(np.allclose(vtu.point_data['normalized_tracer_120s_s_m3'],C120) and np.allclose(vtu.cell_data_dict['velocity_RT0_transport_m_s']['tetra'],UR))
   elif suffix=='_refined':raise ValueError('Missing combined primary VTU')
   key=layout+'_'+name
   checks[key]={'mesh_nodes':len(p),'tetrahedra':len(t),'independent_fluid_volume_m3':float(volume.sum()),'geometry_approx_volume_error_m3':float(volume.sum()-9.18282),'minimum_tetra_volume_m3':float(volume.min()),'RT0_internal_normal_flux_jump_max_m3_s':float(np.max(abs(normaljump))),'RT0_independent_net_boundary_flux_m3_s':net_flux,'RT0_independent_hood_flow_m3_s':hood_flux,'independent_steady_tracer_total_outflow_reference_units_s':tracer_out,'independent_steady_inventory_source_seconds':float(mass_weights@c),'tracer_transient_mass_balance_max_reference_units':float(np.max(abs(H[:,1]+H[:,2]+H[:,3]+H[:,4]-H[:,5]))),'direct_heat_balance_residual_W':flow['discrete_heat_balance_residual_W'],'FEM_nonlinear_relative_residual':flow['nonlinear_relative_residual'],'FEM_thermal_relative_residual':flow['thermal_relative_residual'],'FEM_divergence_gradient_ratio':flow['divergence_gradient_norm_ratio'],'RT0_cell_velocity_relative_RMS_correction':smoke['velocity_transport']['relative_RMS_cell_velocity_correction'],'CSV_matches_NPZ':nodal_consistency,'combined_VTU_matches_NPZ':vtu_ok,'min_tracer':float(min(c.min(),d['concentration_frames'].min())),'min_temperature_C':float(d['temperature'].min()),'max_density_reduction_buoyancy_model':(flow['max_nodal_temperature_C']-25)/(flow['max_nodal_temperature_C']+273.15),'converged':flow['converged']}
   assert flow['converged'] and flow['heat_input_W']==4000 and flow['nonlinear_relative_residual']<1e-3 and flow['thermal_relative_residual']<1e-3
   assert nodal_consistency and vtu_ok and volume.min()>0 and abs(volume.sum()-9.18282)<2e-6
   assert abs(net_flux)<1e-8 and abs(hood_flux-flow['fan_flow_m3_s'])<1e-8 and np.max(abs(normaljump))<1e-8
   assert abs(tracer_out-1)<1e-7 and abs(mass_weights@c-smoke['steady_equivalent_residence_time_s'])<1e-8
   assert abs(mass_weights@C120-smoke['retained_fraction_at_120s']*120)<1e-8
   assert checks[key]['tracer_transient_mass_balance_max_reference_units']<1e-7
   assert abs(flow['discrete_heat_balance_residual_W'])<.01 and d['temperature'].min()>=25-1e-7
   assert checks[key]['min_tracer']>=-1e-10 and np.isfinite(c).all()
(ROOT/'validation.json').write_text(json.dumps({'method':'Independent exported tetrahedron volumes, RT0 shared-face normal continuity, boundary volume/tracer flux, inventory and time-history mass balance, and CSV/NPZ/VTU consistency. No CFD or transport matrix reused.','checks':checks},ensure_ascii=False,indent=2))
print('Validated',len(checks),'CFD/tracer datasets')
