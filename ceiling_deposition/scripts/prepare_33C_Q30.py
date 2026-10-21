from pathlib import Path
import shutil,sys,re,json
import numpy as np
R=Path(__file__).resolve().parents[1];sys.path.insert(0,str(R/'scripts'))
from foam_mesh import mesh_geometry,list_body
S=R/'scenarios/ceiling33_Q30';S.mkdir(parents=True,exist_ok=True)
# Only reorder external faces to expose the ceiling patch. All cells, internal
# faces, vertices, neighbour addressing and physical geometry remain unchanged.
def prepare(layout,doors,fine=False,end=240):
 name=f'{layout}_{doors}_T33_Q30'+('_fine' if fine else '')
 out=R/'flow_sensitivity'/name
 if (out/'boundary_update.json').exists():return out
 if out.exists():shutil.rmtree(out)
 source=Path('/workspace/kitchen_openfoam')/('sensitivity' if fine else 'cases')/(f'{layout}_both_closed_grid' if fine else f'{layout}_{doors}')
 out.mkdir(parents=True)
 for d in ['constant','system']:shutil.copytree(source/d,out/d)
 shutil.copytree(Path('/workspace/kitchen_openfoam/cases')/f'{layout}_{doors}'/'0',out/'0')
 shutil.copy2(source/'mesh_reference.npz',out/'mesh_reference.npz')
 if (source/'case_manifest.json').exists():shutil.copy2(source/'case_manifest.json',out/'case_manifest.json')
 g=mesh_geometry(out);ids=set(g['ceiling_ids']);n,s=list_body(out/'constant/polyMesh/neighbour');perm=np.r_[np.arange(n),[i for i in range(n,len(g['faces'])) if i not in ids],sorted(ids)].astype(int)
 pm=out/'constant/polyMesh'
 def head(obj,cls):return f'FoamFile {{ version 2.0; format ascii; class {cls}; object {obj}; }}\n'
 (pm/'faces').write_text(head('faces','faceList')+str(len(perm))+'\n(\n'+'\n'.join(str(len(g['faces'][i]))+'('+' '.join(map(str,g['faces'][i]))+')' for i in perm)+'\n)\n')
 (pm/'owner').write_text(head('owner','labelList')+str(len(perm))+'\n(\n'+'\n'.join(map(str,g['owner'][perm]))+'\n)\n')
 patches=g['patches'];bb=head('boundary','polyBoundaryMesh')+str(len(patches)+1)+'\n(\n';start=n
 for b in patches:
  count=sum(i not in ids for i in range(b['start'],b['start']+b['count']))
  bb+=f'{b["name"]}\n{{ type {"wall" if b["wall"] else "patch"}; nFaces {count}; startFace {start}; }}\n';start+=count
 bb+=f'ceiling\n{{ type wall; nFaces {len(ids)}; startFace {start}; }}\n)\n';(pm/'boundary').write_text(bb)
 # Remove stale face sets/zones whose labels refer to original ordering.
 for f in ['faceZones','cellZones','pointZones']:
  if (pm/f).exists():(pm/f).unlink()
 if (pm/'sets').exists():shutil.rmtree(pm/'sets')
 for f in (out/'0').iterdir():
  t=f.read_text();m=re.search(r'\bwalls\s*\{([^{}]*)\}',t);assert m,f
  bc=m[1] if f.name!='T' else 'type fixedValue; value uniform 306.15;'
  t=t[:t.rfind('}')]+f'ceiling {{ {bc} }}\n'+t[t.rfind('}'):]
  if f.name=='U':t=re.sub(r'\bhood\s*\{[^{}]*\}', 'hood { type flowRateOutletVelocity; volumetricFlowRate constant 0.5; extrapolateProfile no; value uniform (0 0 1.85185185185); }',t)
  if f.name=='p_rgh':t=re.sub(r'\bhood\s*\{[^{}]*\}', 'hood { type fixedFluxPressure; value uniform 101325; }',t)
  if fine and f.name=='p':
   z=np.load(out/'mesh_reference.npz')['centres'][:,2];v=101325*np.exp(-9.81*z/(287.05*298.15))
   t=re.sub(r'internalField.*?;', 'internalField nonuniform List<scalar> '+str(len(z))+' (\n'+'\n'.join(map(str,v))+'\n);',t,count=1,flags=re.S)
  f.write_text(t)
 p=out/'constant/aerosolProperties';t=p.read_text();t=re.sub(r'hoodTargetQ\s+[^;]+;','hoodTargetQ 0.5;',t);t+='\nprescribedHoodFlow true;\n';p.write_text(t)
 p=out/'system/controlDict';t=p.read_text();t=re.sub(r'application\s+[^;]+;','application kitchenCeilingAirFoam;',t);t=t.replace('startFrom latestTime','startFrom startTime');t=re.sub(r'endTime\s+[^;]+;',f'endTime {end};',t);t=re.sub(r'writeInterval\s+[^;]+;','writeInterval 2;',t);
 if doors!='both_closed':
  t=re.sub(r'maxCo\s+[^;]+;','maxCo 6;',t);t=re.sub(r'maxDeltaT\s+[^;]+;','maxDeltaT 0.1;',t)
 p.write_text(t)
 gg=mesh_geometry(out)
 assert np.array_equal(g['points'],gg['points']);assert np.array_equal(g['owner'][:n],gg['owner'][:n]);assert abs(gg['area'][gg['ceiling_ids']].sum()-4.245)<1e-8
 (out/'boundary_update.json').write_text(json.dumps(dict(source=str(source),ceiling_temperature_C=33,prescribed_flow_m3_min=30,flow_mode='prescribed outlet velocity; required pressure diagnosed, NOT a promise that 1300 Pa fan achieves this flow',old_ceiling_face_ids=list(map(int,sorted(ids))),new_ceiling_face_ids=gg['ceiling_ids'].tolist(),points_unchanged=True,cell_topology_unchanged=True,internal_faces_unchanged=True,end_s=end),indent=2))
 return out
if __name__=='__main__':
 for layout in ['corner','center']:
  for doors in ['both_open','both_closed','left_open','right_open']:print(prepare(layout,doors))
 # Fine only used for a matched low-cost test, never presented as full event.
 print(prepare('corner','both_closed',True,20))
