"""Conforming Gmsh tetrahedra with explicitly modeled 45-degree hood baffle."""
from pathlib import Path
import argparse,json,math,numpy as np
import gmsh,meshio
from skfem import MeshTet

def build(folder,h=.18):
 folder=Path(folder);c=json.loads((folder/'config.json').read_text());r=c['room'];k=c['counter'];f=c['hood'];b=c['baffle'];d=c['doors']
 L,D,H=r['length'],r['depth'],r['height'];x0=k['x0'];dep=k['depth'];z0=k['height'];ud=c['upper_cabinet_depth_m'];hz=f['z']
 gmsh.initialize();gmsh.option.setNumber('General.Terminal',0);gmsh.model.add('kitchen_'+c['layout']);o=gmsh.model.occ
 air=o.addBox(0,0,0,L,D,H)
 solids=[o.addBox(0,0,0,x0,dep,H),o.addBox(x0,0,0,L-x0,dep,z0),o.addBox(x0,0,hz,L-x0,ud,H-hz)]
 a=math.radians(b['down_angle_deg']);dy,dz=math.cos(a),-math.sin(a);ny,nz=math.sin(a),math.cos(a)
 half=b['thickness_m_assumed']/2;rooty=f['y1']-.003*dy;rootz=hz-.003*dz;tipy=f['y1']+b['length_m']*dy;tipz=hz+b['length_m']*dz
 yz=[(rooty+half*ny,rootz+half*nz),(rooty-half*ny,rootz-half*nz),(tipy-half*ny,tipz-half*nz),(tipy+half*ny,tipz+half*nz)]
 ps=[o.addPoint(f['x0'],y,z) for y,z in yz];ls=[o.addLine(ps[i],ps[(i+1)%4]) for i in range(4)]
 sf=o.addPlaneSurface([o.addCurveLoop(ls)]);ex=o.extrude([(2,sf)],f['x1']-f['x0'],0,0)
 solids.extend(tag for dim,tag in ex if dim==3)
 cut,_=o.cut([(3,air)],[(3,s) for s in solids],removeObject=True,removeTool=True)
 def rectangle(points):
  pp=[o.addPoint(*p) for p in points];ll=[o.addLine(pp[i],pp[(i+1)%4]) for i in range(4)]
  return o.addPlaneSurface([o.addCurveLoop(ll)])
 tools=[]
 for x in [0,L]:
  for za,zb in [(0,.05),(.05,d['homogenized_inlet_height_m']),(d['homogenized_inlet_height_m'],H)]:
   tools.append((2,rectangle([(x,d['y0'],za),(x,d['y1'],za),(x,d['y1'],zb),(x,d['y0'],zb)])))
 tools.append((2,rectangle([(f['x0'],f['y0'],hz),(f['x1'],f['y0'],hz),(f['x1'],f['y1'],hz),(f['x0'],f['y1'],hz)])))
 domain,_=o.fragment(cut,tools,removeObject=True,removeTool=True);o.synchronize()
 volumes=[tag for dim,tag in domain if dim==3]
 pg=gmsh.model.addPhysicalGroup(3,volumes);gmsh.model.setPhysicalName(3,pg,'air')
 gmsh.option.setNumber('Mesh.MeshSizeMax',h);gmsh.option.setNumber('Mesh.MeshSizeMin',.008)
 gmsh.option.setNumber('Mesh.MeshSizeFromCurvature',0);gmsh.option.setNumber('Mesh.MeshSizeExtendFromBoundary',0)
 gmsh.option.setNumber('Mesh.MeshSizeFromPoints',0);gmsh.option.setNumber('Mesh.Algorithm3D',1)
 # Airflow source neighborhoods require at least ~8 cm elements.
 fs=[]
 for bx,by,bz in c['burners']:
  fld=gmsh.model.mesh.field.add('Ball');gmsh.model.mesh.field.setNumber(fld,'XCenter',bx);gmsh.model.mesh.field.setNumber(fld,'YCenter',by);gmsh.model.mesh.field.setNumber(fld,'ZCenter',bz+.13)
  for name,val in [('Radius',.24),('Thickness',.15),('VIn',min(.09,.5*h)),('VOut',h)]:gmsh.model.mesh.field.setNumber(fld,name,val)
  fs.append(fld)
 fld=gmsh.model.mesh.field.add('Box')
 vals={'XMin':f['x0']-.06,'XMax':f['x1']+.06,'YMin':f['y1']-.08,'YMax':tipy+.08,'ZMin':tipz-.08,'ZMax':hz+.04,'VIn':min(.07,.4*h),'VOut':h,'Thickness':.12}
 for name,val in vals.items():gmsh.model.mesh.field.setNumber(fld,name,val)
 fs.append(fld)
 mf=gmsh.model.mesh.field.add('Min');gmsh.model.mesh.field.setNumbers(mf,'FieldsList',fs);gmsh.model.mesh.field.setAsBackgroundMesh(mf)
 gmsh.model.mesh.generate(3);gmsh.model.mesh.optimize('Netgen')
 tags,xyz,_=gmsh.model.mesh.getNodes();xyz=np.array(xyz).reshape(-1,3);tag_to_index={int(tag):i for i,tag in enumerate(tags)}
 types,_,nodes=gmsh.model.mesh.getElements(3);tt=[]
 for typ,nd in zip(types,nodes):
  if typ==4:tt.append(np.array([tag_to_index[int(i)] for i in nd]).reshape(-1,4))
 tetra=np.vstack(tt);used=np.unique(tetra);remap=np.full(len(xyz),-1);remap[used]=np.arange(len(used));points=xyz[used];tetra=remap[tetra]
 gmsh.write(str(folder/'results'/f'model_h{h:.2f}.msh'));gmsh.finalize()
 mesh=MeshTet(points.T,tetra.T.astype(np.int32));volume=float(np.sum(np.abs(np.linalg.det(points[tetra[:,1:]]-points[tetra[:,0],None,:]))/6))
 np.savez_compressed(folder/'results'/f'mesh_h{h:.2f}.npz',points=mesh.p.T,cells=mesh.t.T)
 meshio.write(folder/'results'/f'geometry_mesh_h{h:.2f}.vtu',meshio.Mesh(mesh.p.T,[('tetra',mesh.t.T)]))
 meta={'layout':c['layout'],'mesh_target_m':h,'nodes':int(mesh.nvertices),'tetrahedra':int(mesh.nelements),'fluid_volume_m3':volume,'baffle_thickness_m_assumed':b['thickness_m_assumed'],'baffle_width_m':f['x1']-f['x0'],'baffle_tip_yz_m':[tipy,tipz],'method':'Gmsh OCC conforming tetrahedra with solid subtraction and partitioned hood/door boundaries'}
 (folder/'results'/f'mesh_h{h:.2f}_metadata.json').write_text(json.dumps(meta,ensure_ascii=False,indent=2));print(json.dumps(meta,ensure_ascii=False),flush=True)
 return mesh
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('--layout',choices=['corner','center'],required=True);ap.add_argument('--h',type=float,default=.18);a=ap.parse_args();build(Path(__file__).resolve().parent/a.layout,a.h)
