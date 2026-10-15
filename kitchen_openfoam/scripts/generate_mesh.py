"""Gmsh conforming kitchen air mesh with zero-thickness embedded hood baffle.
Same confirmed dimensions as existing FEM model; no 2 mm assumption required.
"""
from pathlib import Path
import sys,json,math,numpy as np,gmsh
ROOT=Path(__file__).resolve().parents[1]

def generate(layout,h=.14):
    out=ROOT/'meshes';out.mkdir(exist_ok=True)
    path=out/f'{layout}_h{h:.2f}.npz'
    if path.exists():return path
    c=json.loads((ROOT.parent/'kitchen_compare'/layout/'config.json').read_text())
    gmsh.initialize();gmsh.option.setNumber('General.Terminal',0);gmsh.model.add(layout);o=gmsh.model.occ
    air=o.addBox(0,0,0,3.3,1.65,2.2)
    solids=[o.addBox(0,0,0,.7,.6,2.2),o.addBox(.7,0,0,2.6,.6,.9),o.addBox(.7,0,1.6,2.6,.3,.6)]
    air,_=o.cut([(3,air)],[(3,s) for s in solids])
    def rect(ps):
        pts=[o.addPoint(*p) for p in ps];ls=[o.addLine(pts[i],pts[(i+1)%4]) for i in range(4)]
        return (2,o.addPlaneSurface([o.addCurveLoop(ls)]))
    hd=c['hood'];xa,xb=hd['x0'],hd['x1'];tip=.1/math.sqrt(2)
    tools=[rect([(xa,.3,1.6),(xb,.3,1.6),(xb,.3+tip,1.6-tip),(xa,.3+tip,1.6-tip)]),rect([(xa,0,1.6),(xb,0,1.6),(xb,.3,1.6),(xa,.3,1.6)])]
    for x in [0,3.3]:
        for za,zb in [(0,.1),(.1,2.2)]:tools.append(rect([(x,.7,za),(x,1.6,za),(x,1.6,zb),(x,.7,zb)]))
    dom,_=o.fragment(air,tools);o.synchronize()
    vs=[tag for dim,tag in dom if dim==3];pg=gmsh.model.addPhysicalGroup(3,vs);gmsh.model.setPhysicalName(3,pg,'air')
    for k,v in {'Mesh.MeshSizeMax':h,'Mesh.MeshSizeMin':.03,'Mesh.MeshSizeFromCurvature':0,'Mesh.MeshSizeExtendFromBoundary':0,'Mesh.MeshSizeFromPoints':0,'Mesh.Algorithm3D':1}.items():gmsh.option.setNumber(k,v)
    fs=[]
    for x,y,z in c['burners']:
        f=gmsh.model.mesh.field.add('Ball')
        for k,v in {'XCenter':x,'YCenter':y,'ZCenter':z+.13,'Radius':.25,'Thickness':.15,'VIn':.5*h,'VOut':h}.items():gmsh.model.mesh.field.setNumber(f,k,v)
        fs.append(f)
    f=gmsh.model.mesh.field.add('Box')
    for k,v in {'XMin':xa-.08,'XMax':xb+.08,'YMin':.22,'YMax':.46,'ZMin':1.44,'ZMax':1.66,'VIn':.5*h,'VOut':h,'Thickness':.12}.items():gmsh.model.mesh.field.setNumber(f,k,v)
    fs.append(f);mf=gmsh.model.mesh.field.add('Min');gmsh.model.mesh.field.setNumbers(mf,'FieldsList',fs);gmsh.model.mesh.field.setAsBackgroundMesh(mf)
    gmsh.model.mesh.generate(3);gmsh.model.mesh.optimize('')
    tags,xyz,_=gmsh.model.mesh.getNodes();xyz=np.asarray(xyz).reshape(-1,3);mp={int(t):i for i,t in enumerate(tags)}
    types,_,nds=gmsh.model.mesh.getElements(3);tt=[]
    for typ,ns in zip(types,nds):
        if typ==4:tt.append(np.array([mp[int(i)] for i in ns]).reshape(-1,4))
    t=np.vstack(tt);used=np.unique(t);remap=np.full(len(xyz),-1);remap[used]=np.arange(len(used));p=xyz[used];t=remap[t]
    gmsh.write(str(out/f'{layout}_h{h:.2f}.msh'));gmsh.finalize()
    np.savez_compressed(path,points=p,cells=t)
    print('Gmsh',layout,h,len(p),len(t),flush=True);return path

if __name__=='__main__':
    for l in ['corner','center']:generate(l,float(sys.argv[1]) if len(sys.argv)>1 else .14)
