"""Finite-volume-oriented hex/prism mesh: recombined Gmsh y-z slice extruded x.

Retains the exact 45 degree embedded baffle and all confirmed dimensions.
Gmsh quadrilaterals are extruded to hexahedra; residual triangles to prisms.
"""
from pathlib import Path
import json,numpy as np,gmsh
ROOT=Path(__file__).resolve().parents[1]

def generate(layout,h=.09,gap_layers=2):
    graded='_graded' if gap_layers>2 else ''
    cache=ROOT/'meshes'/f'{layout}_extruded_gap005_L{gap_layers}{graded}_h{h:.3f}.npz'
    cache.parent.mkdir(parents=True,exist_ok=True)
    if cache.exists():return cache
    cfg=json.loads((ROOT/'inputs'/layout/'config.json').read_text())
    gmsh.initialize();gmsh.option.setNumber('General.Terminal',0);gmsh.model.add('yz');o=gmsh.model.occ
    # 2D coordinates are (y,z); conversion to (x,y,z) is below.
    s=o.addRectangle(0,0,0,1.65,2.2)
    base=o.addRectangle(0,0,0,.6,.9);upper=o.addRectangle(0,1.6,0,.3,.6)
    sf,_=o.cut([(2,s)],[(2,base),(2,upper)])
    def line(a,b):return (1,o.addLine(o.addPoint(*a,0),o.addPoint(*b,0)))
    tip=.1/np.sqrt(2)
    tools=[line((.3,1.6),(.3+tip,1.6-tip)),line((.6,0),(.6,2.2)),line((.7,0),(.7,2.2)),line((1.6,0),(1.6,2.2))]+[line((.6,.005*j/gap_layers),(1.65,.005*j/gap_layers)) for j in range(1,gap_layers+1)]
    if gap_layers>2:tools += [line((.6,z),(1.65,z)) for z in [.0075,.01,.015,.02,.03,.04,.06,.08]]
    dom,_=o.fragment(sf,tools);o.synchronize()
    surfaces=[tag for dim,tag in dom if dim==2]
    pg=gmsh.model.addPhysicalGroup(2,surfaces);gmsh.model.setPhysicalName(2,pg,'yzAir')
    for key,val in {'Mesh.MeshSizeMax':h,'Mesh.MeshSizeMin':h*.45,'Mesh.MeshSizeFromPoints':0,'Mesh.MeshSizeExtendFromBoundary':0,'Mesh.MeshSizeFromCurvature':0,'Mesh.Algorithm':8,'Mesh.RecombineAll':1}.items():gmsh.option.setNumber(key,val)
    fields=[]
    for y,z,radius in [(.3,.98,.22),(.3,1.58,.18)]:
        f=gmsh.model.mesh.field.add('Ball')
        for k,v in {'XCenter':y,'YCenter':z,'ZCenter':0,'Radius':radius,'Thickness':.1,'VIn':h*.55,'VOut':h}.items():gmsh.model.mesh.field.setNumber(f,k,v)
        fields.append(f)
    mf=gmsh.model.mesh.field.add('Min');gmsh.model.mesh.field.setNumbers(mf,'FieldsList',fields);gmsh.model.mesh.field.setAsBackgroundMesh(mf)
    gmsh.model.mesh.generate(2)
    tags,xyz,_=gmsh.model.mesh.getNodes();yz=np.asarray(xyz).reshape(-1,3)[:,:2];mp={int(tag):i for i,tag in enumerate(tags)}
    types,_,nodes=gmsh.model.mesh.getElements(2);polys=[]
    for ty,nd in zip(types,nodes):
        if ty not in [2,3]:continue
        n=3 if ty==2 else 4
        polys.extend(np.array([mp[int(k)] for k in nd]).reshape(-1,n).tolist())
    gmsh.write(str(cache.with_suffix('.msh')));gmsh.finalize()
    cuts=sorted(set([0,.7,3.3,cfg['hood']['x0'],cfg['hood']['x1']]))
    xs=[]
    for xa,xb in zip(cuts[:-1],cuts[1:]):xs.extend(np.linspace(xa,xb,int(np.ceil((xb-xa)/h))+1)[:-1])
    xs=np.r_[xs,3.3];ny=len(yz)
    points=np.array([(x,y,z) for x in xs for y,z in yz])
    cellpoints=[];nverts=[];centres=[];volumes=[]
    for j,(xa,xb) in enumerate(zip(xs[:-1],xs[1:])):
        for poly in polys:
            q=yz[poly];cy,cz=q.mean(0)
            if (xa+xb)/2<.7 and cy<.6:continue
            ar=abs(np.sum(q[:,0]*np.roll(q[:,1],-1)-q[:,1]*np.roll(q[:,0],-1)))/2
            cc=[j*ny+i for i in poly]+[(j+1)*ny+i for i in poly]
            cellpoints.append(cc+[-1]*(8-len(cc)));nverts.append(len(cc));centres.append([(xa+xb)/2,cy,cz]);volumes.append(ar*(xb-xa))
    cp=np.array(cellpoints);used=np.unique(cp[cp>=0]);remap=np.full(len(points),-1);remap[used]=np.arange(len(used));valid=cp>=0;cp[valid]=remap[cp[valid]];points=points[used]
    np.savez_compressed(cache,points=points,cells=cp,nverts=nverts,centres=centres,volumes=volumes)
    print('EXTRUDED',layout,h,len(points),len(cp),np.sum(np.array(nverts)==8),'hexahedra; volume',sum(volumes),flush=True)
    return cache

if __name__=='__main__':
    import argparse
    ap=argparse.ArgumentParser();ap.add_argument('--layout',default='corner');ap.add_argument('--h',type=float,default=.09);a=ap.parse_args();generate(a.layout,a.h)
