"""Minimal reader for the existing ASCII OpenFOAM polyMesh; no remeshing."""
from pathlib import Path
import re,numpy as np

def list_body(p):
 s=Path(p).read_text();m=re.search(r'\n\s*(\d+)\s*\n\s*\(',s)
 if not m:raise ValueError(p)
 return int(m[1]),s[m.end():s.rfind(')')]
def mesh_geometry(case):
 pm=Path(case)/'constant/polyMesh'
 n,s=list_body(pm/'points');points=np.fromstring(s.replace('(',' ').replace(')',' '),sep=' ').reshape(n,3)
 n,s=list_body(pm/'faces');faces=[np.fromstring(m[2],sep=' ',dtype=int) for m in re.finditer(r'(\d+)\(([^)]*)\)',s)]
 assert len(faces)==n
 _,s=list_body(pm/'owner');owner=np.fromstring(s,sep=' ',dtype=int)
 patches=[]
 for m in re.finditer(r'(\w+)\s*\{([^{}]*)\}',(pm/'boundary').read_text()):
  if 'nFaces' not in m[2]:continue
  patches.append(dict(name=m[1],start=int(re.search(r'startFace\s+(\d+)',m[2])[1]),count=int(re.search(r'nFaces\s+(\d+)',m[2])[1]),wall=bool(re.search(r'type\s+wall\s*;',m[2]))))
 cc=np.array([points[f].mean(0) for f in faces]);sv=np.array([np.cross(points[f],np.roll(points[f],-1,axis=0)).sum(0)/2 for f in faces]);area=np.linalg.norm(sv,axis=1)
 ids=np.flatnonzero((abs(cc[:,2]-2.2)<1e-7)&(sv[:,2]/area>.99))
 return dict(points=points,faces=faces,owner=owner,patches=patches,centres=cc,area=area,ceiling_ids=ids)
def write_ceiling_set(case):
 a=mesh_geometry(case);ids=a['ceiling_ids'];pm=Path(case)/'constant/polyMesh';(pm/'sets').mkdir(exist_ok=True)
 (pm/'sets/ceilingFaces').write_text('FoamFile { version 2.0; format ascii; class faceSet; object ceilingFaces; }\n'+f'{len(ids)}\n(\n'+'\n'.join(map(str,ids))+'\n)\n')
 (pm/'faceZones').write_text('FoamFile { version 2.0; format ascii; class regIOobject; object faceZones; }\n1\n(\nceiling\n{\ntype faceZone;\nfaceLabels List<label> '+str(len(ids))+'\n(\n'+'\n'.join(map(str,ids))+'\n);\nflipMap List<bool> '+str(len(ids))+'\n(\n'+'\n'.join(['false']*len(ids))+'\n);\n}\n)\n')
 return len(ids),float(a['area'][ids].sum())
