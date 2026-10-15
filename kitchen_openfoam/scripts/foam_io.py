"""Read native OpenFOAM ASCII cell fields without interpolation."""
from pathlib import Path
import re,numpy as np,csv

# OpenCFD 2406 defaults: NA=6.0221417930e23, k=1.38065e-23.
# Reconstruct density with the solver's constant, rather than modern CODATA,
# so exported cell mass matches native diagnostics to ASCII write precision.
FOAM_RR=6.0221417930e23*1.38065e-23*1000
AIR_MOLAR_MASS=28.965

def read_internal(path,n=None):
    text=Path(path).read_text()
    m=re.search(r'internalField\s+nonuniform\s+List<(\w+)>\s+(\d+)\s*\((.*?)\)\s*;',text,re.S)
    if m:
        values=np.fromstring(m[3].replace('(',' ').replace(')',' '),sep=' ')
        count=int(m[2]);dim=3 if m[1]=='vector' else 1
        assert values.size==count*dim,(path,values.size,count,dim)
        return values.reshape(count,dim) if dim>1 else values
    m=re.search(r'internalField\s+uniform\s+([^;]+);',text)
    if not m:raise ValueError(f'No internal field: {path}')
    a=np.fromstring(m[1].replace('(',' ').replace(')',' '),sep=' ')
    if n is None:raise ValueError('Need cell count for a uniform field')
    return np.repeat(a[None,:],n,axis=0) if len(a)>1 else np.full(n,a[0])

def times(case):
    return sorted(float(p.name) for p in Path(case).iterdir() if p.is_dir() and re.fullmatch(r'\d+(?:\.\d+)?',p.name))

def time_dir(case,t):
    choices=[p for p in Path(case).iterdir() if p.is_dir() and re.fullmatch(r'\d+(?:\.\d+)?',p.name)]
    q=min(choices,key=lambda p:abs(float(p.name)-t))
    if abs(float(q.name)-t)>1e-6:raise ValueError(f'Missing saved time {t}: {case}')
    return q

def load_mesh(case):
    a=dict(np.load(Path(case)/'mesh_reference.npz'))
    # Exact centroid of each extruded polygon; preserve OpenFOAM cell ordering.
    centre=np.empty_like(a['centres'])
    for nv in [6,8]:
        mask=a['nverts']==nv;idx=np.flatnonzero(mask);q=a['points'][a['cells'][mask,:nv//2]]
        yz=q[:,:,1:];yn=np.roll(yz,-1,axis=1);cross=yz[:,:,0]*yn[:,:,1]-yn[:,:,0]*yz[:,:,1]
        area2=cross.sum(1);c=((yz+yn)*cross[:,:,None]).sum(1)/(3*area2[:,None])
        centre[idx,0]=a['points'][a['cells'][mask,:nv],0].mean(1);centre[idx,1:]=c
    a['centres']=centre
    native=Path(case)/'0/C'
    if native.exists():
        c=read_internal(native,len(centre));assert np.max(abs(c-centre))<1e-5
        a['centres']=c
    return a

def history(case):
    # A restart writes a separate CSV. Sort by its starting physical time.
    files=sorted((Path(case)/'diagnostics').glob('aerosol_*.csv'),key=lambda p:float(p.stem.split('_',1)[1]))
    records={}
    for p in files:
        start=float(p.stem.split('_',1)[1])
        # A continuation replaces the old branch after its saved restart time.
        if start>0:records={k:v for k,v in records.items() if k[0]<=start+1e-9}
        for row in csv.DictReader(p.open()):
            try:d={k:float(v) for k,v in row.items() if k is not None}
            except (TypeError,ValueError):continue
            if len(d)<25:continue
            records[(round(d['time_s'],10),d['diameter_um'])]=d
    counts={}
    for tt,d in records:counts[tt]=counts.get(tt,0)+1
    end_match=re.search(r'endTime\s+([\d.]+)',(Path(case)/'system/controlDict').read_text())
    end=float(end_match[1]) if end_match else np.inf
    rows=[records[k] for k in sorted(records) if counts[k[0]]==5 and k[0]<=end+1e-7]
    return rows

def fields(case,t,mesh=None):
    mesh=load_mesh(case) if mesh is None else mesh;n=len(mesh['cells']);d=time_dir(case,t)
    f={k:read_internal(d/k,n) for k in ['U','T','p','p_rgh','k','omega','nut','Yp01','Yp03','Yp1','Yp3','Yp10']}
    f['rho']=f['p']/((FOAM_RR/AIR_MOLAR_MASS)*f['T']);f['speed']=np.linalg.norm(f['U'],axis=1)
    f['C_over_R']=np.array([f['rho']*f[k]/1e-9 for k in ['Yp01','Yp03','Yp1','Yp3','Yp10']]).T
    return f
