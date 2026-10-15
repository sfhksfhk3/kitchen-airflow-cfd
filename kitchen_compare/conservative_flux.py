"""H(div) RT0 flux reconstruction of the finite-element velocity field.
Keep hood/wall normal fluxes, project makeup/internal fluxes to cell-wise balance.
This removes artificial tracer compression from stabilized P1 velocity divergence.
"""
from pathlib import Path
import numpy as np
from scipy.sparse import coo_matrix,diags
from scipy.sparse.linalg import spsolve

def reconstruct(mesh,U,open_faces,hood_faces):
 p=mesh.p.T;t=mesh.t.T;fac=mesh.facets.T;nf=len(fac);ne=len(t)
 v=p[fac];fc=v.mean(axis=1);av=np.cross(v[:,1]-v[:,0],v[:,2]-v[:,0])/2;area=np.linalg.norm(av,axis=1);un=av/area[:,None]
 tc=p[t].mean(axis=1);tf=mesh.t2f.T
 signs=np.sign(np.sum(av[tf]*(fc[tf]-tc[:,None,:]),axis=2))
 B=coo_matrix((signs.ravel(),(np.repeat(np.arange(ne),4),tf.ravel())),shape=(ne,nf)).tocsr()
 raw=np.sum(U.T[fac].mean(axis=1)*av,axis=1);flux=raw.copy()
 isboundary=mesh.f2t[1]<0;allowed=np.zeros(nf,bool);allowed[open_faces]=True;allowed[hood_faces]=True
 flux[isboundary&~allowed]=0
 free=np.flatnonzero(~isboundary|np.isin(np.arange(nf),open_faces))
 own0=mesh.f2t[0];own1=mesh.f2t[1]
 distance=np.abs(np.sum((fc-tc[own0])*un,axis=1))
 inside=own1>=0;distance[inside]+=np.abs(np.sum((fc[inside]-tc[own1[inside]])*un[inside],axis=1))
 w=area[free]/np.maximum(distance[free],1e-8);BB=B[:,free]
 potential=spsolve((BB@diags(w)@BB.T).tocsc(),B@flux)
 flux[free]-=w*(BB.T@potential)
 local=flux[tf]*signs;imbalance=local.sum(axis=1)
 vol=np.abs(np.linalg.det(p[t[:,1:]]-p[t[:,0],None,:]))/6
 # Opposite vertex for each reference face, determined without ordering assumptions.
 opposite=np.array([[int(next(i for i in tet if i not in fac[f])) for f in faces] for tet,faces in zip(t,tf)])
 # RT0: sum q_i (x-v_i)/(3V). Cell-balanced flux makes the x coefficient zero.
 urt=-np.sum(local[:,:,None]*p[opposite],axis=1)/(3*vol[:,None])
 rawmean=U.T[t].mean(axis=1);correction=float(np.sqrt(np.sum(np.sum((urt-rawmean)**2,axis=1)*vol)/np.sum(np.sum(rawmean**2,axis=1)*vol)))
 meta={'method':'RT0 H(div) least-energy normal-flux reconstruction from P1 FEM velocities','cell_flux_max_imbalance_m3_s':float(np.max(abs(imbalance))),'relative_RMS_cell_velocity_correction':correction,'boundary_fixed':'hood normal flux fixed to FEM value, wall flux zero; door and interior fluxes projected','max_raw_cell_divergence_s_inv':float(np.max(abs(B@raw)/vol)),'max_projected_cell_divergence_s_inv':float(np.max(abs(imbalance)/vol))}
 assert np.max(abs(imbalance))<1e-8
 return urt,flux,area,meta
