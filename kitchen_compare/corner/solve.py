"""Reproducible 3D tetrahedral P1/P1 SUPG/PSPG finite-element kitchen model.
Preliminary steady fixed-reference-density incompressible model with nonlinear ideal-gas buoyancy, constant eddy viscosity.
Pressure-natural makeup openings; prescribed hood volume flow; no-slip solids.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
import argparse, json, time, pathlib
import numpy as np
from scipy.sparse import bmat, block_diag
from scipy.sparse.linalg import spsolve
from skfem import MeshTet, Basis, FacetBasis, ElementTetP1, BilinearForm, LinearForm, asm
from skfem.helpers import dot, grad
import meshio
ROOT=pathlib.Path(__file__).resolve().parent
CFG=json.loads((ROOT/'config.json').read_text())
P=CFG['physics']
L,DEP,H=[CFG['room'][k] for k in ['length','depth','height']]
CX,CD,CZ=[CFG['counter'][k] for k in ['x0','depth','height']]
HX0,HX1,HY0,HY1,HZ=[CFG['hood'][k] for k in ['x0','x1','y0','y1','z']]
DY0,DY1=CFG['doors']['y0'],CFG['doors']['y1']
SLOT=CFG['doors']['homogenized_inlet_height_m']
HEAT=P['sensible_heat_each_W']*len(CFG['burners'])
T0=P['ambient_C']

@BilinearForm
def mass(u,v,w): return u*v
@BilinearForm
def diffusion(u,v,w): return dot(grad(u),grad(v))
@BilinearForm
def transport(u,v,w):
    adv=dot(w.vel,grad(u)); test=dot(w.vel,grad(v))
    return w.nu*dot(grad(u),grad(v))+(adv+w.invdt*u)*(v+w.tau*test)
@BilinearForm
def pressure_term(u,v,w):
    return -u*grad(v)[w.i]+w.tau*grad(u)[w.i]*dot(w.vel,grad(v))
@BilinearForm
def continuity_term(u,v,w):
    return v*grad(u)[w.i]+w.tau*grad(v)[w.i]*(dot(w.vel,grad(u))+w.invdt*u)
@BilinearForm
def p_stab(u,v,w): return w.tau*dot(grad(u),grad(v))
@LinearForm
def body_force(v,w): return w.force*(v+w.tau*dot(w.vel,grad(v)))
@LinearForm
def p_force(v,w): return w.tau*dot(grad(v),w.fvec)
@BilinearForm
def backflow(u,v,w): return w.penalty*u*v
@BilinearForm
def grad_div(u,v,w): return .5*grad(u)[w.j]*grad(v)[w.i]
@LinearForm
def heat_load(v,w): return w.source*(v+w.tau*dot(w.vel,grad(v)))
@BilinearForm
def heat_galerkin(u,v,w): return w.nu*dot(grad(u),grad(v))-u*dot(w.vel,grad(v))
@BilinearForm
def heat_boundary(u,v,w): return w.vn*u*v

def thermal_matrix(bs,fb,U,alpha,wallmat):
    vel=np.array([bs.interpolate(ui) for ui in U]);uf=np.array([fb.interpolate(ui) for ui in U])
    vn=np.sum(uf*fb.normals,axis=0)
    return monotone_heat(asm(heat_galerkin,bs,vel=vel,nu=alpha)+asm(heat_boundary,fb,vn=vn)+wallmat)
@LinearForm
def positive_heat_load(v,w):return w.source*v

def monotone_heat(A):
    off=A.tocsr().copy();off.setdiag(0);off.eliminate_zeros()
    pos=off.maximum(0);edge=pos.maximum(pos.T)
    from scipy.sparse import diags
    return (A-edge+diags(np.asarray(edge.sum(axis=1)).ravel())).tocsr()

@LinearForm
def integral_basis(v,w): return v


def axis(points,h):
    vals=[]
    for a,b in zip(points[:-1],points[1:]):
        vals.extend(np.linspace(a,b,max(1,int(np.ceil((b-a)/h)))+1)[:-1])
    return np.array(vals+[points[-1]])

def make_mesh(h):
    path=ROOT/'results'/f'mesh_h{h:.2f}.npz'
    if not path.exists():
        import sys
        sys.path.insert(0,str(ROOT.parent))
        from build_mesh import build
        return build(ROOT,h)
    d=np.load(path)
    return MeshTet(d['points'].T,d['cells'].T)

def faces(m,case):
    f=m.boundary_facets(); c=m.p[:,m.facets[:,f]].mean(axis=1)
    near=lambda a,b:np.abs(a-b)<1e-7
    hood=near(c[2],HZ)&(c[0]>HX0)&(c[0]<HX1)&(c[1]>HY0)&(c[1]<HY1)
    dl=near(c[0],0)&(c[1]>DY0)&(c[1]<DY1)
    dr=near(c[0],L)&(c[1]>DY0)&(c[1]<DY1)
    if case=='both_closed': dl&=c[2]<SLOT;dr&=c[2]<SLOT
    elif case=='left_open':dr[:]=False
    elif case=='right_open':dl[:]=False
    wall=~(hood|dl|dr)
    return {k:f[v] for k,v in [('hood',hood),('left',dl),('right',dr),('walls',wall)]}

def fan_estimate(case):
    # Shared assumed fan curve / orifice model. Neglect the single closed gap
    # in one-open cases (<0.23% of open-door geometric area).
    d=CFG['doors'];da=d['width']*d['height'];ga=d['width']*d['closed_gap_m']
    A={'both_open':2*da,'both_closed':2*ga,'left_open':da,'right_open':da}[case]
    rho=P['rho'];cd=CFG['fan']['discharge_coefficient_assumed'];pmax=CFG['fan']['max_static_Pa'];qmax=CFG['fan']['free_delivery_m3_s']
    R=rho/(2*(cd*A)**2)
    q=np.sqrt(pmax/(R+pmax/qmax**2))
    return q,R*q*q,A

def facet_flux(m,fs,U):
    if len(fs)==0:return 0.,0.,0.
    fb=FacetBasis(m,ElementTetP1(),facets=fs,intorder=3)
    uf=np.array([fb.interpolate(U[i]) for i in range(3)])
    vn=np.sum(uf*fb.normals,axis=0)
    return float(np.sum(vn*fb.dx)),float(np.sum(np.maximum(-vn,0)*fb.dx)),float(np.sum(np.maximum(vn,0)*fb.dx))

def run(case,h=.18,nu_t=None,suffix='',maxiter=300,warm=None):
    start=time.time(); nu_t=P['nu_t_assumed'] if nu_t is None else nu_t
    nu=P['nu_molecular']+nu_t;alpha=P['alpha_molecular']+nu_t/P['Pr_t']
    m=make_mesh(h);bs=Basis(m,ElementTetP1(),intorder=3);n=bs.N;f=faces(m,case)
    M=asm(mass,bs);K=asm(diffusion,bs)
    q,dp,area=fan_estimate(case)
    wn=np.unique(m.facets[:,f['walls']]);hn=np.unique(m.facets[:,f['hood']])
    dn=np.unique(np.r_[m.facets[:,f['left']].ravel(),m.facets[:,f['right']].ravel()])
    velD=np.unique(np.r_[wn,hn]);D=np.concatenate([velD+i*n for i in range(3)])
    free=np.setdiff1d(np.arange(4*n),D)
    u_bc=np.zeros((3,n));x,y=m.p[:2,hn]
    profile=np.maximum(0,np.sin(np.pi*(x-HX0)/(HX1-HX0)))*np.maximum(0,np.sin(np.pi*(y-HY0)/(HY1-HY0)))
    u_bc[2,hn]=profile
    u_bc[:,wn]=0
    raw=facet_flux(m,f['hood'],u_bc)[0];u_bc*=q/raw
    sol=np.r_[u_bc.ravel(),np.zeros(n)];U=sol[:3*n].reshape(3,n).copy();theta=np.zeros(n)
    checkpoint_history=[]; first_iteration=0; elapsed_before=0.
    if warm:
        from scipy.interpolate import LinearNDInterpolator
        w=np.load(warm)
        same_mesh=w['points'].shape==m.p.T.shape and np.allclose(w['points'],m.p.T,rtol=0,atol=1e-12)
        values=np.c_[w['velocity'],w['pressure']/P['rho'],w['temperature']-T0]
        vals=values if same_mesh else LinearNDInterpolator(w['points'],values)(m.p.T)
        vals=np.nan_to_num(vals);U=vals[:,:3].T.copy();U[:,velD]=u_bc[:,velD];theta=np.maximum(vals[:,4],0);theta[dn]=0
        sol=np.r_[U.ravel(),vals[:,3]]
        if same_mesh and 'checkpoint_iteration' in w.files:
            if str(w['checkpoint_case'])!=case or float(w['checkpoint_nu_t'])!=nu_t:
                raise ValueError('Checkpoint case or eddy viscosity does not match')
            if str(w['checkpoint_config'])!=json.dumps(CFG,sort_keys=True,ensure_ascii=False):
                raise ValueError('Checkpoint configuration does not match')
            first_iteration=int(w['checkpoint_iteration'])
            checkpoint_history=w['checkpoint_history'].tolist()
            elapsed_before=float(checkpoint_history[-1][-1]) if checkpoint_history else 0.
    fbw=FacetBasis(m,ElementTetP1(),facets=f['walls'],intorder=3)
    fbo=FacetBasis(m,ElementTetP1(),facets=np.r_[f['left'],f['right']],intorder=3)
    fbt=FacetBasis(m,ElementTetP1(),facets=np.r_[f['hood'],f['left'],f['right']],intorder=3)
    wallmat=asm(mass,fbw)*P['wall_heat_transfer_W_m2K']/(P['rho']*P['cp'])
    coords=bs.global_coordinates();source=np.zeros_like(coords[0]);sources=[]
    for bx,by,bz in CFG['burners']:
        s=np.exp(-((coords[0]-bx)**2+(coords[1]-by)**2)/(2*.085**2)-((coords[2]-(bz+.07))**2)/(2*.055**2))
        s*=((coords[2]>bz)&(coords[2]<bz+.25)&(coords[1]<CD))
        s*=P['sensible_heat_each_W']/np.sum(s*bs.dx)
        sources.append(s);source+=s
    source/=P['rho']*P['cp']
    t_free=np.setdiff1d(np.arange(n),dn)
    history=checkpoint_history;stable=0
    GD=[[asm(grad_div,bs,i=i,j=j) for j in range(3)] for i in range(3)]
    print(f'START {case}{suffix}: nodes={n} tets={m.nelements} Q={q*60:.5f} m3/min',flush=True)
    hq=np.asarray(bs.mesh_parameters())
    for it in range(first_iteration,maxiter):
        oldU=U.copy();oldT=theta.copy()
        vel=np.array([bs.interpolate(ui) for ui in U]);speed=np.sqrt(np.sum(vel**2,axis=0))
        invdt=10. if it<20 else (4. if it<50 else (1. if it<80 else .5))
        tau=1/np.sqrt((2*speed/hq)**2+(12*nu/hq**2)**2)
        force=P['gravity']*bs.interpolate(theta)/(T0+273.15+bs.interpolate(theta))
        kw=dict(vel=vel,tau=tau,nu=nu,invdt=invdt)
        uf=np.array([fbo.interpolate(ui) for ui in U]);vn=np.sum(uf*fbo.normals,axis=0)
        A=asm(transport,bs,**kw)+asm(backflow,fbo,penalty=.5*np.maximum(-vn,0))
        G=[asm(pressure_term,bs,**kw,i=i) for i in range(3)]
        B=[asm(continuity_term,bs,**kw,i=i) for i in range(3)]
        S=asm(p_stab,bs,**kw)
        AA=bmat([[A+GD[0][0],GD[0][1],GD[0][2],G[0]],[GD[1][0],A+GD[1][1],GD[1][2],G[1]],[GD[2][0],GD[2][1],A+GD[2][2],G[2]],[B[0],B[1],B[2],S]],format='csc')
        fvec=invdt*vel.copy();fvec[2]+=force
        rhs=np.r_[*[asm(body_force,bs,force=fvec[i],**kw) for i in range(3)],asm(p_force,bs,fvec=fvec,**kw)]
        new=np.r_[u_bc.ravel(),np.zeros(n)]
        rhsf=rhs[free]-AA[free][:,D]@new[D]
        new[free]=spsolve(AA[free][:,free],rhsf)
        if not np.all(np.isfinite(new)):raise RuntimeError('Nonfinite FEM solution')
        linres=np.linalg.norm((AA@new-rhs)[free])/(np.linalg.norm(rhsf)+1e-15)
        relax=.85
        U=(1-relax)*U+relax*new[:3*n].reshape(3,n)
        sol[:3*n]=U.ravel();sol[3*n:]=(1-relax)*sol[3*n:]+relax*new[3*n:]
        vel=np.array([bs.interpolate(ui) for ui in U]);speed=np.sqrt(np.sum(vel**2,axis=0))
        taut=1/np.sqrt((2*speed/hq)**2+(12*alpha/hq**2)**2)
        AT=thermal_matrix(bs,fbt,U,alpha,wallmat)
        rhst=asm(positive_heat_load,bs,source=source)
        newt=np.zeros(n);newt[t_free]=spsolve(AT[t_free][:,t_free].tocsc(),rhst[t_free])
        theta=.5*theta+.5*newt
        du=np.linalg.norm(U-oldU)/(np.linalg.norm(U)+1e-14)
        dt=np.linalg.norm(theta-oldT)/(np.linalg.norm(theta)+1e-14)
        history.append([it+1,du,dt,linres,float(np.min(theta+T0)),float(np.max(theta+T0)),elapsed_before+time.time()-start])
        if (it+1)%10==0:
            np.savez_compressed(ROOT/'results'/(case+suffix+'_checkpoint.npz'),
                points=m.p.T,cells=m.t.T,velocity=U.T,pressure=sol[3*n:]*P['rho'],temperature=theta+T0,
                checkpoint_iteration=it+1,checkpoint_history=np.array(history),checkpoint_case=case,
                checkpoint_nu_t=nu_t,checkpoint_config=json.dumps(CFG,sort_keys=True,ensure_ascii=False))
        if (it+1)%5==0 or it==0:
            print(f'{case}{suffix} {it+1}: du={du:.3g} dT={dt:.3g} Tmax={np.max(theta+T0):.2f}C elapsed={time.time()-start:.1f}s',flush=True)
        if du<2e-5 and dt<2e-5:stable+=1
        else:stable=0
        if stable>=3:break
    # Final nonlinear residual, including retained nodal pressure from relaxation.
    vel=np.array([bs.interpolate(ui) for ui in U]);speed=np.sqrt(np.sum(vel**2,axis=0))
    tau=1/np.sqrt((2*speed/hq)**2+(12*nu/hq**2)**2);kw=dict(vel=vel,tau=tau,nu=nu,invdt=0.)
    force=P['gravity']*bs.interpolate(theta)/(T0+273.15+bs.interpolate(theta))
    uf=np.array([fbo.interpolate(ui) for ui in U]);vn=np.sum(uf*fbo.normals,axis=0)
    A=asm(transport,bs,**kw)+asm(backflow,fbo,penalty=.5*np.maximum(-vn,0));G=[asm(pressure_term,bs,**kw,i=i) for i in range(3)];B=[asm(continuity_term,bs,**kw,i=i) for i in range(3)];S=asm(p_stab,bs,**kw)
    AA=bmat([[A+GD[0][0],GD[0][1],GD[0][2],G[0]],[GD[1][0],A+GD[1][1],GD[1][2],G[1]],[GD[2][0],GD[2][1],A+GD[2][2],G[2]],[B[0],B[1],B[2],S]],format='csr')
    fvec=np.zeros_like(vel);fvec[2]=force
    rhs=np.r_[np.zeros(2*n),asm(body_force,bs,force=force,**kw),asm(p_force,bs,fvec=fvec,**kw)]
    nonlin=np.linalg.norm((AA@sol-rhs)[free])/(np.linalg.norm(rhs[free]-AA[free][:,D]@sol[D])+1e-15)
    flux={key:facet_flux(m,ff,U) for key,ff in f.items()}
    total_flux=sum(v[0] for v in flux.values())
    # Exact integrated FEM volume, averages, gradient norms.
    vol=float(np.sum(bs.dx));Tq=bs.interpolate(theta)+T0;uq=np.array([bs.interpolate(ui) for ui in U]);spd=np.sqrt(np.sum(uq**2,axis=0))
    gradients=np.array([bs.interpolate(ui).grad for ui in U]);div=sum(gradients[i,i] for i in range(3))
    div_norm=float(np.sqrt(np.sum(div**2*bs.dx)));grad_norm=float(np.sqrt(np.sum(gradients**2*bs.dx[None,None])))
    allpq=bs.interpolate(sol[3*n:]*P['rho']);occ=(coords[1]>CD)&(coords[2]>CZ)&(coords[2]<CZ+.9)
    def avg(a,mask=None):
        wt=bs.dx if mask is None else bs.dx*mask
        return float(np.sum(a*wt)/np.sum(wt))
    # Advective energy outflow; wall loss and door conductive exchange are also exported.
    energy_flux={}
    for name in ['hood','left','right']:
        if len(f[name])==0:energy_flux[name]=0.;continue
        fb=FacetBasis(m,ElementTetP1(),facets=f[name],intorder=3)
        uf=np.array([fb.interpolate(U[i]) for i in range(3)]);vn=np.sum(uf*fb.normals,axis=0)
        tf=fb.interpolate(theta)
        energy_flux[name]=float(P['rho']*P['cp']*np.sum(vn*tf*fb.dx))
    wall_loss=float(P['wall_heat_transfer_W_m2K']*np.sum(fbw.interpolate(theta)*fbw.dx))
    # The fixed ambient-temperature reaction captures physical and algebraic diffusive exchange at doors.
    vel=np.array([bs.interpolate(ui) for ui in U]);speed=np.sqrt(np.sum(vel**2,axis=0));taut=1/np.sqrt((2*speed/hq)**2+(12*alpha/hq**2)**2)
    AT=thermal_matrix(bs,fbt,U,alpha,wallmat);rt=asm(positive_heat_load,bs,source=source)
    thermal_res=AT@theta-rt;thermal_free_res=float(np.linalg.norm(thermal_res[t_free])/(np.linalg.norm(rt[t_free])+1e-15))
    inlet_reaction=float(-np.sum(thermal_res[dn])*P['rho']*P['cp'])
    # Diagnostic only: discrepancy that a nonconservative heat form would introduce. Not subtracted from the conservative balance.
    div_heat=float(P['rho']*P['cp']*np.sum(bs.interpolate(theta)*div*bs.dx))
    heat_closure=sum(energy_flux.values())+wall_loss+inlet_reaction-HEAT
    metric={'layout':CFG['layout'],'layout_label':CFG['layout_label'],'buoyancy_model':'nonlinear ideal-gas density correction with fixed-reference-density incompressible equations','case':case,'label':CFG['case_labels'][case], 'mesh_target_m':h,'nodes':int(n),'tetrahedra':int(m.nelements),'fluid_volume_m3':vol,'nu_t_m2_s':nu_t,'iterations':len(history),'converged':bool(stable>=3),'last_relative_velocity_update':du,'last_relative_temperature_update':dt,'nonlinear_relative_residual':float(nonlin),'thermal_relative_residual':thermal_free_res,'last_linear_relative_residual':linres,'runtime_s':elapsed_before+time.time()-start,'fan_flow_m3_s':q,'fan_flow_m3_min':q*60,'network_makeup_pressure_drop_Pa':dp,'physical_makeup_area_m2':area,'boundary_flux_m3_s':{k:{'net_outward':v[0],'gross_inward':v[1],'gross_outward':v[2]} for k,v in flux.items()},'net_flux_error_m3_s':total_flux,'net_flux_relative_error':abs(total_flux)/q,'divergence_L2_s_inv_m_3_2':div_norm,'divergence_gradient_norm_ratio':div_norm/(grad_norm+1e-15),'max_nodal_speed_m_s':float(np.linalg.norm(U,axis=0).max()),'mean_speed_m_s':avg(spd),'occupied_zone_mean_speed_m_s':avg(spd,occ),'mean_temperature_C':avg(Tq),'max_nodal_temperature_C':float(theta.max()+T0),'min_nodal_temperature_C':float(theta.min()+T0),'occupied_zone_mean_temperature_C':avg(Tq,occ),'pressure_reference':'Natural traction at modeled room-side makeup boundary; leakage network pressure drop reported separately.','FE_pressure_min_Pa':float(allpq.min()),'FE_pressure_max_Pa':float(allpq.max()),'heat_input_W':HEAT,'heat_formulation':'conservative_galerkin_edge_diffusion','advective_heat_flux_W':energy_flux,'wall_heat_loss_W':wall_loss,'door_thermal_reaction_W':inlet_reaction,'weak_divergence_heat_correction_W':0.,'nonconservative_heat_divergence_diagnostic_W':div_heat,'discrete_heat_balance_residual_W':float(heat_closure)}
    out=ROOT/'results';name=case+suffix
    np.savez_compressed(out/(name+'.npz'),points=m.p.T,cells=m.t.T,velocity=U.T,pressure=sol[3*n:]*P['rho'],temperature=theta+T0,boundary_hood=f['hood'],boundary_left=f['left'],boundary_right=f['right'])
    meshio.write(out/(name+'.vtu'),meshio.Mesh(m.p.T,[('tetra',m.t.T)],point_data={'velocity_m_s':U.T,'speed_m_s':np.linalg.norm(U,axis=0),'pressure_room_reference_Pa':sol[3*n:]*P['rho'],'temperature_C':theta+T0}))
    np.savetxt(out/(name+'_nodes.csv'),np.c_[m.p.T,U.T,sol[3*n:]*P['rho'],theta+T0],delimiter=',',header='x_m,y_m,z_m,ux_m_s,uy_m_s,uz_m_s,pressure_room_reference_Pa,temperature_C',comments='')
    np.savetxt(out/(name+'_convergence.csv'),np.array(history),delimiter=',',header='iteration,relative_velocity_update,relative_temperature_update,linear_relative_residual,min_temperature_C,max_temperature_C,elapsed_s',comments='')
    (out/(name+'_metrics.json')).write_text(json.dumps(metric,ensure_ascii=False,indent=2))
    summaries={p.stem.removesuffix('_metrics'):json.loads(p.read_text()) for p in out.glob('*_metrics.json')}
    (out/'metrics.json').write_text(json.dumps(summaries,ensure_ascii=False,indent=2))
    print('DONE '+json.dumps(metric,ensure_ascii=False),flush=True)
    return metric

if __name__=='__main__':
    a=argparse.ArgumentParser();a.add_argument('--case',default='all');a.add_argument('--h',type=float,default=.18);a.add_argument('--nu-t',type=float);a.add_argument('--suffix',default='');a.add_argument('--warm');a.add_argument('--maxiter',type=int,default=300);args=a.parse_args()
    for case in CFG['cases'] if args.case=='all' else [args.case]:run(case,args.h,args.nu_t,args.suffix,args.maxiter,args.warm)
