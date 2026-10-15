"""Continue matched saved states for timestep and late-stage mesh sensitivity.
Does not overwrite completed primary simulations. Positive interpolation plus
an explicit per-bin mass correction is used only for mesh transfer.
"""
from pathlib import Path
import argparse,json,re,shutil,numpy as np
from scipy.spatial import cKDTree
from foam_io import times,time_dir,load_mesh,read_internal,fields,FOAM_RR,AIR_MOLAR_MASS
from build_cases import build
ROOT=Path(__file__).resolve().parents[1]

def replace_internal(path,values):
    values=np.asarray(values);kind='vector' if values.ndim==2 else 'scalar'
    lines=['('+' '.join(f'{x:.14g}' for x in v)+')' for v in values] if kind=='vector' else [f'{x:.14g}' for x in values]
    new=f'internalField nonuniform List<{kind}> {len(values)} (\n'+'\n'.join(lines)+'\n);'
    s=re.sub(r'internalField\s+.*?;',lambda _:new,path.read_text(),count=1,flags=re.S);path.write_text(s)

def prepare(layout,doors,kind):
    src=ROOT/'cases'/f'{layout}_{doors}';pool=[t for t in times(src) if 90<=t<=110]
    if not pool:raise RuntimeError('Need a saved primary state around 100s')
    start=min(pool,key=lambda v:abs(v-100));sd=time_dir(src,start)
    dst=ROOT/'sensitivity'/f'{layout}_{doors}_{kind}'
    if dst.exists():print('EXISTS',dst);return dst
    end=140.
    if kind=='dt':
        dst.mkdir(parents=True)
        for name in ['constant','system']:shutil.copytree(src/name,dst/name)
        shutil.copytree(sd,dst/sd.name);shutil.copy2(src/'mesh_reference.npz',dst/'mesh_reference.npz');shutil.copy2(src/'case_manifest.json',dst/'case_manifest.json')
        c=dst/'system/controlDict';s=c.read_text();s=re.sub(r'maxCo\s+[\d.]+;','maxCo 3;',s);s=re.sub(r'maxDeltaT\s+[\d.]+;','maxDeltaT 0.02;',s);c.write_text(s)
        transfer={'method':'identical saved fields and mesh; no interpolation','mass_relative_transfer_error':[0]*5}
    elif kind=='grid':
        build(layout,doors,h=.08,end=end,case_dir=dst,gap_layers=4,maxdt=.04)
        old=load_mesh(src);new=load_mesh(dst);oldf=fields(src,start,old)
        dist,ix=cKDTree(old['centres']).query(new['centres'],k=4)
        ww=1/np.maximum(dist,1e-12)**2;ww/=ww.sum(1)[:,None]
        mapped={}
        for name in ['U','T','p','p_rgh','k','omega','nut','alphat','Yp01','Yp03','Yp1','Yp3','Yp10']:
            a=read_internal(sd/name,len(old['cells']));mapped[name]=np.einsum('ij,ijk->ik',ww,a[ix]) if a.ndim==2 else np.sum(ww*a[ix],axis=1)
        R=FOAM_RR/AIR_MOLAR_MASS;rn=mapped['p']/(R*mapped['T']);errs=[];factors=[]
        for name in ['Yp01','Yp03','Yp1','Yp3','Yp10']:
            mo=np.sum(oldf['rho']*oldf[name]*old['volumes']);mn=np.sum(rn*mapped[name]*new['volumes']);factor=mo/mn;mapped[name]*=factor;factors.append(float(factor));errs.append(float((np.sum(rn*mapped[name]*new['volumes'])-mo)/mo))
        for name,val in mapped.items():replace_internal(dst/'0'/name,val)
        shutil.copy2(sd/'aerosolState',dst/'0/aerosolState');(dst/'0').rename(dst/sd.name)
        transfer={'method':'positive 4-neighbor interpolation; per-size total-mass rescale; restart flow re-equilibrates','source_cells':len(old['cells']),'target_cells':len(new['cells']),'source_gap_layers':2,'target_gap_layers':4,'mass_rescale_factors':factors,'mass_relative_transfer_error':errs,'limitation':'late-stage matched-state grid sensitivity, not independent full-history grid convergence'}
        c=dst/'system/controlDict';c.write_text(c.read_text().replace('maxCo 3;','maxCo 12;'))
        q=dst/'constant/aerosolProperties';q.write_text(q.read_text()+'\nprojectInitialMassFlux yes;\n')
        q=dst/'system/fvSolution';q.write_text(q.read_text().replace('solvers {','solvers {\n fluxCorrection { solver GAMG; tolerance 1e-12; relTol 0; smoother GaussSeidel; }'))
        transfer['grading']='Near-floor cell heights grow through 1.25,2.5,5,10,20mm layers before bulk mesh'
    c=dst/'system/controlDict';s=re.sub(r'endTime\s+[\d.]+;',f'endTime {end};',c.read_text());s=s.replace('stopAt writeNow','stopAt endTime');c.write_text(s)
    meta={'source_case':str(src),'source_time_s':start,'target_end_s':end,'kind':kind,'transfer':transfer}
    (dst/'sensitivity_manifest.json').write_text(json.dumps(meta,indent=2));print('PREPARED',dst,start,flush=True);return dst

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--layout',default='center');ap.add_argument('--doors',default='both_closed');ap.add_argument('--kind',choices=['dt','grid'],required=True);a=ap.parse_args();prepare(a.layout,a.doors,a.kind)
