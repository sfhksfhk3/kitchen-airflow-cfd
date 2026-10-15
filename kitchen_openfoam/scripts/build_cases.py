"""Generate/reuse Gmsh extruded hex/prism mesh; write native OpenFOAM polyMesh.

Five aerosol fields represent independent monodisperse reference experiments,
not an assumed five-bin mass distribution. SI units throughout.
"""
from pathlib import Path
import argparse, json, math, numpy as np

ROOT=Path(__file__).resolve().parents[1]
OLD=ROOT/'inputs'
DIAMS=[.1,.3,1.,3.,10.]
NAMES=['Yp01','Yp03','Yp1','Yp3','Yp10']

def header(name,cls='dictionary'):
    return f'FoamFile\n{{ version 2.0; format ascii; class {cls}; object {name}; }}\n'

def write(path,body,cls='dictionary'):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(header(path.name,cls)+body)

def make_mesh(case,layout,doors,h,gap_layers=2):
    from generate_extruded_mesh import generate
    a=np.load(generate(layout,h,gap_layers))
    points,tets=a['points'],a['cells'];centres=a['centres'];vol=a['volumes']
    nverts=a['nverts']
    faces={}
    for ci,t in enumerate(tets):
        n=int(nverts[ci])//2
        faceids=[tuple(range(n)),tuple(range(n,2*n))]+[(i,(i+1)%n,(i+1)%n+n,i+n) for i in range(n)]
        for ids in faceids:
            f=[int(t[i]) for i in ids];key=tuple(sorted(f))
            if key in faces: faces[key][2]=ci
            else:
                q=points[f]
                if np.dot(np.cross(q[1]-q[0],q[2]-q[0]),q.mean(0)-centres[ci])<0: f.reverse()
                faces[key]=[f,ci,-1]
    c=json.loads((OLD/layout/'config.json').read_text());hd=c['hood']
    interior=[];patches={k:[] for k in ['hood','left','right','cabinetSide','walls']}
    baffle_faces=0
    for f,ow,ne in faces.values():
        fp=points[f];mid=fp.mean(0)
        is_baffle=(hd['x0']-1e-7<mid[0]<hd['x1']+1e-7 and .3<mid[1]<.3+.1/math.sqrt(2) and np.max(np.abs(fp[:,1]+fp[:,2]-1.9))<1e-7)
        if ne>=0:
            if is_baffle:
                patches['walls'].append((f,ow,-1));patches['walls'].append((f[::-1],ne,-1));baffle_faces+=1
            else:interior.append((f,ow,ne))
            continue
        x,y,z=points[f].mean(0);pn='walls'
        if abs(z-1.6)<1e-7 and hd['x0']<x<hd['x1'] and 0<y<.3: pn='hood'
        elif abs(x-.7)<1e-7 and y<.6: pn='cabinetSide'
        elif .7<y<1.6:
            if abs(x)<1e-7 and (doors in ['both_open','left_open'] or z<.005):pn='left'
            if abs(x-3.3)<1e-7 and (doors in ['both_open','right_open'] or z<.005):pn='right'
        patches[pn].append((f,ow,ne))
    interior.sort(key=lambda f:(f[1],f[2]))
    ordered=interior+[f for ff in patches.values() for f in ff]
    pm=case/'constant/polyMesh';pm.mkdir(parents=True,exist_ok=True)
    write(pm/'points',f'{len(points)}\n(\n'+'\n'.join('('+ ' '.join(f'{v:.14g}' for v in p)+')' for p in points)+'\n)\n','vectorField')
    write(pm/'faces',f'{len(ordered)}\n(\n'+'\n'.join(str(len(f[0]))+'('+' '.join(map(str,f[0]))+')' for f in ordered)+'\n)\n','faceList')
    write(pm/'owner',f'{len(ordered)}\n(\n'+'\n'.join(str(f[1]) for f in ordered)+'\n)\n','labelList')
    write(pm/'neighbour',f'{len(interior)}\n(\n'+'\n'.join(str(f[2]) for f in interior)+'\n)\n','labelList')
    b='5\n(\n';off=len(interior)
    areas={}
    for name,ff in patches.items():
        typ='wall' if name in ['walls','cabinetSide'] else 'patch'
        b+=f'{name}\n{{ type {typ}; nFaces {len(ff)}; startFace {off}; }}\n';off+=len(ff)
        areas[name]=float(sum(np.linalg.norm(np.sum(np.cross(points[f[0]],np.roll(points[f[0]],-1,axis=0)),axis=0))/2 for f in ff))
    write(pm/'boundary',b+')\n','polyBoundaryMesh')
    np.savez_compressed(case/'mesh_reference.npz',points=points,cells=tets,centres=centres,volumes=vol,nverts=nverts)
    assert baffle_faces>10, f'Embedded baffle not found: {baffle_faces}'
    hf=np.array([ff[0] for ff in patches['hood']]);hp=points[hf]
    hcentres=hp.mean(1);hareas=np.linalg.norm(np.sum(np.cross(hp,np.roll(hp,-1,axis=1)),axis=1),axis=1)/2
    return c,areas,(hcentres,hareas)

def field(case,name,dim,initial,bc,vector=False):
    body=f'dimensions [{dim}];\ninternalField uniform {initial};\nboundaryField\n{{\n'
    for pn in ['hood','left','right','cabinetSide','walls']:body+=f'{pn}\n{{ {bc[pn]} }}\n'
    write(case/'0'/name,body+'}\n','volVectorField' if vector else 'volScalarField')

def build(layout,doors,h=.1,end=240.,maxdt=.05,case_dir=None,gap_layers=2):
    case=Path(case_dir) if case_dir else ROOT/'cases'/f'{layout}_{doors}'
    if case.exists():
        from foam_io import times
        if any(t>0 for t in times(case)):raise RuntimeError(f'Existing computed fields: resume this case instead of rebuilding {case}')
    case.mkdir(parents=True,exist_ok=True)
    cfg,areas,(hcentres,hareas)=make_mesh(case,layout,doors,h,gap_layers)
    op=['both_open','left_open','right_open'];openL=doors in ['both_open','left_open'];openR=doors in ['both_open','right_open']
    areaL=.9*(2.2 if openL else .005);areaR=.9*(2.2 if openR else .005)
    rho=101325/(287.05*298.15);Q=math.sqrt(1300/(rho/(2*(.6*(areaL+areaR))**2)+1300/.5**2))
    # Flow is an explicitly assumed network operating point, not a measured fan curve.
    # The actual 5 mm closed-door apertures are now geometric boundary patches.
    ubc={p:'type noSlip;' for p in areas}
    hx=cfg['hood']['x0']
    hw=np.sin(np.pi*(hcentres[:,0]-hx)/.9)*np.sin(np.pi*hcentres[:,1]/.3)
    hw*=Q/np.sum(hw*hareas)
    ubc['hood']='type zeroGradient;'
    for pn,isopen,ar in [('left',openL,areaL),('right',openR,areaR)]:
        ubc[pn]='type pressureInletOutletVelocity; value uniform (0 0 0);'
    field(case,'U','0 1 -1 0 0 0 0','(0 0 0)',ubc,True)
    pbc={pn:'type calculated; value uniform 101325;' for pn in areas}
    field(case,'p','1 -1 -2 0 0 0 0','101325',pbc)
    pressure_path=case/'0/p'
    z=np.load(case/'mesh_reference.npz')['centres'][:,2]
    p0=101325*np.exp(-9.81*z/(287.05*298.15))
    pressure_path.write_text(pressure_path.read_text().replace('internalField uniform 101325;','internalField nonuniform List<scalar> '+str(len(p0))+' (\n'+'\n'.join(f'{v:.12g}' for v in p0)+'\n);'))
    prbc={pn:'type fixedFluxPressure; value uniform 101325;' for pn in areas}
    for pn,isopen in [('left',openL),('right',openR)]:
        prbc[pn]='type prghTotalHydrostaticPressure; value uniform 101325;'
    prbc['hood']='type fixedValue; value uniform 101323;'
    field(case,'p_rgh','1 -1 -2 0 0 0 0','101325',prbc)
    tbc={pn:'type externalWallHeatFluxTemperature; mode coefficient; Ta constant 298.15; h uniform 3; kappaMethod fluidThermo; value uniform 298.15;' for pn in areas}
    for pn in ['hood','left','right']:tbc[pn]='type inletOutlet; inletValue uniform 298.15; value uniform 298.15;'
    field(case,'T','0 0 0 1 0 0 0','298.15',tbc)
    for name,dim,initial,wtype in [('k','0 2 -2 0 0 0 0','0.001','kqRWallFunction'),('omega','0 0 -1 0 0 0 0','1','omegaWallFunction'),('nut','0 2 -1 0 0 0 0','0','nutkWallFunction'),('alphat','1 -1 -1 0 0 0 0','0','compressible::alphatWallFunction')]:
        bc={}
        for pn in areas:
            if pn in ['walls','cabinetSide']:bc[pn]=f'type {wtype}; '+('Prt 0.85; ' if name=='alphat' else '')+f'value uniform {initial};'
            elif name in ['nut','alphat']:bc[pn]=f'type calculated; value uniform {initial};'
            else:bc[pn]=f'type inletOutlet; inletValue uniform {initial}; value uniform {initial};'
        field(case,name,dim,initial,bc)
    for name in NAMES:
        bc={pn:('type zeroGradient;' if pn in ['walls','cabinetSide'] else 'type inletOutlet; phi aerosolPhi; inletValue uniform 0; value uniform 0;') for pn in areas}
        field(case,name,'0 0 0 0 0 0 0','0',bc)
    write(case/'constant/g','dimensions [0 1 -2 0 0 0 0]; value (0 0 -9.81);','uniformDimensionedVectorField')
    write(case/'constant/thermophysicalProperties','''thermoType { type heRhoThermo; mixture pureMixture; transport sutherland; thermo hConst; equationOfState perfectGas; specie specie; energy sensibleEnthalpy; }
mixture { specie { molWeight 28.965; } thermodynamics { Cp 1006; Hf 0; } transport { As 1.458e-6; Ts 110.4; } }
''')
    write(case/'constant/turbulenceProperties','simulationType RAS; RAS { RASModel kOmegaSST; turbulence on; printCoeffs on; }\n')
    write(case/'constant/fvOptions','buoyancyProduction { type buoyancyTurbSource; active yes; selectionMode all; }\n')
    write(case/'constant/radiationProperties','radiation off; radiationModel none;\n')
    write(case/'constant/dynamicMeshDict','dynamicFvMesh staticFvMesh;\n')
    pts=' '.join('('+' '.join(str(x) for x in p)+')' for p in cfg['burners'])
    write(case/'constant/aerosolProperties',f'''diameters (1e-7 3e-7 1e-6 3e-6 1e-5);
fieldNames (Yp01 Yp03 Yp1 Yp3 Yp10);
particleDensity 850;
SchmidtT 0.7;
thermophoreticCoefficient 0.55;
thermophoresisModel Talbot;
particleThermalConductivity 0.17;
hoodTargetQ {Q:.12g};
pressureControllerGain {2000 if doors=='both_closed' else 20};
initialHoodPressureDrop {150 if doors=='both_closed' else 2};
stickingProbability 1;
emissionKgPerSecPerBin 1e-9;
sourceStopTime 120;
heatStopTime 120;
heatEachW 2000;
burners ({pts});
sourceSigmaXY 0.085;
sourceSigmaZ 0.055;
''')
    write(case/'system/controlDict',f'''application kitchenAerosolFoam;
startFrom latestTime; startTime 0; stopAt endTime; endTime {end};
deltaT 0.002; adjustTimeStep yes; maxCo 3; maxDeltaT {maxdt};
writeControl adjustableRunTime; writeInterval 10; purgeWrite 0;
writeFormat ascii; writePrecision 10; writeCompression off;
timeFormat general; timePrecision 10; runTimeModifiable yes;
functions {{
 yplus {{ type yPlus; libs ("libfieldFunctionObjects.so"); executeControl writeTime; writeControl writeTime; }}
}}
''')
    write(case/'system/fvSchemes','''ddtSchemes { default Euler; }
gradSchemes { default cellLimited Gauss linear 1; }
divSchemes {
 default none;
 div(phi,U) Gauss upwind;
 div(phi,h) Gauss limitedLinear 1;
 div(phi,K) Gauss upwind;
 div(phi,k) Gauss upwind;
 div(phi,omega) Gauss upwind;
 div(aerosolPhi,Y) Gauss upwind;
 div(((rho*nuEff)*dev2(T(grad(U))))) Gauss linear;
}
laplacianSchemes { default Gauss linear limited 0.5; }
interpolationSchemes { default linear; }
snGradSchemes { default limited 0.5; }
wallDist { method meshWave; }
fluxRequired { default no; p_rgh; }
''')
    write(case/'system/fvSolution','''solvers {
 "rho.*" { solver PCG; preconditioner DIC; tolerance 1e-10; relTol 0; }
 p_rgh { solver GAMG; tolerance 1e-8; relTol 0.01; smoother GaussSeidel; cacheAgglomeration true; agglomerator faceAreaPair; nCellsInCoarsestLevel 50; }
 p_rghFinal { $p_rgh; relTol 0; }
 "(U|h|k|omega)" { solver PBiCGStab; preconditioner DILU; tolerance 1e-8; relTol 0.01; }
 "(U|h|k|omega)Final" { solver PBiCGStab; preconditioner DILU; tolerance 1e-9; relTol 0; }
 "Yp.*" { solver PBiCGStab; preconditioner DILU; tolerance 1e-18; relTol 1e-9; }
}
PIMPLE { momentumPredictor yes; nOuterCorrectors 2; nCorrectors 2; nNonOrthogonalCorrectors 0; pRefCell 0; pRefValue 101325; }
relaxationFactors { equations { ".*" 1; } }
''')
    write(case/'system/decomposeParDict','numberOfSubdomains 3; method scotch;\n')
    manifest={'layout':layout,'doors':doors,'mesh_h_m':h,'mesh_kind':'Gmsh extruded hexahedra and prisms','Q_assumed_m3_s':Q,'patch_area_m2':areas,'endTime_s':end,'sourceStopTime_s':120,'heatStopTime_s':120,'closed_gap_actual_m':.005,'closed_gap_computational_height_m':.005,'gap_cells_through_height':2,'oil_density_assumed_kg_m3':850,'independent_monodisperse_diameters_um':DIAMS,'reference_emission_each_independent_bin_kg_s':1e-9,'source_geometry_config':cfg,'notes':['finite volume method, not finite elements','dilute one-way particle coupling; two-way temperature-density-momentum coupling','fan operating flow remains a network assumption; actual 5mm gap geometry with two cells over height; gap-jet convergence needs checking','five bins are independent response functions; no measured size distribution','wall deposition closure and SST wall functions need mesh/model validation']}
    manifest['gap_cells_through_height']=gap_layers
    (case/'case_manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2))
    (case/'kitchen.foam').touch()
    print(case.name,Q,len(np.load(case/'mesh_reference.npz')['cells']),flush=True)

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--h',type=float,default=.1);ap.add_argument('--end',type=float,default=240);ap.add_argument('--layout',choices=['corner','center']);ap.add_argument('--doors',choices=['both_open','both_closed','left_open','right_open']);a=ap.parse_args()
    for layout in ([a.layout] if a.layout else ['corner','center']):
        for doors in ([a.doors] if a.doors else ['both_open','both_closed','left_open','right_open']):build(layout,doors,a.h,a.end)
