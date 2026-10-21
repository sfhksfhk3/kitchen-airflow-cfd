"""Reuse immutable original mesh/flow frames; no source mesh generation."""
from pathlib import Path
import sys,json,re,shutil,hashlib,argparse,numpy as np
ROOT=Path(__file__).resolve().parents[1]
FIELDS=['U','T','p','k','nut','phi']
NAMES=['Yp01','Yp03','Yp1','Yp3','Yp10']
def header(name,cls='dictionary'):
 return f'FoamFile\n{{ version 2.0; format ascii; class {cls}; object {name}; }}\n'
def prepare(source,name,start=0,end=240,dt=.1,stride=10,**options):
 source=Path(source).resolve();case=ROOT/'runs'/name
 if case.exists():raise FileExistsError(case)
 case.mkdir(parents=True);(case/'system').mkdir();(case/'constant').mkdir();(case/'frames').mkdir()
 shutil.copytree(source/'constant/polyMesh',case/'constant/polyMesh')
 frames=sorted((float(p.name),p) for p in source.iterdir() if p.is_dir() and re.fullmatch(r'\d+(?:\.\d+)?',p.name) and all((p/f).exists() for f in FIELDS if f!='phi') and ((p/'phi').exists() or float(p.name)==0))
 before=[i for i,(t,p) in enumerate(frames) if t<=start+1e-6]
 after=[i for i,(t,p) in enumerate(frames) if t>=end-1e-6]
 assert before and after,(source,start,end)
 frames=frames[max(before):min(after)+1]
 thin=max(1,round(stride/float(np.median(np.diff([t for t,p in frames])))));frames=[x for i,x in enumerate(frames) if i%thin==0 or i==len(frames)-1]
 assert len(frames)>1 and frames[0][0]<=start+1e-6 and frames[-1][0]>=end-1e-6,(source,start,end,[t for t,p in frames])
 for t,p in frames:
  target=case/'frames'/p.name
  if (p/'phi').exists():target.symlink_to(p,target_is_directory=True)
  else:
   target.mkdir()
   for f in FIELDS:
    if f!='phi':(target/f).symlink_to(p/f)
   b=(source/'constant/polyMesh/boundary').read_text()
   patches=re.findall(r'(\w+)\s*\{[^{}]*nFaces',b)
   body='dimensions [1 0 -1 0 0 0 0];\ninternalField uniform 0;\nboundaryField {\n'+''.join(f'{pn} {{ type calculated; value uniform 0; }}\n' for pn in patches)+'}\n'
   (target/'phi').write_text(header('phi','surfaceScalarField')+body)
 init=case/f'{start:.12g}';init.mkdir()
 for n in NAMES:shutil.copy2(frames[0][1]/n,init/n)
 shape=next(p/'sourceShape' for _,p in frames if (p/'sourceShape').exists())
 shutil.copy2(shape,init/'sourceShape')
 # sourceShape is not always saved in mapped sensitivity cases; normalization remains original.
 shutil.copy2(source/'mesh_reference.npz',case/'mesh_reference.npz')
 text=(source/'system/fvSchemes').read_text();(case/'system/fvSchemes').write_text(text)
 (case/'system/fvSolution').write_text(header('fvSolution')+'''solvers {
 "Yp.*" { solver PBiCGStab; preconditioner DILU; tolerance 1e-20; relTol 1e-10; }
 }
''')
 (case/'system/controlDict').write_text(header('controlDict')+f'''application ceilingReplayFoam;
startFrom startTime; startTime {start}; stopAt endTime; endTime {end}; deltaT {dt};
writeControl runTime; writeInterval 20; purgeWrite 0;
writeFormat ascii; writePrecision 14; writeCompression off;
timeFormat general; timePrecision 12; runTimeModifiable false;
''')
 params=dict(ceilingSticking=1,otherSticking=1,SchmidtT=.7,inertiaFactor=1,turbophoresisFactor=0,impactionFactor=1,impactionMode="resistance",thermophoresisFactor=1,wallTemperatureOffsetK=0,coefficientUpdateSeconds=1,emissionKgPerSecPerBin=1e-9,sourceStopTime=120,wallQuadrature=72,legacyWallModel=False,useSavedInitialAerosols=False)
 params.update(options)
 body='frameTimes ('+' '.join(str(t) for t,p in frames)+');\nframeDirectories ('+' '.join('"'+p.name+'"' for t,p in frames)+');\n'
 body+='fieldNames ('+' '.join(NAMES)+');\ndiameters (1e-7 3e-7 1e-6 3e-6 1e-5);\n'
 for k,v in params.items():body+=f'{k} {str(v).lower() if isinstance(v,bool) else v};\n'
 (case/'constant/ceilingProperties').write_text(header('ceilingProperties')+body)
 mesh_hash={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (case/'constant/polyMesh').iterdir() if p.is_file()}
 manifest=dict(source=str(source),name=name,start_s=start,end_s=end,dt_s=dt,frame_stride_s=stride,frame_times_s=[t for t,p in frames],parameters=params,source_mesh_sha256=mesh_hash,scope='conservative aerosol replay; original coupled air flow not recomputed; native mass flux linearly interpolated')
 from foam_mesh import write_ceiling_set
 manifest['ceiling_face_count'],manifest['ceiling_area_m2']=write_ceiling_set(case)
 (case/'replay_manifest.json').write_text(json.dumps(manifest,indent=2));(case/'ceiling.foam').touch()
 return case
if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('source');a.add_argument('name');a.add_argument('--start',type=float,default=0);a.add_argument('--end',type=float,default=240);a.add_argument('--dt',type=float,default=.1);a.add_argument('--stride',type=float,default=10);a.add_argument('--options',default='{}');v=a.parse_args();print(prepare(v.source,v.name,v.start,v.end,v.dt,v.stride,**json.loads(v.options)))
