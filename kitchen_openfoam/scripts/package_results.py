"""Archive reviewable results and reproducible native cases with SHA-256 checks."""
from pathlib import Path
import json,hashlib,zipfile,re,shutil
from foam_io import times,time_dir
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'deliverables'

def bundle(name,files,base):
    files=sorted(set(Path(p) for p in files if Path(p).is_file()));manifest={str(p.relative_to(base)):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    target=OUT/name
    with zipfile.ZipFile(target,'w',compression=zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for p in files:z.write(p,str(p.relative_to(base)))
        z.writestr('SHA256_MANIFEST.json',json.dumps(manifest,ensure_ascii=False,indent=2))
    with zipfile.ZipFile(target) as z:
        if z.testzip() is not None:raise RuntimeError('ZIP CRC check failed')
        for n,h in manifest.items():
            if hashlib.sha256(z.read(n)).hexdigest()!=h:raise RuntimeError(f'ZIP hash failed: {n}')
    return {'file':target.name,'files':len(files),'bytes':target.stat().st_size,'sha256':hashlib.sha256(target.read_bytes()).hexdigest(),'verified':True}

def main():
    evidence=OUT/'verification';evidence.mkdir(exist_ok=True)
    for name in ['final_validation.json','particle_physics_checks.json','mesh_inventory.json','environment_versions.json','artifact_review.json','checkMesh_corner_primary.log','checkMesh_center_primary.log','checkMesh_corner_graded.log','checkMesh_center_graded.log','checkMesh_corner_primary_basic.log','checkMesh_center_primary_basic.log','checkMesh_corner_graded_basic.log','checkMesh_center_graded_basic.log']:
        p=ROOT/'verification'/name
        if p.exists():shutil.copy2(p,evidence/name)
    refs=OUT/'references';refs.mkdir(exist_ok=True)
    shutil.copy2(ROOT/'references/SOURCES_zh.md',refs/'SOURCES_zh.md')
    shutil.copy2(ROOT/'emission_scenario.json',OUT/'emission_scenario.json')
    results=[p for p in OUT.rglob('*') if p.is_file() and p.suffix.lower() not in ['.zip'] and p.name!='archive_verification.json']
    result=bundle('openfoam_results.zip',results,OUT)
    source=[ROOT/'README_zh.md',ROOT/'PROGRESS.json',ROOT/'PROGRESS_zh.md',ROOT/'COPYING',ROOT/'emission_scenario.json',ROOT/'runtime_timestep_events.json',ROOT/'bin/kitchenAerosolFoam',ROOT/'requirements-analysis.txt',ROOT/'verification/environment_versions.json',ROOT/'verification/mesh_inventory.json']
    source+=[p for p in (ROOT/'scripts').glob('*.py') if p.name not in {'wait_corner_grid.py','generate_mesh.py'}]
    source+=[p for p in (ROOT/'src/kitchenAerosolFoam').rglob('*') if p.is_file() and ('Make/linux' not in str(p)) and p.suffix not in ['.o','.dep']]
    source+=list((ROOT/'inputs').glob('*/config.json'))+[ROOT/'references/SOURCES_zh.md']
    source+=list((ROOT/'verification').glob('checkMesh*.log'))+[ROOT/'verification/particle_physics_checks.json',ROOT/'verification/final_validation.json']
    for rootname in ['cases','sensitivity']:
        for case in (ROOT/rootname).iterdir():
            if not (case/'mesh_reference.npz').exists():continue
            (case/'kitchen.foam').touch(exist_ok=True)
            source+=list((case/'constant').rglob('*'))+list((case/'system').glob('*'))+list((case/'diagnostics').glob('*.csv'))
            source += [case/'mesh_reference.npz',case/'case_manifest.json',case/'run_status.json',case/'sensitivity_manifest.json',case/'kitchen.foam']
            ts=times(case)
            if not ts:continue
            selected={ts[0],ts[-1]}
            if rootname=='cases':
                selected.add(min(ts,key=lambda t:abs(t-100)))
            before=[t for t in ts if 110<=t<=120+1e-7]
            if before:selected.add(max(before))
            for t in selected:source+=list(time_dir(case,t).glob('*'))
    native=bundle('openfoam_cases_sources.zip',source,ROOT)
    (OUT/'archive_verification.json').write_text(json.dumps({'archives':[result,native],'native_checkpoint_selection':'initial, around100s, latest cooking state at/before120s, final; all per-step histories retained','excluded':'rejected pilot fields and unused steady solver; original workspace keeps them for audit'},ensure_ascii=False,indent=2))
    print(json.dumps([result,native],indent=2))

if __name__=='__main__':main()
