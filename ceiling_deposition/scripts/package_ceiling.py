"""Separate readable/reproducible output from large native case data."""
from pathlib import Path
import json,hashlib,zipfile,tarfile,io,datetime,subprocess,argparse
ROOT=Path(__file__).resolve().parents[1]
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  while b:=f.read(1024*1024):h.update(b)
 return h.hexdigest()
def main(out):
 out.mkdir(exist_ok=True,parents=True);commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip()
 note={'status':'research preview; ceiling grid robustness NOT passed','source_commit':commit,'source_branch':'feature/ceiling-grease-deposition','original_release_unchanged':'https://github.com/sfhksfhk3/kitchen-airflow-cfd/releases/tag/v1.0.0','original_flow_snapshot_sha256':'608f2eb30b22f97b8ed542a9f66eecefc3499d8ea66edc32eb0b314fb35fe1f6','created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
 review=out/'ceiling_deposition_review_bundle.zip';files=[]
 excluded={'runs','flow_sensitivity','tail_sources','bin','data','__pycache__'}
 for p in ROOT.rglob('*'):
  rel=p.relative_to(ROOT)
  if not p.is_file() or p.is_symlink():continue
  if any(part in excluded for part in rel.parts):continue
  if any(part.startswith('linux') for part in rel.parts) or p.suffix in ['.log','.pyc','.o','.dep']:continue
  files.append(p)
 manifest={str(p.relative_to(ROOT.parent)):{'bytes':p.stat().st_size,'sha256':sha(p)} for p in files}
 with zipfile.ZipFile(review,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
  for p in files:z.write(p,p.relative_to(ROOT.parent))
  z.writestr('ceiling_deposition/REVIEW_MANIFEST.json',json.dumps({'publication':note,'files':manifest},ensure_ascii=False,indent=2))
 with zipfile.ZipFile(review) as z:
  assert z.testzip() is None
  for name,m in manifest.items():assert hashlib.sha256(z.read(name)).hexdigest()==m['sha256']
 print('REVIEW_BUNDLE_VERIFIED',review.stat().st_size,len(files),flush=True)
 native=out/'ceiling_deposition_native_cases.tar.gz'
 cases=[p for p in (ROOT/'runs').iterdir() if (p/'RUN_COMPLETE.json').exists() or '_clearance_' in p.name or p.name=='pilot_final']
 sources=cases+[ROOT/'flow_sensitivity',ROOT/'tail_sources',ROOT/'bin']
 logs=[p for p in (ROOT/'verification').glob('*.log')];sources+=logs
 native_manifest={};symlinks={}
 for source in sources:
  for p in ([source] if source.is_file() else source.rglob('*')):
   name=str(p.relative_to(ROOT.parent))
   if p.is_symlink():symlinks[name]=str(p.readlink());continue
   if p.is_file() and p.name not in ['core']:
    native_manifest[name]={'bytes':p.stat().st_size,'sha256':sha(p)}
 meta=json.dumps({'publication':note,'case_count':len(cases),'notes':['New computed cases only; old flow restored separately from v1.0.0.','Absolute read-only frame symlinks preserved. Use relink_frames.py if project or flow root changes.','Very large native tail-case timestamps are algebraic regularization parameters, NOT physical observation time.','Failed exploration logs retained for diagnostics, never as accepted scientific results.'],'files':native_manifest,'symlinks':symlinks},ensure_ascii=False,indent=2).encode()
 with tarfile.open(native,'w:gz',compresslevel=1,dereference=False) as tf:
  for source in sources:tf.add(source,arcname=str(source.relative_to(ROOT.parent)),filter=lambda info:None if Path(info.name).name=='core' else info)
  info=tarfile.TarInfo('ceiling_deposition/NATIVE_MANIFEST.json');info.size=len(meta);tf.addfile(info,io.BytesIO(meta))
 with tarfile.open(native) as tf:
  for m in tf:
   assert not m.name.startswith('/') and '..' not in Path(m.name).parts,m.name
   if m.isfile() and m.name in native_manifest:
    h=hashlib.sha256();f=tf.extractfile(m)
    while b:=f.read(1024*1024):h.update(b)
    assert h.hexdigest()==native_manifest[m.name]['sha256'],m.name
 print('NATIVE_BUNDLE_VERIFIED',native.stat().st_size,len(native_manifest),len(symlinks),flush=True)
 assets=[review,native]
 (out/'SHA256SUMS').write_text(''.join(f'{sha(p)}  {p.name}\n' for p in assets))
 (out/'DELIVERY_INDEX.json').write_text(json.dumps({**note,'assets':[{'name':p.name,'bytes':p.stat().st_size,'sha256':sha(p)} for p in assets],'native_case_count':len(cases),'review_file_count':len(files),'native_regular_files':len(native_manifest),'native_readonly_symlinks':len(symlinks),'review_and_native_contents_verified':True},ensure_ascii=False,indent=2))
 print('ALL_PACKAGING_COMPLETE',flush=True)
if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('--out',type=Path,default=ROOT.parents[1]/'ceiling-deliverables');v=a.parse_args();main(v.out)
