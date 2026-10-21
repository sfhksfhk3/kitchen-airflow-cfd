from pathlib import Path
import re,shutil,subprocess,json,concurrent.futures
ROOT=Path(__file__).resolve().parents[1]
def run(layout):
 src=ROOT/'flow_sensitivity'/f'{layout}_both_closed_Qtarget120';dst=ROOT/'flow_sensitivity'/f'{layout}_both_closed_control'
 if not dst.exists():
  dst.mkdir()
  for d in ['constant','system']:shutil.copytree(src/d,dst/d)
  init=json.loads((src/'intervention.json').read_text())['matched_initial_checkpoint_s'];idir=min((p for p in src.iterdir() if p.is_dir() and re.fullmatch(r'\d+(?:\.\d+)?',p.name)),key=lambda p:abs(float(p.name)-init));shutil.copytree(idir,dst/idir.name)
  for n in ['mesh_reference.npz','case_manifest.json']:shutil.copy2(src/n,dst/n)
  p=dst/'constant/aerosolProperties';old=json.loads((src/'intervention.json').read_text())['original_target_m3_s'];p.write_text(re.sub(r'hoodTargetQ\s+[^;]+;',f'hoodTargetQ {old:.14g};',p.read_text()))
 if (dst/'140/U').exists():return str(dst)
 with (ROOT/'verification'/f'flow_{layout}_matched_control.log').open('w') as f:
  subprocess.run(['docker','exec','kitchen-openfoam-2406','bash','-lc',f'/workspace/kitchen_openfoam/bin/kitchenAerosolFoam -case {dst}'],stdout=f,stderr=subprocess.STDOUT,check=True)
 return str(dst)
if __name__=='__main__':
 with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
  for c in pool.map(run,['corner','center']):print(c,flush=True)
