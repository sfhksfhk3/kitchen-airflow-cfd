"""Coupled-air intervention from saved 100 s checkpoint, preserving original cases."""
from pathlib import Path
import shutil,json,re
ROOT=Path(__file__).resolve().parents[1]
def prepare(layout):
 src=Path('/workspace/kitchen_openfoam/cases')/f'{layout}_both_closed';dst=ROOT/'flow_sensitivity'/f'{layout}_both_closed_Qtarget120'
 if dst.exists():return dst
 dst.mkdir(parents=True)
 checkpoint=min((p for p in src.iterdir() if p.is_dir() and re.fullmatch(r'\d+(?:\.\d+)?',p.name)),key=lambda p:abs(float(p.name)-100))
 for d in ['constant','system']:shutil.copytree(src/d,dst/d)
 shutil.copytree(checkpoint,dst/checkpoint.name)
 shutil.copy2(src/'mesh_reference.npz',dst/'mesh_reference.npz')
 shutil.copy2(src/'case_manifest.json',dst/'case_manifest.json')
 p=dst/'constant/aerosolProperties';text=p.read_text();old=float(re.search(r'hoodTargetQ\s+([\d.eE+-]+)',text)[1]);text=re.sub(r'hoodTargetQ\s+[^;]+;',f'hoodTargetQ {old*1.2:.14g};',text);p.write_text(text)
 p=dst/'system/controlDict';text=p.read_text();text=re.sub(r'endTime\s+[^;]+;', 'endTime 140;',text);text=re.sub(r'writeInterval\s+[^;]+;', 'writeInterval 5;',text);p.write_text(text)
 (dst/'intervention.json').write_text(json.dumps({'source':str(src),'matched_initial_checkpoint_s':float(checkpoint.name),'end_s':140,'target_flow_factor':1.2,'original_target_m3_s':old,'new_target_m3_s':old*1.2,'pressure_cap_Pa':1300,'note':'Actual achieved Q must be measured; 1300 Pa cap may prevent full 20% target increase. Matched 100–140 s test, not full meal.'},indent=2))
 return dst
def prepare_heat_off():
 src=Path('/workspace/kitchen_openfoam/cases/corner_both_closed');dst=ROOT/'flow_sensitivity/corner_both_closed_heatOff100'
 if dst.exists():return dst
 dst.mkdir(parents=True)
 for d in ['constant','system','100']:shutil.copytree(src/d,dst/d)
 for n in ['mesh_reference.npz','case_manifest.json']:shutil.copy2(src/n,dst/n)
 p=dst/'constant/aerosolProperties';p.write_text(re.sub(r'heatEachW\s+[^;]+;', 'heatEachW 0;',p.read_text()))
 p=dst/'system/controlDict';text=re.sub(r'endTime\s+[^;]+;', 'endTime 140;',p.read_text());text=re.sub(r'writeInterval\s+[^;]+;', 'writeInterval 5;',text);p.write_text(text)
 (dst/'intervention.json').write_text(json.dumps({'start_s':100,'end_s':140,'change':'heat to air off from100; emissions until120 unchanged'},indent=2))
 return dst
if __name__=='__main__':
 for l in ['corner','center']:print(prepare(l))
 print(prepare_heat_off())
