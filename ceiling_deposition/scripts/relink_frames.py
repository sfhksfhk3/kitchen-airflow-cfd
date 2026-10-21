"""Relink archived read-only frame symlinks after moving the project."""
from pathlib import Path
import os,argparse
ROOT=Path(__file__).resolve().parents[1]
if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('--flow-root',type=Path,default=Path('/workspace/kitchen_openfoam'));args=a.parse_args();count=0
 for p in ROOT.rglob('*'):
  if not p.is_symlink():continue
  old=Path(os.readlink(p));new=None
  for prefix,target in [(Path('/workspace/kitchen-airflow-cfd/ceiling_deposition'),ROOT),(Path('/workspace/kitchen_openfoam'),args.flow_root.resolve())]:
   if old.is_absolute() and old.is_relative_to(prefix):new=target/old.relative_to(prefix);break
  if new is not None and new!=old:p.unlink();p.symlink_to(new);count+=1
 print('Updated read-only frame links:',count)
