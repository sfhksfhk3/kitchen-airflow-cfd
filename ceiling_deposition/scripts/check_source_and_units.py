from pathlib import Path
import sys,json,hashlib,re,numpy as np,pandas as pd
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT.parent/'kitchen_openfoam/scripts'))
from foam_io import read_internal
rows=[]
for p in (ROOT/'runs').glob('*/RUN_COMPLETE.json'):
 c=p.parent;m=json.loads((c/'replay_manifest.json').read_text());source=Path(m['source']);unchanged=True
 for name,sha in m['source_mesh_sha256'].items():
  assert hashlib.sha256((source/'constant/polyMesh'/name).read_bytes()).hexdigest()==sha,(c,name)
 shape=read_internal(c/f"{m['start_s']:.12g}"/'sourceShape',len(np.load(c/'mesh_reference.npz')['volumes']));vol=np.load(c/'mesh_reference.npz')['volumes'];integral=float(shape@vol);assert abs(integral-1)<1e-8,(c,integral)
 rows.append({'case':c.name,'source_shape_integral':integral,'source_mesh_original_hashes_unchanged':unchanged})
pd.DataFrame(rows).to_csv(ROOT/'verification/source_integrals_and_preservation.csv',index=False)
flow=[]
for c in (ROOT/'flow_sensitivity').iterdir():
 csvs=list((c/'diagnostics').glob('aerosol_*.csv'))
 if not csvs:continue
 d=pd.concat(pd.read_csv(p) for p in csvs);d=d[(d.diameter_um==1)&(d.time_s>=110)]
 flow.append({'case':c.name,'Q_late_m3_s':float(d[d.time_s>=120].Qhood_m3_s.mean()),'reconstructed_energy_residual_RMS_W':float(np.sqrt(np.mean(d.energy_balance_W**2))),'air_mass_residual_RMS_kg_s':float(np.sqrt(np.mean(d.air_mass_residual_kg_s**2))),'max_T_K':float(d.Tmax_K.max()),'pressure_cap_observed_max_Pa':float(d.hood_pressure_drop_Pa.max())})
pd.DataFrame(flow).to_csv(ROOT/'verification/coupled_interventions_quality.csv',index=False)
print('Source profiles/unchanged original meshes verified:',len(rows));print(pd.DataFrame(flow).to_string(index=False))
