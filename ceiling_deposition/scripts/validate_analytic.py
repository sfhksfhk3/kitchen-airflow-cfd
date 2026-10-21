"""Meaningful closed-form limits and discrete conservation tests."""
from pathlib import Path
import numpy as np,pandas as pd,json
from scipy.special import exprel
ROOT=Path(__file__).resolve().parents[1]
def conductance(v,R):
 x=np.asarray(v)*R
 # exprel(-x) = [exp(-x)-1]/(-x), hence v/[1-exp(-x)]
 with np.errstate(over='ignore'):return 1/(R*exprel(-x))
R=1e5
assert abs(conductance(0,R)-1/R)<1e-15
assert abs(conductance(1e-3,R)-1e-3)<1e-15
assert conductance(-1e-3,R)<1e-40
vs=np.array([-1e-5,0,1e-5]);assert np.all(np.diff(conductance(vs,R))>0)
# Implicit Euler well-mixed box: exactly conserve input + remaining + removal.
m=0.;dep=0.;emit=0.;dt=.1;lam=.03
for i in range(1000):
 q=1e-9 if i<500 else 0;m=(m+dt*q)/(1+dt*lam);dep+=dt*lam*m;emit+=dt*q
assert abs(m+dep-emit)<1e-20
# Identical surfaces must not be assigned an unsupported numerical material ranking.
cfg=json.loads((ROOT/'config/long_term.json').read_text());assert all(x['sticking_probabilities']==[.1,.5,1.] for x in cfg['materials'])
result={'zero_drift_limit':'PASS','advective_limit':'PASS','opposing_gravity_suppression':'PASS','monotonic_normal_drift':'PASS','implicit_box_mass_conservation_error_kg':m+dep-emit,'materials_no_invented_alpha_ranking':'PASS','units':{'Y':'dimensionless','phi':'kg/s','J_dep':'kg/(m2 s)','M_dep_native':'kg/m2','M_dep_CSV':'mg/m2','normalized_M':'(mg/m2)/(mg emitted)'}}
(ROOT/'verification/analytic_validation.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
