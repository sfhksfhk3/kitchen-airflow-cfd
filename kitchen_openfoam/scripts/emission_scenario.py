"""Literature-informed oil-MASS scenario; not a measured number distribution."""
from pathlib import Path
import json,numpy as np
from scipy.special import ndtr
ROOT=Path(__file__).resolve().parents[1]
DIAMETERS=np.array([.1,.3,1,3,10.])
def weights(median_um=2.5,gsd=1.5):
    edges=np.r_[0,np.sqrt(DIAMETERS[:-1]*DIAMETERS[1:]),np.inf]
    with np.errstate(divide='ignore'):cdf=ndtr(np.log(edges/median_um)/np.log(gsd))
    return np.diff(cdf)

def config():
    return {'label_zh':'中式炒菜油滴總質量參考情境（非住家實測）','diameters_um':DIAMETERS.tolist(),'mass_weights':weights().tolist(),'mass_median_um':2.5,'mass_geometric_standard_deviation_assumed':1.5,'rate_per_burner_mg_min':7.2,'total_rate_mg_min':14.4,'total_rate_ug_s':240.,'interpretation':'7.2 mg/min/頭取自商用中式廚房CFD文獻的0.12 mg/s油滴源。2.5µm質量中位徑參考六道中式菜餚研究；GSD1.5及離散化權重為本次工程假設。五個單粒徑響應先分別求解，再線性重組。這不是完整超細粒子數分布，也不是PM2.5濃度。','rate_sensitivity_factors':[.25,1,4],'width_sensitivity_gsd':[1.2,1.5,2.]}

if __name__=='__main__':
    p=ROOT/'emission_scenario.json';p.write_text(json.dumps(config(),ensure_ascii=False,indent=2));print(p)
