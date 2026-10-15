from pathlib import Path
import json,numpy as np
from visualize import TetraField,slice_data,draw_solids,plt
root=Path(__file__).resolve().parent
cases=['both_open','both_open_refined']
metrics={c:json.loads((root/'results'/(c+'_metrics.json')).read_text()) for c in cases}
fields={c:TetraField(dict(np.load(root/'results'/(c+'.npz')))) for c in cases}
slices={c:slice_data(fields[c])['burner'] for c in cases}
fig,axs=plt.subplots(1,2,figsize=(13,4.7),constrained_layout=True)
for ax,c in zip(axs,cases):
 s=slices[c];m=metrics[c]
 im=ax.pcolormesh(s['u'],s['v'],s['temperature'],shading='auto',cmap='inferno',vmin=25,vmax=max(metrics[x]['max_nodal_temperature_C'] for x in cases))
 draw_solids(ax);ax.set_title(f"h = {m['mesh_target_m']:.2f} m · {m['tetrahedra']:,} 四面體\n最高節點溫度 {m['max_nodal_temperature_C']:.2f}°C",fontsize=11)
fig.colorbar(im,ax=axs.tolist(),shrink=.75,label='空氣溫度／°C')
fig.suptitle('網格加密比較｜兩門全開\n最高溫度仍對網格敏感；不能宣稱已達網格無關',fontsize=14)
(root/'figures').mkdir(exist_ok=True)
fig.savefig(root/'figures/mesh_comparison.png',dpi=180,bbox_inches='tight');plt.close(fig)
base,ref=[metrics[c] for c in cases]
summary={'base_h_m':base['mesh_target_m'],'refined_h_m':ref['mesh_target_m'],'occupied_zone_speed_relative_change_percent':(ref['occupied_zone_mean_speed_m_s']/base['occupied_zone_mean_speed_m_s']-1)*100,'max_temperature_change_C':ref['max_nodal_temperature_C']-base['max_nodal_temperature_C'],'max_temperature_rise_relative_change_percent':((ref['max_nodal_temperature_C']-25)/(base['max_nodal_temperature_C']-25)-1)*100,'refined_divergence_gradient_ratio':ref['divergence_gradient_norm_ratio'],'grid_independence_demonstrated':False}
(root/'results/mesh_comparison.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2))
print(json.dumps(summary,ensure_ascii=False,indent=2))
