"""Integrate the air-facing high-cabinet side at x=0.7 m.

Temperature is the Robin interface value in an AIR-ONLY model, not a
validated conjugate solid temperature. Direct flame/pot radiation is absent.
"""
from pathlib import Path
import json
import numpy as np
from skfem import MeshTet
import meshio

ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'deliverables'
OUT.mkdir(exist_ok=True)

def analyze(path, name):
    d = np.load(path)
    m = MeshTet(d['points'].T, d['cells'].T)
    fs = m.boundary_facets()
    faces = m.facets[:, fs].T
    vertices = d['points'][faces]
    centers = vertices.mean(axis=1)
    choose = np.all(np.abs(vertices[:, :, 0] - 0.7) < 1e-8, axis=1)
    choose &= (centers[:, 1] >= -1e-8) & (centers[:, 1] <= 0.6 + 1e-8)
    faces = faces[choose]
    vertices = vertices[choose]
    centers = centers[choose]
    area = np.linalg.norm(np.cross(vertices[:, 1] - vertices[:, 0],
                                   vertices[:, 2] - vertices[:, 0]), axis=1) / 2
    temp = d['temperature'][faces]
    mean = temp.mean(axis=1)
    heat_flux = 3.0 * (mean - 25.0)
    power = area * heat_flux
    nodes, local = np.unique(faces, return_inverse=True)
    local = local.reshape(-1, 3)
    maximum = nodes[np.argmax(d['temperature'][nodes])]
    result = dict(name=name, source_npz=str(path),
                  wall='air-facing cabinet side x=0.7 m, y=0–0.6 m',
                  area_m2=float(area.sum()), total_net_heat_into_wall_W=float(power.sum()),
                  area_average_temperature_C=float(np.dot(area, mean) / area.sum()),
                  minimum_interface_temperature_C=float(temp.min()),
                  maximum_interface_temperature_C=float(temp.max()),
                  maximum_temperature_position_m=d['points'][maximum].tolist(),
                  area_average_heat_flux_W_m2=float(power.sum() / area.sum()),
                  maximum_nodal_heat_flux_W_m2=float(3 * (temp.max() - 25)),
                  robin_sink_temperature_C=25.0, robin_conductance_W_m2K=3.0,
                  temperature_scope='Robin interface temperature; no solid conduction or direct radiation solved',
                  triangular_facets=len(faces), boundary_nodes=len(nodes))
    (OUT / f'{name}_wall_metrics.json').write_text(json.dumps(result, ensure_ascii=False, indent=2))
    np.savetxt(OUT / f'{name}_wall_faces.csv',
               np.c_[faces, centers, area, mean, heat_flux, power], delimiter=',',
               header='node0,node1,node2,x_m,y_m,z_m,area_m2,mean_interface_temperature_C,net_heat_flux_W_m2,net_heat_into_wall_W', comments='')
    np.savetxt(OUT / f'{name}_wall_nodes.csv',
               np.c_[nodes, d['points'][nodes], d['temperature'][nodes], 3*(d['temperature'][nodes]-25)],
               delimiter=',', header='global_node_id,x_m,y_m,z_m,interface_temperature_C,net_robin_heat_flux_W_m2', comments='')
    meshio.write(OUT / f'{name}_wall.vtu', meshio.Mesh(d['points'][nodes], [('triangle', local)],
                 point_data={'interface_temperature_C':d['temperature'][nodes],
                             'net_robin_heat_flux_W_m2':3*(d['temperature'][nodes]-25)},
                 cell_data={'area_m2':[area], 'net_heat_into_wall_W':[power]}))
    return result, d['points'][nodes], local, d['temperature'][nodes]

def main():
    jobs = [('4000W_coarse', ROOT.parent/'corner/results/both_closed.npz'),
            ('4000W_refined', ROOT.parent/'corner/results/both_closed_refined.npz'),
            ('9000W_coarse', ROOT/'results/both_closed_9000W_coarse.npz'),
            ('9000W_refined', ROOT/'results/both_closed_9000W_refined.npz'),
            ('9000W_fine', ROOT/'results/both_closed_9000W_fine.npz')]
    outputs = {name: analyze(path, name) for name, path in jobs}
    metrics = {name: data[0] for name, data in outputs.items()}
    fine, coarse = metrics['9000W_fine'], metrics['9000W_refined']
    validation = {key:dict(coarse=coarse[key], refined=fine[key],
                          relative_change_percent=100*(fine[key]-coarse[key])/abs(fine[key]))
                  for key in ['total_net_heat_into_wall_W', 'area_average_temperature_C',
                              'maximum_interface_temperature_C']}
    validation['temperature_rise_average_relative_change_percent'] = 100*(fine['area_average_temperature_C']-coarse['area_average_temperature_C'])/(fine['area_average_temperature_C']-25)
    metrics['grid_comparison_h014_to_h010'] = validation
    metrics['primary'] = '9000W_fine'
    (OUT/'wall_summary.json').write_text(json.dumps(metrics, ensure_ascii=False, indent=2))
    import csv
    with (OUT/'wall_summary.csv').open('w') as f:
        fields=['name','area_m2','total_net_heat_into_wall_W','area_average_temperature_C',
                'maximum_interface_temperature_C','area_average_heat_flux_W_m2','maximum_nodal_heat_flux_W_m2']
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
        for name in outputs:w.writerow({k:metrics[name][k] for k in fields})
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.tri as mtri
    plt.rcParams['font.family']='Noto Sans CJK JP'
    p, tri, t = outputs['9000W_fine'][1:]
    fig, ax = plt.subplots(figsize=(7,9))
    tr=mtri.Triangulation(p[:,1],p[:,2],tri)
    plot=ax.tripcolor(tr,t,shading='gouraud',cmap='inferno')
    ax.set_aspect('equal');ax.set_xlabel('距後牆 y（m）');ax.set_ylabel('高度 z（m）')
    ax.set_title('靠左雙爐・兩門關閉・總放熱 9000 W\n高櫃側面：簡化 Robin 邊界溫度')
    fig.colorbar(plot,ax=ax,label='模型邊界溫度（°C）')
    pos=fine['maximum_temperature_position_m']
    ax.plot(pos[1],pos[2],marker='x',color='cyan',markersize=9)
    fig.text(.12,.025,f"受熱面積 {fine['area_m2']:.4f} m²；總淨受熱 {fine['total_net_heat_into_wall_W']:.1f} W\n平均 {fine['area_average_temperature_C']:.1f}°C；最高 {fine['maximum_interface_temperature_C']:.1f}°C（青色 ×）\n未計固體導熱及火焰／鍋具直接輻射",fontsize=10)
    fig.tight_layout(rect=[0,.10,1,1]);fig.savefig(OUT/'wall_temperature_9000W.png',dpi=170);plt.close(fig)
    import plotly.graph_objects as go
    views=[]
    for name in ['4000W_refined','9000W_fine']:
        r, p, tri, t = outputs[name]
        views.append(go.Mesh3d(x=p[:,0],y=p[:,1],z=p[:,2],i=tri[:,0],j=tri[:,1],k=tri[:,2],
                     intensity=t,intensitymode='vertex',colorscale='Inferno',
                     cmin=25,cmax=fine['maximum_interface_temperature_C'],
                     colorbar=dict(title='邊界 °C'),visible=name=='9000W_fine',
                     customdata=np.c_[3*(t-25)],
                     hovertemplate='x=%{x:.3f} m<br>y=%{y:.3f} m<br>z=%{z:.3f} m<br>邊界溫度 %{intensity:.1f}°C<br>淨熱通量 %{customdata[0]:.1f} W/m²<extra></extra>', name=name))
    figure=go.Figure(views)
    figure.update_layout(title='高櫃右側面受熱：模型邊界溫度（不含固體導熱／直接輻射）',
                         scene=dict(xaxis_title='x (m)',yaxis_title='y (m)',zaxis_title='z (m)',aspectmode='data',camera=dict(eye=dict(x=1.8,y=1.2,z=.7))),
                         updatemenus=[dict(buttons=[dict(label='原工況 4000 W',method='update',args=[{'visible':[True,False]}]),
                                                  dict(label='新工況 9000 W',method='update',args=[{'visible':[False,True]}])])],
                         margin=dict(l=0,r=0,b=0,t=85))
    figure.write_html(OUT/'wall_interactive_3d.html',include_plotlyjs=True,full_html=True)
    print(json.dumps(metrics,ensure_ascii=False,indent=2))

if __name__=='__main__':main()
