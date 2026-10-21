from pathlib import Path
import json,numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.colors import LogNorm
from matplotlib.backends.backend_pdf import PdfPages
from foam_mesh import mesh_geometry
from process_ceiling import ROOT,OUT
from long_term import cycle_units
FIG=ROOT/'figures';FIG.mkdir(exist_ok=True)
from matplotlib import font_manager
font_path=Path('/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc')
if font_path.exists():font_manager.fontManager.addfont(str(font_path))
CJK_FONT=font_manager.FontProperties(fname=str(font_path)).get_name() if font_path.exists() else 'sans-serif'
plt.rcParams.update({'font.family':CJK_FONT,'axes.unicode_minus':False,'font.size':10,'savefig.dpi':220})
DOORS=['both_open','both_closed','left_open','right_open'];ZH={'both_open':'兩門全開','both_closed':'兩門全關','left_open':'左開右關','right_open':'右開左關'}
def label(case):return ('靠左' if case.startswith('corner') else '置中')+'・'+next((ZH[d] for d in DOORS if d in case),case)
def main_cases():return [f'{l}_{d}_parallel' for d in DOORS for l in ['corner','center']]
def save(fig,name):
 engine=fig.get_layout_engine()
 if engine and hasattr(engine,'set'):engine.set(rect=(0,.055,1,.945))
 fig.text(.5,.006,'未通過天花沉積網格敏感性；僅供模型檢查，不作實際積油量或材料壽命判定',ha='center',fontsize=7,color='#9c3030')
 for ext in ['png','pdf']:fig.savefig(FIG/f'{name}.{ext}',bbox_inches=None)
 plt.close(fig)
def heatmap(ax,case,norm):
 fp=OUT/case/('ceiling_complete_event_faces.csv' if (OUT/case/'ceiling_complete_event_faces.csv').exists() else 'ceiling_faces.csv');f=pd.read_csv(fp);g=mesh_geometry(ROOT/'runs'/case);ids=f.original_face_id.to_numpy(dtype=int);polys=[g['points'][g['faces'][i],:2] for i in ids];val=f.M_unit_mg_m2_per_mg.to_numpy()
 summary=json.loads((OUT/case/('complete_event_summary.json' if (OUT/case/'complete_event_summary.json').exists() else 'summary.json')).read_text())
 note=f"參考總量 {summary['M_total_reference_mg']:.2e} mg"
 if 'both_closed' not in case:note+='\n數值近零，不定位實際污漬'
 ax.text(.015,.98,note,transform=ax.transAxes,va='top',fontsize=8,bbox={'facecolor':'white','alpha':.85,'edgecolor':'none'})
 pc=PolyCollection(polys,array=np.maximum(val,norm.vmin),cmap='inferno',norm=norm,edgecolors='none');ax.add_collection(pc)
 ax.set(xlim=(0,3.3),ylim=(0,1.65),xlabel='沿檯面 x（m）',ylabel='距後牆 y（m）',title=label(case));ax.set_aspect('equal')
 ax.set_facecolor('#d9dde2');x0=.7 if case.startswith('corner') else 1.55
 ax.plot([x0,x0+.9],[.3,.3],color='cyan',lw=1.8,label='煙機前緣')
 xburn=[.95,1.35] if case.startswith('corner') else [1.8,2.2]
 ax.scatter(xburn,[.3,.3],c='lime',s=28,marker='x')
 j=np.argmax(val);ax.scatter([f.iloc[j].x_m],[f.iloc[j].y_m],s=70,marker='*',c='deepskyblue',edgecolors='white')
 return pc

def make_plots():
 existing=[c for c in main_cases() if (OUT/c/'summary.json').exists()];summ={c:json.loads((OUT/c/'summary.json').read_text()) for c in existing}
 for c in existing:
  if (OUT/c/'complete_event_summary.json').exists():summ[c].update(json.loads((OUT/c/'complete_event_summary.json').read_text()))
 for door in DOORS:
  pair=[f'{l}_{door}_parallel' for l in ['corner','center']]
  if not all(c in existing for c in pair):continue
  maxv=max(summ[c]['M_max_unit_mg_m2_per_mg'] for c in pair);maxv=max(maxv,1e-30);norm=LogNorm(maxv*1e-8,maxv)
  fig,axes=plt.subplots(2,1,figsize=(9,9),layout='constrained')
  for ax,c in zip(axes,pair):pc=heatmap(ax,c,norm)
  fig.colorbar(pc,ax=axes,shrink=.8,label='累積沉積 mg/m²／每排放 1 mg 油滴（對數色階）')
  fig.suptitle(f'天花積油俯視・{ZH[door]}\n120 s 排放＋240 s 瞬態＋末端流場清空積分；同色階',fontsize=15)
  fig.text(.04,.012,'灰色：高櫃／煙機佔用區；綠×：爐頭投影；藍★：最大網格面。色階下限不是檢出限。',fontsize=8)
  save(fig,'01_ceiling_'+door+'_1to1')
 if existing:
  fig,ax=plt.subplots(figsize=(9,9),layout='constrained');labels=[label(c) for c in existing];v=[summ[c]['M_total_unit_mg_per_mg'] for c in existing]
  ax.barh(labels,np.maximum(v,1e-300),color=['#3b82a0' if c.startswith('corner') else '#d77447' for c in existing]);ax.set_xscale('log');ax.set_xlabel('天花總沉積 mg／每排放 1 mg 油滴');ax.set_title('八工況天花總沉積（已形成油滴；非氣相冷凝）')
  save(fig,'02_total_deposition_1to1')
  fig,ax=plt.subplots(figsize=(9,9),layout='constrained');keys=['hood_fraction','left_fraction','right_fraction','other_wall_fraction','ceiling_fraction','airborne_fraction'];names=['抽油煙機','左門淨流出','右門淨流出','其他表面','天花','末端懸浮'];left=np.zeros(len(existing))
  for k,n in zip(keys,names):
   vv=np.array([summ[c][k] for c in existing]);ax.barh(labels,vv,left=left,label=n);left+=vv
  ax.set_xlim(0,1.02);ax.set_xlabel('總輸入油滴質量比例');ax.set_title('質量去向（很小的天花份額可能在堆疊圖中不可見）');ax.legend(loc='lower center',bbox_to_anchor=(.5,1.02),ncol=3)
  save(fig,'03_mass_balance_1to1')
 closed=[c for c in existing if 'both_closed' in c]
 if closed:
  fig,axes=plt.subplots(2,1,figsize=(9,9),layout='constrained')
  for ax,events in zip(axes,[1,2]):
   for c in closed:
    m=summ[c]['M_mean_unit_mg_m2_per_mg'];days=np.arange(1,366)
    ax.plot(days,m*events*days,label=label(c)+'・不清潔')
    ax.plot(days,[m*cycle_units(int(d),events,30,.9)[0] for d in days],ls='--',label=label(c)+'・30日清潔90%')
   ax.set(title=f'每日 {events} 次；每次排放 1 mg 的條件式累積',xlabel='日數',ylabel='天花平均 mg/m²');ax.legend(fontsize=9)
  fig.suptitle('長期積油情境：沒有預測黃變、老化或厚度',fontsize=14);save(fig,'04_long_term_normalized_1to1')
 if existing:
  fig,ax=plt.subplots(figsize=(9,9),layout='constrained')
  for c in closed:
   f=pd.read_csv(OUT/c/'ceiling_complete_event_faces.csv') if (OUT/c/'ceiling_complete_event_faces.csv').exists() else pd.read_csv(OUT/c/'ceiling_faces.csv');vals=[float(f[f'M_{d:g}um_unit_mg_m2_per_mg']@f.area_m2) for d in [.1,.3,1,3,10]];ax.loglog([.1,.3,1,3,10],np.maximum(vals,1e-300),marker='o',label=label(c))
  ax.set(xlabel='獨立油滴粒徑（µm）',ylabel='天花沉積／該粒徑輸入質量',title='五粒徑沉積響應（不是粒子數分布）');ax.legend();ax.grid(alpha=.2);save(fig,'05_size_response_1to1')
 modelrows=[]
 for c in existing:
  other=c.removesuffix('_parallel')
  if (OUT/other/'summary.json').exists():modelrows.append((c,json.loads((OUT/other/'summary.json').read_text())))
 if modelrows:
  fig,ax=plt.subplots(figsize=(9,9),layout='constrained');idx=np.arange(len(modelrows))
  ax.barh(idx-.18,[max(r['M_total_unit_mg_per_mg'],1e-300) for c,r in modelrows],height=.35,label='漂移＋等效阻力閉合')
  ax.barh(idx+.18,[max(json.loads((OUT/c/'summary.json').read_text())['M_total_unit_mg_per_mg'],1e-300) for c,r in modelrows],height=.35,label='慣性接觸另計（主篩查情境）')
  ax.set_yticks(idx,[label(c) for c,r in modelrows]);ax.set_xscale('log');ax.set_xlabel('天花沉積 mg／每排放 1 mg');ax.set_title('0–240 s 壁面閉合差異：非實驗驗證或嚴格上下界');ax.legend();save(fig,'06_wall_model_uncertainty_1to1')
 with PdfPages(FIG/'all_ceiling_figures_1to1.pdf') as pdf:
  import matplotlib.image as mpimg
  for p in sorted(FIG.glob('*_1to1.png')):
   fig=plt.figure(figsize=(9,9));ax=fig.add_axes([0,0,1,1]);ax.imshow(mpimg.imread(p));ax.axis('off');pdf.savefig(fig);plt.close(fig)
 return existing

if __name__=='__main__':print('Plotted',make_plots())
