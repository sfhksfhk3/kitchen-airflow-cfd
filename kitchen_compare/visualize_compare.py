"""Offline scientific visualization of computed FEM tracer/layout comparisons."""
from pathlib import Path
import json,html,csv,math
import numpy as np
from scipy.spatial import cKDTree
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.patches import Rectangle
from plotly.offline import get_plotlyjs
from plotly.colors import get_colorscale
ROOT=Path(__file__).resolve().parent
CASES=['both_open','both_closed','left_open','right_open'];LAYOUTS=['corner','center']
LABEL={'both_open':'兩門全開','both_closed':'兩門全關','left_open':'左開右關','right_open':'右開左關','corner':'靠左角落','center':'檯面中央'}
plt.rcParams.update({'font.family':'Noto Sans CJK JP','axes.unicode_minus':False,'font.size':10})
FACES=np.array([[0,2,1],[0,3,2],[4,5,6],[4,6,7],[0,1,5],[0,5,4],[1,2,6],[1,6,5],[2,3,7],[2,7,6],[3,0,4],[3,4,7]])
EDGES=[(0,1),(1,2),(2,3),(3,0),(4,5),(5,6),(6,7),(7,4),(0,4),(1,5),(2,6),(3,7)]
def vertices(b):
 x0,x1,y0,y1,z0,z1=b
 return np.array([[x0,y0,z0],[x1,y0,z0],[x1,y1,z0],[x0,y1,z0],[x0,y0,z1],[x1,y0,z1],[x1,y1,z1],[x0,y1,z1]])
def boxes(c):
 h=c['hood'];x0=c['counter']['x0'];L=c['room']['length'];H=c['room']['height'];ud=c['upper_cabinet_depth_m']
 out=[('高櫃',[0,x0,0,.6,0,H],'#aabcc4'),('檯面／下櫃',[x0,L,0,.6,0,.9],'#c7d4da'),('煙機',[h['x0'],h['x1'],0,.3,1.6,H],'#294e65')]
 if h['x0']>x0+1e-5:out.append(('左上櫃',[x0,h['x0'],0,ud,1.6,H],'#b5c8d0'))
 if h['x1']<L-1e-5:out.append(('右上櫃',[h['x1'],L,0,ud,1.6,H],'#b5c8d0'))
 return out
def baffle(c):
 h=c['hood'];b=c['baffle'];a=np.pi/4;half=b['thickness_m_assumed']/2
 y0,z0=.3,1.6;y1,z1=y0+b['length_m']*np.cos(a),z0-b['length_m']*np.sin(a)
 yz=[(y0+half*np.sin(a),z0+half*np.cos(a)),(y0-half*np.sin(a),z0-half*np.cos(a)),(y1-half*np.sin(a),z1-half*np.cos(a)),(y1+half*np.sin(a),z1+half*np.cos(a))]
 return np.array([[x,y,z] for x in [h['x0'],h['x1']] for y,z in yz])
class Field:
 def __init__(self,d):
  self.p=d['points'];self.c=d['cells'];self.a=self.p[self.c[:,0]];self.inv=np.linalg.inv(np.stack([self.p[self.c[:,k]]-self.a for k in [1,2,3]],axis=-1));self.tree=cKDTree(self.p[self.c].mean(axis=1))
  self.val=np.c_[d['velocity'],d['temperature'],d['pressure'],d['concentration_steady'],d['concentration_frames'].T]
  self.rt=d['velocity_RT0_cells'];self.times=d['frame_times_s']
 def at(self,q):
  q=np.asarray(q).reshape(-1,3);out=np.full((len(q),self.val.shape[1]),np.nan);rt=np.full((len(q),3),np.nan);remaining=np.arange(len(q))
  for count in [24,128]:
   if not len(remaining):break
   _,ix=self.tree.query(q[remaining],k=min(count,len(self.c)));b=np.einsum('nkij,nkj->nki',self.inv[ix],q[remaining,None,:]-self.a[ix]);w=np.concatenate([1-b.sum(2,keepdims=True),b],axis=2)
   ok=(w>=-1e-7).all(2)&(w<=1+1e-7).all(2);found=ok.any(1);rr=np.flatnonzero(found);jj=np.argmax(ok[found],axis=1);ci=ix[rr,jj]
   out[remaining[found]]=np.einsum('ni,nij->nj',w[rr,jj],self.val[self.c[ci]]);rt[remaining[found]]=self.rt[ci];remaining=remaining[~found]
  return out,rt

def serial(x):
 if isinstance(x,np.ndarray):return serial(x.tolist())
 if isinstance(x,dict):return {k:serial(v) for k,v in x.items()}
 if isinstance(x,list):return [serial(v) for v in x]
 if isinstance(x,(float,np.floating)):return float(f'{float(x):.7g}') if np.isfinite(x) else None
 if isinstance(x,np.integer):return int(x)
 return x

def model(c,case):
 traces=[]
 for name,box,color in boxes(c):
  p=vertices(box);traces.append(dict(type='mesh3d',x=p[:,0],y=p[:,1],z=p[:,2],i=FACES[:,0],j=FACES[:,1],k=FACES[:,2],color=color,opacity=.32,name=name,flatshading=True,hoverinfo='name'))
  e=np.array([p[i] if i is not None else [np.nan]*3 for edge in EDGES for i in [*edge,None]])
  traces.append(dict(type='scatter3d',mode='lines',x=e[:,0],y=e[:,1],z=e[:,2],line=dict(color='#6b8797',width=1.5),hoverinfo='skip',showlegend=False))
 p=baffle(c)
 # Prism uses yz-corner order on each x end; make all six faces explicitly.
 fs=np.array([[0,1,2],[0,2,3],[4,6,5],[4,7,6],[0,4,5],[0,5,1],[1,5,6],[1,6,2],[2,6,7],[2,7,3],[3,7,4],[3,4,0]])
 traces.append(dict(type='mesh3d',x=p[:,0],y=p[:,1],z=p[:,2],i=fs[:,0],j=fs[:,1],k=fs[:,2],color='#e5a037',opacity=.95,name='擋板 10 cm／45°／僅煙機90 cm範圍',hoverinfo='name',flatshading=True))
 p=vertices([0,3.3,0,1.65,0,2.2]);e=np.array([p[i] if i is not None else [np.nan]*3 for edge in EDGES for i in [*edge,None]])
 traces.append(dict(type='scatter3d',mode='lines',x=e[:,0],y=e[:,1],z=e[:,2],line=dict(color='#8099a7',width=2,dash='dot'),name='房間邊界',hoverinfo='skip'))
 for side,x in [('left',0),('right',3.3)]:
  op=case=='both_open' or case==side+'_open';z0=0 if op else .005
  traces.append(dict(type='mesh3d',x=[x]*4,y=[.7,1.6,1.6,.7],z=[z0,z0,2.2,2.2],i=[0,0],j=[1,2],k=[2,3],color='#28ae99' if op else '#8e9da8',opacity=.10 if op else .4,name=('左門' if side=='left' else '右門')+('全開' if op else '關閉／底縫5mm'),hoverinfo='name'))
  if case=='both_closed':traces.append(dict(type='mesh3d',x=[x]*4,y=[.7,1.6,1.6,.7],z=[0,0,.1,.1],i=[0,0],j=[1,2],k=[2,3],color='#ebb73f',opacity=.7,name='等效補風带0.10m（非實際門縫）',showlegend=side=='left',hoverinfo='name'))
 for i,(x,y,z) in enumerate(c['burners']):
  a=np.linspace(0,2*np.pi,45);traces.append(dict(type='scatter3d',mode='lines',x=x+.085*np.cos(a),y=y+.085*np.sin(a),z=np.full(len(a),z+.01),line=dict(color='#dc613e',width=6),name='燃氣爐2000 W／個',showlegend=i==0,hoverinfo='name'))
 return serial(traces)

def slices(field,c):
 out={}
 for key,axis,value in [('horizontal','z',1.5),('front','y',.45),('side','x',c['layout_center_x_m']-.2)]:
  u=np.linspace(0,3.3 if axis!='x' else 1.65,65);v=np.linspace(0,1.65 if axis=='z' else 2.2,39);UU,VV=np.meshgrid(u,v)
  q=np.c_[UU.ravel(),VV.ravel(),np.full(UU.size,value)] if axis=='z' else (np.c_[UU.ravel(),np.full(UU.size,value),VV.ravel()] if axis=='y' else np.c_[np.full(UU.size,value),UU.ravel(),VV.ravel()])
  values,rt=field.at(q);shape=UU.shape
  out[key]=dict(axis=axis,value=value,u=u,v=v,temperature=values[:,3].reshape(shape),pressure=values[:,4].reshape(shape),speed=np.linalg.norm(rt,axis=1).reshape(shape),concentration=[values[:,6+j].reshape(shape) for j in range(len(field.times))])
 return out

def write_obj(c):
 lines=['# metre units; corrected three-view geometry including 0.9m-width baffle'];idx=0
 for name,box,_ in boxes(c):
  p=vertices(box);lines.append('o cabinet_'+str(idx));lines.extend('v '+' '.join(map(str,q)) for q in p);lines.extend('f '+' '.join(str(idx+int(v)+1) for v in f) for f in FACES);idx+=8
 p=baffle(c);lines.append('o hood_baffle_90cm');lines.extend('v '+' '.join(map(str,q)) for q in p)
 fs=[[0,1,2],[0,2,3],[4,6,5],[4,7,6],[0,4,5],[0,5,1],[1,5,6],[1,6,2],[2,6,7],[2,7,3],[3,7,4],[3,4,0]]
 lines.extend('f '+' '.join(str(idx+v+1) for v in f) for f in fs)
 (ROOT/c['layout']/'geometry.obj').write_text('\n'.join(lines)+'\n')

def plots(data,metrics,configs):
 f=ROOT/'figures';f.mkdir(exist_ok=True)
 fig,axes=plt.subplots(4,2,figsize=(12,12),constrained_layout=True)
 for row,case in enumerate(CASES):
  for col,layout in enumerate(LAYOUTS):
   ax=axes[row,col];d=data[layout+'_'+case];s=d['slices']['horizontal'];ci=d['times'].index(120.);h=configs[layout]['hood']
   im=ax.pcolormesh(s['u'],s['v'],np.maximum(s['concentration'][ci],1e-7),cmap='magma',norm=LogNorm(1e-7,10),shading='auto')
   ax.add_patch(Rectangle((0,0),.7,.6,facecolor='#becbd2',edgecolor='#69818d'))
   ax.add_patch(Rectangle((h['x0'],0),.9,.3,fill=False,edgecolor='#54cad1',linestyle='--',lw=1.5))
   for x,y,_ in configs[layout]['burners']:ax.plot(x,y,'o',color='#5bd5d1',markersize=3)
   m=metrics[layout+'_'+case]['tracer'];ax.set_title(f"{LABEL[layout]} · {LABEL[case]} · 滯留 {m['retained_fraction_at_120s']*100:.2f}%",fontsize=11)
   ax.set(xlim=(0,3.3),ylim=(0,1.65),xlabel='x／m',ylabel='y／m');ax.set_aspect('equal')
 fig.suptitle('油煙示蹤物空間分布｜釋放120秒\n呼吸高度水平剖面 z=1.50m · 同源強、同色階（對數）· 非實測油煙濃度',fontsize=15)
 fig.colorbar(im,ax=axes.ravel().tolist(),shrink=.75,label='歸一化濃度 c/R（s/m³）')
 fig.savefig(f/'smoke_eight_cases.png',dpi=180,bbox_inches='tight');plt.close(fig)
 fig,axes=plt.subplots(1,2,figsize=(12,5),constrained_layout=True);xx=np.arange(4)
 for j,(key,title,factor) in enumerate([('steady_equivalent_residence_time_s','全室平均滯留時間 τ（s）',1),('occupied_zone_normalized_concentration_at_120s','走道濃度指數（c/R ×1000）',1000)]):
  for i,l in enumerate(LAYOUTS):
   vals=[metrics[l+'_'+case]['tracer'][key]*factor for case in CASES];axes[j].bar(xx+(i-.5)*.36,vals,.36,label=LABEL[l],color=['#2589a1','#df9450'][i])
  axes[j].set(xticks=xx,xticklabels=[LABEL[c] for c in CASES],ylabel=title,title=title);axes[j].legend();axes[j].grid(axis='y',alpha=.15)
 fig.suptitle('兩種位置的油煙滯留比較｜雙爐合計4000 W\n全室滯留與走道濃度是不同指標',fontsize=15)
 fig.savefig(f/'retention_comparison.png',dpi=180,bbox_inches='tight');plt.close(fig)
 fig,axes=plt.subplots(1,2,figsize=(12,4.8),constrained_layout=True)
 colors=['#308bb1','#bf633b','#62946b','#815bad']
 for i,l in enumerate(LAYOUTS):
  for case,color in zip(CASES,colors):
   a=np.array(data[l+'_'+case]['history']);at120=a[np.argmin(abs(a[:,0]-120)),1];sel=(a[:,0]>=120)&(a[:,0]<=180)
   axes[i].semilogy(a[sel,0]-120,np.maximum(a[sel,1]/at120,1e-8),color=color,label=LABEL[case])
  axes[i].axhline(.05,color='gray',ls='--',label='剩餘5%');axes[i].set(title=LABEL[l],xlabel='停止示蹤物排放後時間／s',ylabel='室內存量／停止時存量',ylim=(1e-5,1));axes[i].legend(fontsize=8);axes[i].grid(alpha=.2)
 fig.suptitle('清除曲線｜只停止示蹤物，風機及4000 W熱負荷保持不變',fontsize=14)
 fig.savefig(f/'clearance_curves.png',dpi=180,bbox_inches='tight');plt.close(fig)
 fig,axes=plt.subplots(1,2,figsize=(12,4.8),constrained_layout=True)
 for i,l in enumerate(LAYOUTS):
  for case,color in zip(CASES,colors):
   a=np.array(data[l+'_'+case]['history']);at120=a[np.argmin(abs(a[:,0]-120)),7];sel=(a[:,0]>=120)&(a[:,0]<=420)
   axes[i].semilogy(a[sel,0]-120,np.maximum(a[sel,7]/at120,1e-8),color=color,label=LABEL[case])
  axes[i].axhline(.05,color='gray',ls='--',label='剩餘5%');axes[i].set(title=LABEL[l],xlabel='停止示蹤物排放後时间／s',ylabel='走道平均濃度／停止時濃度',ylim=(1e-4,1));axes[i].legend(fontsize=8);axes[i].grid(alpha=.2)
 fig.suptitle('走道尾端清除｜濃度衰減對混合及網格較敏感，非實際關火秒數',fontsize=14)
 fig.savefig(f/'walk_clearance_curves.png',dpi=180,bbox_inches='tight');plt.close(fig)


HTML=r'''<!doctype html><html lang="zh-Hant"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>廚房油煙滯留｜角落與中央</title><style>
*{box-sizing:border-box}body{margin:0;background:#edf3f6;color:#1c394d;font-family:system-ui,"Noto Sans CJK TC",sans-serif}header{padding:30px 4%;background:#16364d;color:white}h1{font-size:29px;margin:10px 0}header p{max-width:1060px;line-height:1.7;font-size:14px;color:#c5d8e2}.eyebrow{font-size:11px;letter-spacing:3px;color:#80d8cb}main{max-width:1520px;margin:auto;padding:22px}.controls{background:white;border:1px solid #d8e3e9;border-radius:10px;padding:17px;display:flex;gap:14px;flex-wrap:wrap;align-items:end}label{display:block;font-size:11px;color:#658091;margin-bottom:6px}select{padding:9px 12px;background:#f5f8fa;border:1px solid #cad9e3;border-radius:5px;color:#224354;font-size:13px}#time{min-width:215px}.cards{display:grid;grid-template-columns:repeat(4,1fr);gap:14px;margin:17px 0}.card,.panel{background:white;border:1px solid #d8e3e9;border-radius:10px}.card{padding:16px}.card small{color:#6a8392}.number{font-size:27px;font-weight:650;margin:6px 0}.unit{font-size:12px;font-weight:400}.detail{font-size:11px;color:#7a8e9b}.plots{display:grid;grid-template-columns:1fr 1fr;gap:17px}.panel h2{font-size:15px;margin:17px 18px 7px}.panel p{font-size:11px;color:#718895;margin:5px 18px 12px;line-height:1.6}.plot{height:530px}.wide{margin-top:17px}.tablewrap{padding:12px 18px 18px;overflow:auto}table{border-collapse:collapse;width:100%;font-size:12px}th,td{padding:10px;border-bottom:1px solid #e3ebef;text-align:right;white-space:nowrap}td:first-child,th:first-child{text-align:left}th{background:#f4f8fa;color:#668090;font-weight:500}.active{background:#edf8f6}.notes{padding:15px 21px;font-size:13px;line-height:1.8}.notes ul{padding-left:18px}.notes strong{color:#9a5521}.foot{font-size:11px;color:#718b9a;padding:16px 0}.badge{color:#b4eed8;border:1px solid #70aa9c;border-radius:20px;padding:5px 10px;font-size:12px;margin-left:12px}@media(max-width:950px){.plots{grid-template-columns:1fr}.cards{grid-template-columns:1fr 1fr}.plot{height:480px}}@media(max-width:600px){main{padding:12px}header{padding:22px}h1{font-size:24px}.cards{gap:8px}.card{padding:11px}.number{font-size:22px}}
</style><script>__PLOTLY__</script></head><body><header><div class="eyebrow">KITCHEN · FINITE ELEMENT / TRACER STUDY</div><h1>油煙滯留｜靠左角落 vs 檯面中央 <span class="badge">條件式數值比較</span></h1><p>依新三視圖建模：頂吸煙機90×30 cm、距灶面70 cm，前緣僅90 cm範圍設10 cm／45°擋板；双燃氣爐各向空氣放熱2000 W，合計4000 W。中央方案整體中心x=2.0 m。油煙以等量被動示蹤物計算，濃度按排放率歸一化。</p></header><main><div class="controls"><div><label>煙機與爐頭位置</label><select id="layout"><option value="corner">靠左角落 · 中心1.15 m</option><option value="center">檯面中央 · 中心2.00 m</option></select></div><div><label>門狀態</label><select id="case"><option value="both_open">兩門全開</option><option value="both_closed">兩門全關 · 各5 mm底縫</option><option value="left_open">左開右關</option><option value="right_open">右開左關</option></select></div><div><label>顯示物理量</label><select id="quantity"><option value="concentration">油煙示蹤物濃度</option><option value="speed">守恆輸運風速</option><option value="temperature">空氣溫度</option><option value="pressure">FEM相對壓力</option></select></div><div><label>剖面</label><select id="plane"><option value="horizontal">呼吸高度 · z=1.50 m</option><option value="front">爐前垂直 · y=0.45 m</option><option value="side">左爐側面剖面</option></select></div><div><label>示蹤物時間 · 0–120s排放，之後停止</label><select id="time"></select></div></div><div class="cards"><div class="card"><small>全室平均滯留時間</small><div class="number" id="tau"></div><div class="detail">穩態室內存量 ÷ 相同源強</div></div><div class="card"><small>120s後仍留在室內</small><div class="number" id="retained"></div><div class="detail">佔120s累積釋放量</div></div><div class="card"><small>走道相對濃度指數</small><div class="number" id="exposure"></div><div class="detail">c/R ×1000 · 不是實際mg/m³</div></div><div class="card"><small>走道清除95%示蹤物</small><div class="number" id="clear"></div><div class="detail">熱負荷及風機保持不變的測試</div></div></div><div class="plots"><section class="panel"><h2>01 · 3D 幾何與計算濃度</h2><p>拖曳旋轉、滾輪縮放；點雲是原始FEM節點濃度，不是合成煙霧；橙色為45°擋板。</p><div id="scene" class="plot"></div></section><section class="panel"><h2 id="sliceTitle">02 · 油煙濃度剖面</h2><p id="sliceNote"></p><div id="slice" class="plot"></div></section></div><section class="panel wide"><h2>03 · 室內存量與走道尾端滯留</h2><p>停止示蹤物排放後，保持4000 W熱負荷及風機；此測試不代表熄火後的真實清除秒數。</p><div id="history" style="height:300px"></div></section><section class="panel wide"><h2>04 · 八工況比較</h2><div class="tablewrap"><table><thead><tr><th>情境</th><th>排風 m³/min</th><th>平均滯留 s</th><th>120s滯留 %</th><th>走道 c/R×1000</th><th>走道清除95% s</th></tr></thead><tbody id="table"></tbody></table></div></section><section class="panel wide notes"><h2>如何解讀</h2><p style="font-size:13px;color:#546f80">__CONCLUSION__</p><ul><li>實際油煙排放率未知：若總排放率為R，真實濃度需將本頁歸一化濃度乘以R；沒有提供健康限值判斷。</li><li>穩態下所有示蹤物最終可能都流向煙機，這不等於實機捕煙率100%。請比較室內存量、走道濃度及有限時間清除。</li><li>每門底縫5 mm，全關流量依假設風機曲線估算；室內補氣带鋪展至0.10 m，未解析門縫噴流。</li><li>無室外風壓，未解析瞬時湍流、油滴沉降或牆面沉積；氣流仍是固定基準密度近似，局部高溫區的熱膨脹未完整建模。</li><li>輸運速度經RT0逐單元守恆重建；原始P1 FEM速度另存於VTU。網格／擴散敏感性與重建幅度見中文報告。</li></ul></section><div class="foot">離線檔案 · 原始FEM + RT0守恆輸運 + 被動示蹤物 · 完整數據及假設見 report.pdf</div></main><script>
const DATA=__DATA__;const els=['layout','case','quantity','plane','time'].map(x=>document.getElementById(x));const pal=DATA.colormaps;const cam={eye:{x:1.5,y:1.8,z:1.2},up:{x:0,y:0,z:1}};for(let i=0;i<DATA.times.length;i++){let t=DATA.times[i],o=document.createElement('option');o.value=i;o.textContent=t+' s'+(t>120?' · 停止後'+(t-120)+' s':' · 排放中');o.selected=t===120;els[4].append(o)}
const base={paper_bgcolor:'white',plot_bgcolor:'white',font:{family:'system-ui,Noto Sans CJK TC,sans-serif',size:11,color:'#36556a'},margin:{l:55,r:65,t:25,b:45}};const desc={concentration:['油煙示蹤物濃度','c/R · s/m³'],speed:['守恆輸運風速','m/s'],temperature:['空氣溫度','°C'],pressure:['FEM相對壓力','Pa']};
function range(q){return DATA.ranges[q]}function scale(q){return pal[q]}function colorVals(q,v){return q==='concentration'?v.map(row=>Array.isArray(row)?row.map(x=>x===null?null:Math.log10(Math.max(x,1e-7))):row===null?null:Math.log10(Math.max(row,1e-7))):v}
function render(){let [l,c,q,p,fi]=els.map(x=>x.value);fi=Number(fi);let key=l+'_'+c,d=DATA.cases[key],m=d.tracer,s=d.slices[p],vals=q==='concentration'?s.concentration[fi]:s[q],limits=range(q),color=colorVals(q,vals),unit=desc[q][1];document.getElementById('tau').innerHTML=m.steady_equivalent_residence_time_s.toFixed(2)+' <span class="unit">s</span>';document.getElementById('retained').innerHTML=(100*m.retained_fraction_at_120s).toFixed(2)+' <span class="unit">%</span>';document.getElementById('exposure').innerHTML=(m.occupied_zone_normalized_concentration_at_120s*1000).toFixed(3);document.getElementById('clear').innerHTML=m.clearance_occupied_zone_t95_s+' <span class="unit">s</span>';document.getElementById('sliceTitle').textContent='02 · '+desc[q][0]+'剖面';document.getElementById('sliceNote').textContent=DATA.labels[l]+' · '+DATA.labels[c]+' · '+s.axis+'='+s.value.toFixed(2)+' m · '+DATA.times[fi]+' s';
 let cloud=d.cloud;let cv=q==='concentration'?cloud.concentration[fi]:cloud[q],cc=colorVals(q,cv);let active=cv.map((v,i)=>v!==null&&(q!=='concentration'||v>1e-7)?i:-1).filter(i=>i>=0),picker=v=>active.map(i=>v[i]);let trace={type:'scatter3d',mode:'markers',x:picker(cloud.x),y:picker(cloud.y),z:picker(cloud.z),marker:{size:q==='concentration'?3.2:2.5,color:picker(cc),colorscale:scale(q),cmin:limits[0],cmax:limits[1],opacity:.65},customdata:picker(cv),name:desc[q][0],hovertemplate:'x=%{x:.2f} m<br>y=%{y:.2f} m<br>z=%{z:.2f} m<br>'+desc[q][0]+'=%{customdata:.5g} '+unit+'<extra></extra>'};let scene={...base,margin:{l:0,r:0,t:0,b:0},uirevision:'kitchen-camera',scene:{xaxis:{title:'x / m',range:[0,3.3]},yaxis:{title:'y / m',range:[0,1.65]},zaxis:{title:'z / m',range:[0,2.2]},aspectmode:'data',camera:cam},legend:{orientation:'h',x:0,y:-.04,font:{size:9}}};Plotly.react('scene',d.model.concat([trace]),scene,{responsive:true,displaylogo:false});let cb={title:{text:unit},thickness:13,len:.8};if(q==='concentration'){cb.tickvals=[-7,-5,-3,-1,0,1];cb.ticktext=['1e-7','1e-5','1e-3','0.1','1','10']}let heat={type:'heatmap',x:s.u,y:s.v,z:color,colorscale:scale(q),zmin:limits[0],zmax:limits[1],customdata:vals,colorbar:cb,hoverongaps:false,connectgaps:false,hovertemplate:'%{customdata:.5g} '+unit+'<extra></extra>'};let shapes=[];if(p==='horizontal')shapes=[{type:'rect',x0:0,x1:.7,y0:0,y1:.6,line:{color:'#8195a0'},fillcolor:'#bbcbd3'}];if(p==='front')shapes=DATA.boxes[l].filter(b=>b[2]<=s.value&&b[3]>=s.value).map(b=>({type:'rect',x0:b[0],x1:b[1],y0:b[4],y1:b[5],line:{color:'#8195a0'},fillcolor:'#bbcbd3'}));let xl=p==='side'?'y / m':'x / m',yl=p==='horizontal'?'y / m':'z / m';Plotly.react('slice',[heat],{...base,xaxis:{title:xl,range:[0,p==='side'?1.65:3.3],constrain:'domain'},yaxis:{title:yl,range:[0,p==='horizontal'?1.65:2.2],scaleanchor:'x',scaleratio:1,constrain:'domain'},shapes},{responsive:true,displaylogo:false});
 let hh=d.history;Plotly.react('history',[{type:'scatter',mode:'lines',x:hh.map(a=>a[0]),y:hh.map(a=>a[1]),line:{color:'#2286a3',width:2},name:'室內存量（源強等效秒）'},{type:'scatter',mode:'lines',x:hh.map(a=>a[0]),y:hh.map(a=>a[7]/hh.find(a=>a[0]===120)[7]),yaxis:'y2',line:{color:'#d38a48',width:2,dash:'dot'},name:'走道濃度／停止時濃度'}],{...base,margin:{l:60,r:30,t:15,b:45},xaxis:{title:'時間／s',range:[0,300]},yaxis:{title:'室內示蹤物存量／R（s）'},yaxis2:{title:'走道濃度／120s時濃度',overlaying:'y',side:'right',range:[0,1.1]},legend:{orientation:'h',y:1.13},shapes:[{type:'line',x0:120,x1:120,y0:0,y1:1,yref:'paper',line:{dash:'dot',color:'#c18750'}}]},{responsive:true,displaylogo:false});document.getElementById('table').innerHTML=Object.entries(DATA.cases).map(([k,a])=>{let v=a.tracer;return '<tr class="'+(k===key?'active':'')+'"><td>'+DATA.labels[k.split('_')[0]]+'／'+DATA.labels[v.case]+'</td><td>'+a.flow.fan_flow_m3_min.toFixed(2)+'</td><td>'+v.steady_equivalent_residence_time_s.toFixed(2)+'</td><td>'+(v.retained_fraction_at_120s*100).toFixed(2)+'</td><td>'+(v.occupied_zone_normalized_concentration_at_120s*1000).toFixed(3)+'</td><td>'+v.clearance_occupied_zone_t95_s+'</td></tr>'}).join('');window.COMPARE_READY=true}
els.forEach(e=>e.addEventListener('change',render));render();
</script></body></html>'''

def main():
 configs={l:json.loads((ROOT/l/'config.json').read_text()) for l in LAYOUTS};data={};metrics={};probe_rows=[]
 for l in LAYOUTS:
  write_obj(configs[l])
  for case in CASES:
   name=case+'_refined';r=ROOT/l/'results';flow=json.loads((r/(name+'_metrics.json')).read_text());smoke=json.loads((r/(name+'_tracer_metrics.json')).read_text())
   if not flow['converged']:raise ValueError('Wait for converged FEM')
   field=Field(dict(np.load(r/(name+'_tracer.npz'))));s=slices(field,configs[l]);history=np.loadtxt(r/(name+'_tracer_history.csv'),delimiter=',',skiprows=1)[::4]
   # Original nodal concentrations, consistently sampled for the point cloud.
   ii=np.linspace(0,len(field.p)-1,min(1300,len(field.p))).astype(int);rv,_=field.at(field.p[ii]);rt,_=field.at(field.p[ii])
   vv=field.val[ii];cloud=dict(x=field.p[ii,0],y=field.p[ii,1],z=field.p[ii,2],temperature=vv[:,3],pressure=vv[:,4],concentration=[vv[:,6+j] for j in range(len(field.times))])
   _,UR=field.at(field.p[ii]);cloud['speed']=np.linalg.norm(UR,axis=1)
   key=l+'_'+case;metrics[key]=dict(flow=flow,tracer=smoke);data[key]=dict(flow=flow,tracer=smoke,model=model(configs[l],case),slices=s,cloud=cloud,history=history,times=field.times.tolist())
   cc=configs[l]['layout_center_x_m'];points=[('操作者呼吸點示例',cc,.85,1.5),('爐前空氣',cc,.55,1.25),('房間走道中點',1.65,1.15,1.5)]
   vals,_=field.at([p[1:] for p in points]);fi=list(field.times).index(120.)
   for p,row in zip(points,vals):probe_rows.append([l,case,*p,row[5],row[6+fi]])
   print('visualize',l,case,flush=True)
 plots(data,metrics,configs)
 with (ROOT/'probes.csv').open('w',newline='',encoding='utf-8-sig') as f:
  w=csv.writer(f);w.writerow(['layout','door_state','probe','x_m','y_m','z_m','normalized_c_steady_s_m3','normalized_c_120s_s_m3']);w.writerows(probe_rows)
 times=data['corner_both_open']['times'];assert all(d['times']==times for d in data.values())
 payload=dict(cases=data,times=times,labels=LABEL,boxes={l:[b for _,b,_ in boxes(c)] for l,c in configs.items()},ranges={'concentration':[-7,1],'speed':[0,5],'temperature':[25,max(m['flow']['max_nodal_temperature_C'] for m in metrics.values())],'pressure':[min(m['flow']['FE_pressure_min_Pa'] for m in metrics.values()),0]},colormaps={k:get_colorscale(n) for k,n in [('concentration','Magma'),('speed','Turbo'),('temperature','Inferno'),('pressure','RdBu')]})
 concise='本次無室外風壓的條件式模型，初步傾向保留靠左位置並讓兩門補風：中央對總滯留量的優勢較小且未穩定，但走道濃度普遍較高。單開左／右門的尾端清除對網格與混合假設敏感，精確秒數尚未確認；全關工況補風不足。'
 page=HTML.replace('__PLOTLY__',get_plotlyjs()).replace('__DATA__',json.dumps(serial(payload),ensure_ascii=False,separators=(',',':')).replace('</','<\\/')).replace('__CONCLUSION__',html.escape(concise))
 (ROOT/'interactive_comparison.html').write_text(page)
 (ROOT/'comparison_metrics.json').write_text(json.dumps(metrics,ensure_ascii=False,indent=2))
 print('offline visualization ready',flush=True)
if __name__=='__main__':main()
