#!/usr/bin/env python3
"""Create offline views directly from solved tetrahedral P1 FEM nodal data.
No synthetic smoke; streamlines integrate interpolated computed velocity.
"""
from pathlib import Path
import argparse,csv,json,html,math
import numpy as np
from scipy.spatial import cKDTree
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import plotly.graph_objects as go
from plotly.offline import get_plotlyjs
from plotly.colors import get_colorscale

ROOT=Path(__file__).resolve().parent
CASES=['both_open','both_closed','left_open','right_open']
LABELS={'both_open':'兩門全開','both_closed':'兩門全關','left_open':'左開右關','right_open':'右開左關'}
BOXES=[('左高櫃',[0,.6,0,.6,0,2.2],'#a5b4bc'),('下櫃',[.6,3.3,0,.6,0,.9],'#bbc9ce'),('頂吸煙機',[.6,1.5,0,.6,1.6,2.2],'#354f62'),('上櫃',[1.5,3.3,0,.6,1.6,2.2],'#b1c1c7')]
plt.rcParams.update({'font.family':'Noto Sans CJK JP','axes.unicode_minus':False,'figure.facecolor':'#f4f7f9','axes.facecolor':'#ffffff','axes.edgecolor':'#cbd5df','text.color':'#183447','axes.labelcolor':'#183447','xtick.color':'#526779','ytick.color':'#526779','font.size':10})

class TetraField:
    """Find original tetrahedra, not a convex hull interpolant across solids."""
    def __init__(self,d):
        self.p=d['points'];self.c=d['cells'];self.v=d['velocity'];self.t=d['temperature'];self.pressure=d['pressure']
        self.a=self.p[self.c[:,0]]
        self.inv=np.linalg.inv(np.stack([self.p[self.c[:,k]]-self.a for k in (1,2,3)],axis=-1))
        self.tree=cKDTree(self.p[self.c].mean(axis=1))
    def at(self,q):
        q=np.asarray(q).reshape(-1,3)
        out=np.full((len(q),5),np.nan)
        alive=np.flatnonzero(np.isfinite(q).all(axis=1))
        for k in [16,64]:
            if not len(alive):break
            _,ix=self.tree.query(q[alive],k=min(k,len(self.c)))
            b=np.einsum('nkij,nkj->nki',self.inv[ix],q[alive,None,:]-self.a[ix])
            weights=np.concatenate([1-b.sum(axis=2,keepdims=True),b],axis=2)
            valid=(weights>=-2e-7).all(axis=2)&(weights<=1+2e-7).all(axis=2)
            ok=valid.any(axis=1);rr=np.flatnonzero(ok);ci=np.argmax(valid[ok],axis=1)
            ids=ix[rr,ci];w=weights[rr,ci];nodes=self.c[ids]
            out[alive[ok],:3]=np.einsum('ni,nij->nj',w,self.v[nodes])
            out[alive[ok],3]=np.einsum('ni,ni->n',w,self.t[nodes])
            out[alive[ok],4]=np.einsum('ni,ni->n',w,self.pressure[nodes])
            alive=alive[~ok]
        return out

def serial(a):
    if isinstance(a,np.ndarray):return serial(a.tolist())
    if isinstance(a,list):return [serial(x) for x in a]
    if isinstance(a,dict):return {k:serial(v) for k,v in a.items()}
    if isinstance(a,(np.floating,float)):return float(a) if np.isfinite(a) else None
    if isinstance(a,np.integer):return int(a)
    return a

def box_vertices(b):
    x0,x1,y0,y1,z0,z1=b
    return np.array([[x0,y0,z0],[x1,y0,z0],[x1,y1,z0],[x0,y1,z0],[x0,y0,z1],[x1,y0,z1],[x1,y1,z1],[x0,y1,z1]])
FACES=np.array([[0,2,1],[0,3,2],[4,5,6],[4,6,7],[0,1,5],[0,5,4],[1,2,6],[1,6,5],[2,3,7],[2,7,6],[3,0,4],[3,4,7]])
EDGES=[(0,1),(1,2),(2,3),(3,0),(4,5),(5,6),(6,7),(7,4),(0,4),(1,5),(2,6),(3,7)]

def model_traces(case):
    traces=[]
    for name,b,color in BOXES:
        p=box_vertices(b)
        traces.append(dict(type='mesh3d',x=p[:,0].tolist(),y=p[:,1].tolist(),z=p[:,2].tolist(),i=FACES[:,0].tolist(),j=FACES[:,1].tolist(),k=FACES[:,2].tolist(),color=color,opacity=.50,name=name,showlegend=True,hoverinfo='name',flatshading=True))
        edges=np.array([p[j] if j is not None else [np.nan]*3 for e in EDGES for j in [e[0],e[1],None]])
        traces.append(dict(type='scatter3d',mode='lines',x=serial(edges[:,0]),y=serial(edges[:,1]),z=serial(edges[:,2]),line=dict(color='#54707f',width=2),showlegend=False,hoverinfo='skip'))
    # Kitchen outline, not an invented enclosed shell.
    p=box_vertices([0,3.3,0,1.65,0,2.2]);edges=np.array([p[j] if j is not None else [np.nan]*3 for e in EDGES for j in [e[0],e[1],None]])
    traces.append(dict(type='scatter3d',mode='lines',x=serial(edges[:,0]),y=serial(edges[:,1]),z=serial(edges[:,2]),line=dict(color='#869ca9',width=2,dash='dot'),name='房間邊界',hoverinfo='name'))
    traces.append(dict(type='mesh3d',x=[0,3.3,3.3,0],y=[0,0,1.65,1.65],z=[0]*4,i=[0,0],j=[1,2],k=[2,3],color='#dae3e8',opacity=.35,name='地面',showlegend=False,hoverinfo='skip'))
    # Actual 5 mm gaps shown geometrically; coarse inlet bands are distinct.
    for side,x in [('left',0),('right',3.3)]:
        opened=case=='both_open' or case==side+'_open'
        z0=0 if opened else .005
        traces.append(dict(type='mesh3d',x=[x]*4,y=[.7,1.6,1.6,.7],z=[z0,z0,2.2,2.2],i=[0,0],j=[1,2],k=[2,3],color='#40c7ae' if opened else '#87929c',opacity=.11 if opened else .48,name=('左' if side=='left' else '右')+('門：全開' if opened else '門：關閉，底縫 5 mm'),showlegend=True,hoverinfo='name'))
        if not opened and case=='both_closed':
            traces.append(dict(type='mesh3d',x=[x]*4,y=[.7,1.6,1.6,.7],z=[0,0,.1,.1],i=[0,0],j=[1,2],k=[2,3],color='#e4aa3a',opacity=.75,name='等效補氣帶 0.10 m（非真實門縫高度）',showlegend=side=='left',hoverinfo='name'))
    for i,x in enumerate([.85,1.25]):
        theta=np.linspace(0,2*np.pi,50)
        traces.append(dict(type='scatter3d',mode='lines',x=(x+.095*np.cos(theta)).tolist(),y=(.3+.095*np.sin(theta)).tolist(),z=[.905]*50,line=dict(color='#e46c38',width=7),name='爐頭：假設向空氣放熱 600 W／個',showlegend=i==0,hoverinfo='name'))
    traces.append(dict(type='mesh3d',x=[.6,1.5,1.5,.6],y=[0,0,.6,.6],z=[1.599]*4,i=[0,0],j=[1,2],k=[2,3],color='#e1b24c',opacity=.75,name='煙機吸入口 0.9 × 0.6 m',showlegend=True,hoverinfo='name'))
    return traces

def streamline_data(field,case):
    seeds=[]
    # Integrate steady velocity direction with a 1 cm spatial midpoint step.
    for x in [.85,1.25]:
        for dx,dy in [(-.06,-.04),(0,0),(.06,.04),(0,.09)]:seeds.append([x+dx,.3+dy,.96])
    for side,x in [('left',.018),('right',3.282)]:
        opened=case=='both_open' or case==side+'_open'
        if opened:
            for y in [.85,1.1,1.4]:
                for z in [.35,.8,1.3,1.8]:seeds.append([x,y,z])
        elif case=='both_closed':
            for y in [.85,1.1,1.4]:seeds.append([x,y,.05])
    pos=np.array(seeds,float);paths=[[p.copy()] for p in pos];speeds=[[] for _ in pos];active=np.ones(len(pos),bool)
    for _ in range(650):
        ix=np.flatnonzero(active)
        if not len(ix):break
        r=field.at(pos[ix]);v=r[:,:3];s=np.linalg.norm(v,axis=1);ok=np.isfinite(s)&(s>2e-4)
        active[ix[~ok]]=False;ix=ix[ok];v=v[ok];s=s[ok]
        if not len(ix):break
        mid=pos[ix]+.005*v/s[:,None];r2=field.at(mid);v2=r2[:,:3];s2=np.linalg.norm(v2,axis=1);ok=np.isfinite(s2)&(s2>2e-4)
        active[ix[~ok]]=False;ix=ix[ok];v2=v2[ok];s2=s2[ok]
        new=pos[ix]+.01*v2/s2[:,None]
        valid=np.isfinite(field.at(new)).all(axis=1)
        for j,p,sp,keep in zip(ix,new,s2,valid):
            if keep:paths[j].append(p.copy());speeds[j].append(float(sp));pos[j]=p
            else:active[j]=False
        # Terminate stagnating/recurrent paths; no claim of full transport.
        for j in ix:
            if len(paths[j])>120 and np.linalg.norm(pos[j]-paths[j][-100])<.06:active[j]=False
    x=[];y=[];z=[];sp=[]
    for p,s in zip(paths,speeds):
        if len(p)<4:continue
        q=np.asarray(p);speed=s+[s[-1]]
        x.extend(q[:,0].tolist()+[None]);y.extend(q[:,1].tolist()+[None]);z.extend(q[:,2].tolist()+[None]);sp.extend(speed+[0])
    return dict(type='scatter3d',mode='lines',x=x,y=y,z=z,line=dict(color=sp,colorscale=get_colorscale('Turbo'),cmin=0,cmax=1.4,width=4),name='速度場積分流線',hoverinfo='skip',showlegend=True),paths

def slice_data(field):
    slices={}
    for key,axis,value,aa,bb,na,nb in [('burner','y',.3,(0,3.3),(0,2.2),111,77),('aisle','y',1.05,(0,3.3),(0,2.2),111,77),('horizontal','z',1.2,(0,3.3),(0,1.65),111,63)]:
        u=np.linspace(*aa,na);v=np.linspace(*bb,nb);U,V=np.meshgrid(u,v)
        if axis=='y':q=np.column_stack([U.ravel(),np.full(U.size,value),V.ravel()])
        else:q=np.column_stack([U.ravel(),V.ravel(),np.full(U.size,value)])
        r=field.at(q);speed=np.linalg.norm(r[:,:3],axis=1)
        # Mask exact solid boundary too: no visual extrapolation through cabinets.
        for _,b,_ in BOXES:
            inside=(q[:,0]>=b[0]-1e-8)&(q[:,0]<=b[1]+1e-8)&(q[:,1]>=b[2]-1e-8)&(q[:,1]<=b[3]+1e-8)&(q[:,2]>=b[4]-1e-8)&(q[:,2]<=b[5]+1e-8)
            speed[inside]=np.nan;r[inside]=np.nan
        arrays=dict(speed=speed.reshape(U.shape),temperature=r[:,3].reshape(U.shape),pressure=r[:,4].reshape(U.shape))
        stride=6;Q=q.reshape(*U.shape,3)[::stride,::stride].reshape(-1,3);R=r.reshape(*U.shape,5)[::stride,::stride].reshape(-1,5)
        uu=R[:,0];vv=R[:,2] if axis=='y' else R[:,1]
        slices[key]=dict(axis=axis,value=value,u=u,v=v,**arrays,arrows=dict(x=Q[:,0],y=Q[:,2] if axis=='y' else Q[:,1],u=uu,v=vv))
    return slices

def write_obj(out):
    lines=['# Kitchen geometry in metres; x width, y depth, z elevation.','# Solid simplifications, not measured as-built geometry.']
    n=0
    for name,b,_ in BOXES:
        lines.append('o '+{'左高櫃':'tall_cabinet','下櫃':'base_cabinet','頂吸煙機':'hood_envelope','上櫃':'upper_cabinet'}[name])
        for p in box_vertices(b):lines.append('v '+' '.join(f'{v:.6f}' for v in p))
        for f in FACES:lines.append('f '+' '.join(str(n+int(i)+1) for i in f))
        n+=8
    (out/'geometry.obj').write_text('\n'.join(lines)+'\n')

def draw_solids(ax,plane='burner'):
    if plane=='burner':
        for _,b,color in BOXES:ax.add_patch(Rectangle((b[0],b[4]),b[1]-b[0],b[5]-b[4],facecolor=color,edgecolor='#6e8592',lw=.7,zorder=4))
        ax.plot([.6,1.5],[1.6,1.6],color='#c99727',lw=3,zorder=5)
        for x in [.85,1.25]:ax.plot(x,.915,'o',color='#dc693b',markersize=5,zorder=5)
    ax.set(xlim=(0,3.3),ylim=(0,2.2),xlabel='x／m',ylabel='z／m');ax.set_aspect('equal')

def static_plots(out,data,metrics):
    fig,axes=plt.subplots(2,2,figsize=(13.5,9.2),constrained_layout=True)
    for ax,case in zip(axes.ravel(),CASES):
        s=data[case]['slices']['burner'];im=ax.pcolormesh(s['u'],s['v'],s['speed'],cmap='turbo',vmin=0,vmax=1.4,shading='auto',rasterized=True)
        a=s['arrows'];speed=np.hypot(a['u'],a['v']);nz=np.isfinite(speed)&(speed>1e-3)
        ax.quiver(a['x'][nz],a['y'][nz],a['u'][nz],a['v'][nz],color='#152b3b',scale=10,width=.0026,alpha=.8)
        draw_solids(ax)
        m=metrics[case];ax.set_title(f"{LABELS[case]}　Q = {m['fan_flow_m3_min']:.2f} m³/min",fontsize=12,loc='left',pad=10)
    fig.suptitle('廚房氣流｜四種開門情境\n爐頭中心剖面 y = 0.30 m · 原始 FEM 場插值 · 共用色階',fontsize=17,fontweight='bold')
    cb=fig.colorbar(im,ax=axes.ravel().tolist(),shrink=.8,pad=.02,extend='max');cb.set_label('速度大小／m/s')
    fig.savefig(out/'four_cases_velocity.png',dpi=180,bbox_inches='tight');plt.close(fig)
    fig,axes=plt.subplots(2,2,figsize=(13.5,9.2),constrained_layout=True)
    tmax=max(float(np.nanmax(data[c]['slices']['burner']['temperature'])) for c in CASES)
    for ax,case in zip(axes.ravel(),CASES):
        s=data[case]['slices']['burner'];im=ax.pcolormesh(s['u'],s['v'],s['temperature'],cmap='inferno',vmin=25,vmax=tmax,shading='auto',rasterized=True);draw_solids(ax);ax.set_title(LABELS[case],loc='left',fontsize=12)
    fig.suptitle('空氣溫度估計｜雙爐各向空氣放熱 600 W\n爐頭中心剖面 y = 0.30 m · 非火焰／鍋具溫度',fontsize=17,fontweight='bold')
    cb=fig.colorbar(im,ax=axes.ravel().tolist(),shrink=.8,pad=.02);cb.set_label('空氣溫度／°C')
    fig.savefig(out/'four_cases_temperature.png',dpi=180,bbox_inches='tight');plt.close(fig)
    fig=plt.figure(figsize=(12,8));ax=fig.add_subplot(111,projection='3d')
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    for name,b,color in BOXES:
        v=box_vertices(b);ax.add_collection3d(Poly3DCollection(v[FACES],facecolor=color,edgecolor='#76909e',linewidths=.4,alpha=.23))
    for path in data['both_open']['paths']:
        p=np.asarray(path)
        if len(p)>3:ax.plot(*p.T,color='#258cba',alpha=.85,lw=1.4)
    for x in [.85,1.25]:ax.scatter(x,.3,.915,c='#e67b38',s=45)
    ax.set(xlim=(0,3.3),ylim=(0,1.65),zlim=(0,2.2),xlabel='x／m',ylabel='y／m',zlabel='z／m',title='3D 幾何與計算流線｜兩門全開\n流線由穩態速度場積分，非油煙顆粒軌跡')
    ax.set_box_aspect((3.3,1.65,2.2));ax.view_init(25,-60);fig.tight_layout();fig.savefig(out/'model_3d.png',dpi=180,bbox_inches='tight');plt.close(fig)

def write_probes(out,fields):
    probes=[('左爐上方20cm',.85,.3,1.1),('右爐上方20cm',1.25,.3,1.1),('煙機入口下方5cm',1.05,.3,1.55),('操作者呼吸區示例',1.05,.95,1.5),('走道中央',1.65,1.125,1.1),('右側走道',2.8,1.125,1.1)]
    with (out/'probes.csv').open('w',newline='',encoding='utf-8-sig') as f:
        w=csv.writer(f);w.writerow(['case','label','probe','x_m','y_m','z_m','u_m_s','v_m_s','w_m_s','speed_m_s','air_temperature_C','FE_relative_pressure_Pa'])
        for case,field in fields.items():
            r=field.at([p[1:] for p in probes])
            for p,row in zip(probes,r):w.writerow([case,LABELS[case],*p,*row[:3],np.linalg.norm(row[:3]),*row[3:]])

HTML=r'''<!doctype html><html lang="zh-Hant"><head><meta charset="UTF-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>廚房氣流 · FEM 情境研究</title><style>
:root{font-family:system-ui,"Noto Sans CJK TC",sans-serif;color:#143244;background:#eef3f6}*{box-sizing:border-box}body{margin:0}header{background:#153449;color:white;padding:32px max(24px,calc((100vw - 1440px)/2))}header .eyebrow{color:#81d4cf;font-size:12px;letter-spacing:3px}h1{margin:12px 0;font-size:30px;font-weight:650}header p{max-width:900px;color:#c6d8e2;line-height:1.65;font-size:14px;margin:0}main{max-width:1490px;margin:0 auto;padding:24px}.badge{display:inline-block;font-size:12px;padding:5px 10px;border:1px solid #739fa8;border-radius:30px;color:#b7f0dc;margin-left:12px;vertical-align:middle}.toolbar{display:flex;gap:18px;align-items:end;flex-wrap:wrap;background:white;padding:16px 20px;border-radius:12px;border:1px solid #dbe5eb}.control label{display:block;font-size:11px;text-transform:uppercase;letter-spacing:1px;color:#637a89;margin:0 0 7px}select{font-size:14px;background:#f5f8fa;color:#18394e;border:1px solid #cbd8df;border-radius:6px;padding:9px 35px 9px 11px;min-width:170px}.metrics{display:grid;grid-template-columns:repeat(4,1fr);gap:14px;margin:18px 0}.metric{background:white;border:1px solid #dbe5eb;border-radius:10px;padding:16px 20px}.metric small{font-size:12px;color:#637a89}.value{font-size:27px;margin:7px 0 2px;font-weight:650}.value span{font-size:13px;color:#637a89;font-weight:400}.sub{font-size:11px;color:#748998}.grid{display:grid;grid-template-columns:1.13fr 1fr;gap:18px}.panel{background:white;border:1px solid #dbe5eb;border-radius:12px;overflow:hidden}.panel h2{margin:18px 20px 0;font-size:16px}.panel .note{margin:8px 20px;color:#6e8390;font-size:12px;line-height:1.6}.plot{height:540px}.footnote{padding:0 20px 18px;color:#748998;font-size:11px;line-height:1.6}.panel.wide{margin-top:18px;padding-bottom:20px}.tablewrap{overflow:auto;padding:15px 20px}table{border-collapse:collapse;width:100%;font-size:13px}th,td{text-align:right;padding:12px;border-bottom:1px solid #e5ecf0;white-space:nowrap}th:first-child,td:first-child{text-align:left}th{color:#65808f;font-weight:500;font-size:11px;background:#f5f8fa}.active-row{background:#edf8f6}.notes{display:grid;grid-template-columns:1fr 1fr;gap:25px;padding:12px 22px;font-size:13px;line-height:1.8;color:#536e7e}.notes ul{padding-left:19px;margin:0}.notes strong{color:#244c60}footer{padding:20px 0;color:#78909d;font-size:11px}.link{color:#1a797c}button{border:1px solid #b5c9d3;background:#f7fafb;color:#31586b;border-radius:6px;padding:9px 12px;cursor:pointer}.spacer{flex:1}#error{display:none;background:#ffe8dc;padding:14px;border-radius:8px;margin:15px 0} @media(max-width:950px){.grid{grid-template-columns:1fr}.metrics{grid-template-columns:1fr 1fr}.plot{height:490px}.notes{grid-template-columns:1fr}h1{font-size:24px}}@media(max-width:550px){main{padding:12px}.metrics{gap:8px}.metric{padding:12px}.value{font-size:23px}.toolbar{gap:12px}.control{flex:1}select{width:100%;min-width:120px}.plot{height:420px}header{padding:24px}.badge{margin-left:0}.notes{font-size:12px}}
</style><script>__PLOTLY__</script></head><body><header><div class="eyebrow">KITCHEN AIRFLOW / FINITE ELEMENT STUDY</div><h1>廚房氣流 · 四種開門情境 <span class="badge">條件式初步估計</span></h1><p>依草圖建立 3.3 × 1.65 × 2.2 m 廚房模型，雙燃氣爐同時使用、頂吸煙機距灶面 70 cm。直接檢視有限元求解的速度、空氣溫度與相對壓力；所有圖表可離線使用。</p></header><main><div class="toolbar"><div class="control"><label for="case">門的狀態</label><select id="case"><option value="both_open">兩門全開</option><option value="both_closed">兩門全關 · 底縫 5 mm</option><option value="left_open">左開右關</option><option value="right_open">右開左關</option></select></div><div class="control"><label for="quantity">顯示物理量</label><select id="quantity"><option value="speed">速度大小 · m/s</option><option value="temperature">空氣溫度 · °C</option><option value="pressure">FE 相對壓力 · Pa</option></select></div><div class="control"><label for="plane">切片位置</label><select id="plane"><option value="burner">爐頭中心 · y = 0.30 m</option><option value="aisle">走道 · y = 1.05 m</option><option value="horizontal">水平面 · z = 1.20 m</option></select></div><div class="spacer"></div><button id="reset">重設 3D 視角</button></div><div id="error"></div><div class="metrics"><div class="metric"><small>估計排風量</small><div class="value" id="flow"></div><div class="sub">假設二次風機曲線與補氣阻力</div></div><div class="metric"><small>補氣阻力網路壓降</small><div class="value" id="drop"></div><div class="sub">單獨計算，非圖中 FE 壓力</div></div><div class="metric"><small>走道區平均風速</small><div class="value" id="speed"></div><div class="sub">區域平均，非峰值</div></div><div class="metric"><small>走道區平均空氣溫度</small><div class="value" id="temp"></div><div class="sub">每爐向空氣放熱 600 W 假設</div></div></div><div class="grid"><section class="panel"><h2>01 · 3D 幾何與計算流線</h2><p class="note">拖曳旋轉、滾輪縮放；點選圖例可隱藏櫃體／門／流線／切片。</p><div id="scene" class="plot"></div><div class="footnote">流線由 FEM 速度場逐段積分，色彩表示速度，並非油煙顆粒軌跡。灰色為櫃體，金色為吸入口；關門底縫 5 mm 與數值補氣帶 0.10 m 分開表示。</div></section><section class="panel"><h2 id="sliceTitle">02 · 定量切片</h2><p class="note" id="sliceNote"></p><div id="slice" class="plot"></div><div class="footnote" id="fieldNote"></div></section></div><section class="panel wide"><h2>03 · 情境比較</h2><p class="note">相同幾何、相同熱負荷、相同數值模型；守恆及收斂狀態取自求解器。共用色階便於比較。</p><div class="tablewrap"><table><thead><tr><th>情境</th><th>排風量<br>m³/min</th><th>補氣網路壓降<br>Pa</th><th>走道平均速度<br>m/s</th><th>走道平均溫度<br>°C</th><th>質量通量相對誤差</th><th>迭代收斂</th></tr></thead><tbody id="comparison"></tbody></table></div></section><section class="panel wide"><h2>04 · 如何解讀</h2><div class="notes"><div><strong>這份模型回答的範圍</strong><ul><li>比較補氣門開關對室內大尺度氣流、熱氣上升及排風量的影響。</li><li>30 m³/min 與 1300 Pa 分別視為自由風量與最大靜壓端點，不同時施加；曲線使用 p = 1300 [1 − (Q/0.5)²]。</li><li>壓力切片為室內 FEM 的相對壓力，基準是數值補氣邊界自然牽引条件；與補氣阻力網路壓降分別報告。</li><li>本模型未求解污染物濃度、燃燒或捕集率，不能據此判定燃气使用安全。</li></ul></div><div><strong>主要假設與資料缺口</strong><ul>__ASSUMPTIONS__</ul></div></div><p class="note">數值方法：P1/P1 四面體有限元；穩態流動與熱浮力耦合，固定有效渦黏性。切片使用原始四面體單元重心座標插值，不跨越櫃體填色。圖中箭頭為該切面的速度分量，顏色為三維速度大小。</p></section><footer>離線交付 · 原始節點與四面體：results/*.npz / *.vtu · 探針：probes.csv · 幾何：geometry.obj · 全部長度以 m 計。</footer></main><script>
const DATA=__DATA__;const camera={eye:{x:1.55,y:1.7,z:1.15},up:{x:0,y:0,z:1}};
const Q={speed:{label:'速度大小',unit:'m/s',colors:DATA.colormaps.speed,note:'速度大小為三維向量長度。切面箭頭顯示平面內的速度方向與大小；固定色階可能截斷局部峰值，懸停仍顯示真實數值。'},temperature:{label:'空氣溫度',unit:'°C',colors:DATA.colormaps.temperature,note:'雙爐各向空氣放熱 600 W，環境 25°C。熱源以有限體積區域分布，溫度不是火焰、鍋具或額定爐火溫度。'},pressure:{label:'FE 相對壓力',unit:'Pa',colors:DATA.colormaps.pressure,note:'此壓力是室內有限元場相對於數值补氣邊界的值，不含完整排煙管、濾網與風機壓力。請勿把此色階與上方補氣阻力網路壓降直接相加。'}};
const layoutBase={font:{family:'system-ui,Noto Sans CJK TC,sans-serif',color:'#36576b',size:11},paper_bgcolor:'white',plot_bgcolor:'white',margin:{l:20,r:20,t:15,b:25}};
const controls=['case','quantity','plane'].map(id=>document.getElementById(id));let initial=true;
function fmt(x,n=3){return Number(x).toFixed(n)}
function arrows(s){let x=[],y=[];const a=s.arrows;for(let i=0;i<a.x.length;i++){if(a.u[i]===null||a.v[i]===null)continue;let u=a.u[i],v=a.v[i],m=Math.hypot(u,v);if(m<.004)continue;let scale=.12/Math.max(m,.4),dx=u*scale,dy=v*scale,xx=a.x[i],yy=a.y[i],endx=xx+dx,endy=yy+dy,theta=Math.atan2(dy,dx),len=Math.min(.028,.3*Math.hypot(dx,dy));x.push(xx,endx,null,endx-len*Math.cos(theta-.45),endx,endx-len*Math.cos(theta+.45),null);y.push(yy,endy,null,endy-len*Math.sin(theta-.45),endy,endy-len*Math.sin(theta+.45),null)}return {type:'scatter',mode:'lines',x,y,line:{color:'rgba(16,39,53,.70)',width:1},hoverinfo:'skip',showlegend:false};}
function planeTrace(s,key,range){let x=[],y=[],z=[],c=[];for(let j=0;j<s.v.length;j++)for(let i=0;i<s.u.length;i++){let v=s[key][j][i];if(v===null)continue;x.push(s.u[i]);y.push(s.axis==='y'?s.value:s.v[j]);z.push(s.axis==='y'?s.v[j]:s.value);c.push(v)}return{type:'scatter3d',mode:'markers',x,y,z,marker:{size:2.8,color:c,colorscale:Q[key].colors,cmin:range[0],cmax:range[1],opacity:.7},name:'定量切片',showlegend:true,visible:'legendonly',hovertemplate:'x=%{x:.2f} m<br>y=%{y:.2f} m<br>z=%{z:.2f} m<br>'+Q[key].label+'=%{marker.color:.3f} '+Q[key].unit+'<extra></extra>'};}
function shapes(plane){if(plane!=='burner')return[];return DATA.boxes.map(b=>({type:'rect',x0:b[0],x1:b[1],y0:b[4],y1:b[5],line:{color:'#8399a5',width:1},fillcolor:'#c6d2d8'}));}
function render(){const [cs,key,pn]=controls.map(x=>x.value),d=DATA.cases[cs],m=d.metrics,s=d.slices[pn],q=Q[key],range=DATA.ranges[key];document.getElementById('flow').innerHTML=fmt(m.fan_flow_m3_min,2)+' <span>m³/min</span>';document.getElementById('drop').innerHTML=fmt(m.network_makeup_pressure_drop_Pa,m.network_makeup_pressure_drop_Pa>10?1:3)+' <span>Pa</span>';document.getElementById('speed').innerHTML=fmt(m.occupied_zone_mean_speed_m_s,3)+' <span>m/s</span>';document.getElementById('temp').innerHTML=fmt(m.occupied_zone_mean_temperature_C,2)+' <span>°C</span>';document.getElementById('sliceTitle').textContent='02 · '+q.label+'切片';document.getElementById('sliceNote').textContent=DATA.labels[cs]+' · '+(s.axis==='y'?'垂直剖面 y = ':'水平剖面 z = ')+s.value.toFixed(2)+' m · 懸停讀取數值';document.getElementById('fieldNote').textContent=q.note;
let scene=d.model.concat([d.stream,planeTrace(s,key,range)]);let sl={...layoutBase,margin:{l:0,r:0,t:0,b:0},uirevision:'kitchen-camera',scene:{xaxis:{title:'x / m',range:[0,3.3],backgroundcolor:'#f2f6f8',showbackground:true},yaxis:{title:'y / m',range:[0,1.65],backgroundcolor:'#f2f6f8',showbackground:true},zaxis:{title:'z / m',range:[0,2.2],backgroundcolor:'#f2f6f8',showbackground:true},aspectmode:'data',camera},legend:{orientation:'h',x:0,y:-.03,font:{size:10}}};Plotly.react('scene',scene,sl,{responsive:true,displaylogo:false,toImageButtonOptions:{filename:'kitchen_3d_'+cs,scale:2}});
let hm={type:'heatmap',x:s.u,y:s.v,z:s[key],colorscale:q.colors,zmin:range[0],zmax:range[1],zsmooth:false,connectgaps:false,colorbar:{title:{text:q.unit},thickness:13,len:.8,x:1.02},hoverongaps:false,hovertemplate:'x=%{x:.2f} m<br>'+(s.axis==='y'?'z':'y')+'=%{y:.2f} m<br>'+q.label+'=%{z:.3f} '+q.unit+'<extra></extra>'};let ll={...layoutBase,margin:{l:52,r:70,t:30,b:50},xaxis:{title:'x / m',range:[0,3.3],constrain:'domain'},yaxis:{title:(s.axis==='y'?'z':'y')+' / m',range:[0,s.axis==='y'?2.2:1.65],scaleanchor:'x',scaleratio:1,constrain:'domain'},shapes:shapes(pn),showlegend:false};Plotly.react('slice',[hm,arrows(s)],ll,{responsive:true,displaylogo:false,toImageButtonOptions:{filename:'kitchen_'+key+'_'+cs,scale:2}});
document.getElementById('comparison').innerHTML=Object.keys(DATA.cases).map(c=>{let m=DATA.cases[c].metrics;return '<tr class="'+(c===cs?'active-row':'')+'"><td>'+DATA.labels[c]+'</td><td>'+fmt(m.fan_flow_m3_min,2)+'</td><td>'+fmt(m.network_makeup_pressure_drop_Pa,3)+'</td><td>'+fmt(m.occupied_zone_mean_speed_m_s,3)+'</td><td>'+fmt(m.occupied_zone_mean_temperature_C,2)+'</td><td>'+m.net_flux_relative_error.toExponential(2)+'</td><td>'+(m.converged?'已收斂':'未收斂')+'</td></tr>'}).join('');window.KITCHEN_READY=true;}
controls.forEach(el=>el.addEventListener('change',render));document.getElementById('reset').onclick=()=>Plotly.relayout('scene',{'scene.camera':camera});try{render()}catch(err){document.getElementById('error').style.display='block';document.getElementById('error').textContent='視覺化載入錯誤：'+err.message;throw err}
</script></body></html>'''

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--allow-coarse',action='store_true');args=ap.parse_args()
    config=json.loads((ROOT/'config.json').read_text());results=ROOT/'results';out=ROOT;fields={};data={};metrics={}
    for case in CASES:
        m=json.loads((results/(case+'_metrics.json')).read_text())
        if not m.get('converged'):raise ValueError(f'{case} is not converged')
        if m.get('heat_formulation')!='conservative_galerkin_edge_diffusion':raise ValueError(f'{case}: conservative thermal results are required')
        if not args.allow_coarse and abs(m['mesh_target_m']-.16)>1e-6:raise ValueError(f'{case}: waiting for final h=.16 results')
        metrics[case]=m;field=TetraField(dict(np.load(results/(case+'.npz'))));fields[case]=field
        print(f'Processing {case}: {len(field.p)} nodes',flush=True)
        stream,paths=streamline_data(field,case)
        data[case]=dict(metrics=m,model=model_traces(case),stream=stream,paths=paths,slices=slice_data(field))
    write_obj(out);write_probes(out,fields);static_plots(out,data,metrics)
    # Fixed common scale for fair visual comparisons, nodal extremes retained in data.
    ranges={'speed':[0,1.4],'temperature':[25,max(float(np.nanmax(data[c]['slices']['burner']['temperature'])) for c in CASES)],'pressure':[min(float(np.nanmin(fields[c].pressure)) for c in CASES),max(float(np.nanmax(fields[c].pressure)) for c in CASES)]}
    payload=dict(cases={c:{k:v for k,v in d.items() if k!='paths'} for c,d in data.items()},labels=LABELS,ranges=ranges,colormaps={k:get_colorscale(v) for k,v in [('speed','Turbo'),('temperature','Inferno'),('pressure','RdBu')]},boxes=[b for _,b,_ in BOXES])
    report=HTML.replace('__PLOTLY__',get_plotlyjs()).replace('__DATA__',json.dumps(serial(payload),ensure_ascii=False,separators=(',',':')).replace('</','<\\/')).replace('__ASSUMPTIONS__',''.join('<li>'+html.escape(x)+'</li>' for x in config['assumptions_zh']))
    (out/'interactive_3d.html').write_text(report,encoding='utf-8')
    (out/'visualization_metadata.json').write_text(json.dumps({'source_files':[str(Path('results')/(c+'.npz')) for c in CASES],'interpolation':'Barycentric interpolation on original FEM tetrahedra; missing/solid regions masked.','streamlines':'Midpoint integration of computed steady velocity, 0.01 m spatial step; direction only, not particle tracking.','common_display_ranges':ranges,'temperature_is':'Air temperature, not flame or pan temperature.','pressure_is':'FE room-side relative field; network makeup pressure drop is separate.'},ensure_ascii=False,indent=2))
    print('Visualization complete',flush=True)
if __name__=='__main__':main()
