"""Square scientific figures and a self-contained, offline 3D data viewer."""
from pathlib import Path
import json,csv,numpy as np
from scipy.spatial import cKDTree
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.patches import Rectangle
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from plotly.offline import get_plotlyjs
from emission_scenario import config,weights,DIAMETERS
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'deliverables';FIG=OUT/'figures_1to1';FIG.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'Noto Sans CJK JP','axes.unicode_minus':False,'font.size':11,'savefig.facecolor':'white'})
DOORS=['both_open','both_closed','left_open','right_open'];LABELS=['兩門全開','兩門全關','左開右關','右開左關'];COLORS=['#208c7e','#db5b43','#427db3','#8b68ab']
LAYOUTS=['corner','center'];LN={'corner':'靠左','center':'置中'}

def save(fig,name):
    fig.savefig(FIG/(name+'.png'),dpi=300);fig.savefig(FIG/(name+'.pdf'));fig.savefig(FIG/(name+'.svg'));plt.close(fig)

def frames(name):return dict(np.load(OUT/'data'/name/'visualization_frames.npz'))

def geometry_figure():
    fig=plt.figure(figsize=(12,12));fig.subplots_adjust(top=.91,bottom=.08,hspace=.12)
    def box(ax,a,b,color,alpha=.3):
        p=np.array([[a[0],a[1],a[2]],[b[0],a[1],a[2]],[b[0],b[1],a[2]],[a[0],b[1],a[2]],[a[0],a[1],b[2]],[b[0],a[1],b[2]],[b[0],b[1],b[2]],[a[0],b[1],b[2]]])
        ff=[[0,1,2,3],[4,5,6,7],[0,1,5,4],[1,2,6,5],[2,3,7,6],[3,0,4,7]]
        ax.add_collection3d(Poly3DCollection(p[ff],facecolors=color,edgecolors='#677a85',linewidths=.5,alpha=alpha))
    for j,lay in enumerate(LAYOUTS):
        ax=fig.add_subplot(2,1,j+1,projection='3d');hx=.7 if lay=='corner' else 1.55
        box(ax,[0,0,0],[.7,.6,2.2],'#7b8b96',.35);box(ax,[.7,0,0],[3.3,.6,.9],'#b1c0c9',.28);box(ax,[.7,0,1.6],[3.3,.3,2.2],'#c2cdd3',.2)
        box(ax,[hx,0,1.596],[hx+.9,.3,1.605],'#268caa',.6)
        q=.1/np.sqrt(2);b=np.array([[hx,.3,1.6],[hx+.9,.3,1.6],[hx+.9,.3+q,1.6-q],[hx,.3+q,1.6-q]])
        ax.add_collection3d(Poly3DCollection([b],facecolors='#2f7085',alpha=.8))
        ax.scatter([hx+.25,hx+.65],[.3,.3],[.9,.9],c='#dc7638',s=65)
        for x in [0,3.3]:ax.plot([x]*5,[.7,1.6,1.6,.7,.7],[0,0,2.2,2.2,0],c='#258e71',lw=1.7)
        ax.set(xlim=(0,3.3),ylim=(0,1.65),zlim=(0,2.2),xlabel='x (m)',ylabel='y (m)',zlabel='z (m)',title=f'{LN[lay]}｜煙機及爐頭整體中心 x={hx+.45:g} m');ax.set_box_aspect([3.3,1.65,2.2]);ax.view_init(elev=22,azim=55)
    fig.suptitle('計算幾何｜同一廚房、兩種灶台位置',fontsize=19)
    fig.text(.08,.025,'煙機 90×30 cm，距灶面 70 cm；前緣擋板長 10 cm、向下 45°，僅在煙機寬度內。\n圖示門洞為全開範圍；關門時僅保留底部 5 mm 開口。',fontsize=10)
    save(fig,'00_geometry_comparison_1to1')

def plane(a,tt,z,cols):
    t=a['times_s'];j=np.argmin(abs(t-tt));p=a['centres'];x=np.linspace(0,3.3,265);y=np.linspace(0,1.65,133);X,Y=np.meshgrid(x,y)
    dist,ix=cKDTree(p).query(np.c_[X.ravel(),Y.ravel(),np.full(X.size,z)])
    v=(a['values'][j,:,cols] if np.isscalar(cols) else a['values'][j][:,cols]@np.array(weights())*240)[ix].reshape(X.shape)
    mask=((X<.7)&(Y<.6)) if z>.9 else Y<.6
    v=np.ma.array(v,mask=mask);return X,Y,v,float(t[j])

def oil_figures(available):
    for door,label in zip(DOORS,LABELS):
        if not all(f'{lay}_{door}' in available for lay in LAYOUTS):continue
        fig,axes=plt.subplots(2,1,figsize=(12,12));fig.subplots_adjust(top=.88,bottom=.16,left=.08,right=.96,hspace=.3)
        for ax,lay in zip(axes,LAYOUTS):
            a=frames(f'{lay}_{door}');j=np.flatnonzero(a['target_times_s']==120)[0];tt=a['times_s'][j]
            X,Y,c,actual=plane(a,tt,1.5,list(range(3,8)))
            image=ax.pcolormesh(X,Y,np.maximum(c,1e-5),norm=LogNorm(1e-4,3000),cmap='magma',shading='auto',rasterized=True)
            ax.add_patch(Rectangle((0,0),.7,.6,facecolor='#c7ced1',edgecolor='#56646c'))
            ax.add_patch(Rectangle((.7,0),2.6,.6,fill=False,ls='--',lw=1,ec='#92a3ac'))
            hx=.7 if lay=='corner' else 1.55
            ax.add_patch(Rectangle((hx,0),.9,.3,fill=False,lw=1.5,ec='#47bec0'))
            ax.scatter([hx+.25,hx+.65],[.3,.3],marker='x',c='#55c83c',s=55,label='爐頭投影')
            for xx,op in [(0,door in ['both_open','left_open']),(3.3,door in ['both_open','right_open'])]:ax.plot([xx,xx],[.7,1.6],c='#27936c' if op else '#a54c43',lw=6)
            ax.set(xlim=(0,3.3),ylim=(0,1.65),xlabel='沿檯面 x (m)',ylabel='距後牆 y (m)',title=f'{LN[lay]}｜z=1.50 m；實際快照 {actual:g} s',aspect='equal')
            ax.legend(loc='upper right',fontsize=9,framealpha=.85)
        fig.suptitle(f'中式炒菜油滴質量濃度｜{label}',fontsize=19,y=.965)
        cb=fig.colorbar(image,ax=axes,orientation='horizontal',fraction=.035,pad=.12,aspect=45,extend='both');cb.set_label('參考情境油滴質量濃度 (µg/m³，對數色階；非 PM2.5)')
        fig.text(.07,.025,'兩爐共 4000 W；油滴參考排放共 14.4 mg/min。分布中位徑 2.5 µm，GSD 1.5 為假設。\n場快照取保存資料中最接近 120 s 的烹調中時刻；定量比較採 90–120 s 時序統計。',fontsize=10)
        save(fig,f'0{DOORS.index(door)+1}_oil_{door}_1to1')

def comparison_figures(summary):
    by={c['case']:c for c in summary['cases']}
    if len(by)<8:return
    fig,axes=plt.subplots(2,1,figsize=(12,12));fig.subplots_adjust(hspace=.34,bottom=.14,top=.91)
    xx=np.arange(4)
    for lay,off,col in [('corner',-.18,'#356c91'),('center',.18,'#d47b48')]:
        vals=[by[f'{lay}_{d}'] for d in DOORS]
        axes[0].bar(xx+off,[v['retention_integral_observed_s'] for v in vals],.34,label=LN[lay],color=col)
        axes[1].bar(xx+off,[v['walkway_90_120_oil_ug_m3'] for v in vals],.34,label=LN[lay],color=col)
    axes[0].set_ylabel('觀測滯留積分 (s)');axes[0].set_title('釋放後至 240 s 的存量積分／總釋放量（完整平均滯留時間的下限）')
    axes[1].set_ylabel('走道平均油滴 (µg/m³)');axes[1].set_yscale('log');axes[1].set_title('烹調末段 90–120 s，檯面前方走道平均濃度')
    for ax in axes:ax.set_xticks(xx,LABELS);ax.grid(axis='y',alpha=.2);ax.legend()
    fig.suptitle('位置與門況：全室滯留、局部走道濃度',fontsize=19)
    fig.text(.08,.035,'走道區：y>0.6 m，z=0.9–1.8 m；尚有殘留時，240 s 的滯留積分不是完整清除時間。\nµg/m³為文獻參考排放情境的總油滴質量，非住家實測／PM2.5。',fontsize=10)
    save(fig,'05_inventory_exposure_comparison_1to1')
    for normalized in [False,True]:
        fig,axes=plt.subplots(2,2,figsize=(12,12));fig.subplots_adjust(hspace=.27,wspace=.23,bottom=.16,top=.91)
        for col,lay in enumerate(LAYOUTS):
            for door,label,color in zip(DOORS,LABELS,COLORS):
                rows=list(csv.DictReader((OUT/'data'/f'{lay}_{door}'/'mixture_history_reference_scenario.csv').open()));t=np.array([float(r['time_s']) for r in rows]);sel=t>=120
                for row,key in [(0,'airborne_oil_mg'),(1,'walkway_oil_ug_m3')]:
                    y=np.array([float(r[key]) for r in rows]);value=y[sel]/max(np.interp(120,t,y),1e-30) if normalized else y[sel]
                    axes[row,col].semilogy(t[sel]-120,np.maximum(value,1e-16),label=label,color=color,lw=1.6)
            axes[0,col].set_title(LN[lay]);axes[1,col].set_xlabel('停止放熱與排放後時間 (s)')
        for ax in axes.flat:ax.grid(alpha=.2);ax.set_xlim(0,120)
        axes[0,0].set_ylabel('全室存量 / 停止時存量' if normalized else '全室空氣油滴存量 (mg)')
        axes[1,0].set_ylabel('走道濃度 / 停止時濃度' if normalized else '走道平均油滴 (µg/m³)')
        for ax in axes[0]:ax.set_ylim((1e-8,2) if normalized else (1e-8,10))
        for ax in axes[1]:ax.set_ylim((1e-8,3) if normalized else (1e-8,300))
        handles,labels=axes[0,0].get_legend_handles_labels();fig.legend(handles,labels,loc='lower center',ncol=4,bbox_to_anchor=(.5,.075))
        fig.suptitle('關火後油滴衰減｜'+('各自起點歸一化' if normalized else '實際量級的參考情境'),fontsize=19)
        fig.text(.07,.025,'120 s 時同時停止放熱及油滴排放，保持排風。較低曲線可由煙機排出、門外流出或沉積造成。\n歸一化曲線各自起點不同；局部比例更低不代表絕對濃度更低。',fontsize=10)
        save(fig,'07_normalized_clearance_1to1' if normalized else '06_clearance_four_panels_1to1')
    fig,axes=plt.subplots(2,1,figsize=(12,12));fig.subplots_adjust(hspace=.3,bottom=.12,top=.91)
    for ax,lay in zip(axes,LAYOUTS):
        for door,label,color in zip(DOORS,LABELS,COLORS):
            b=sorted([b for b in summary['size_cases'] if b['case']==f'{lay}_{door}'],key=lambda b:b['diameter_um'])
            ax.loglog(DIAMETERS,[100*v['deposited_fraction_end'] for v in b],'-o',label=label,color=color)
        ax.set(title=LN[lay]+'｜沉積占釋放量的比例',xlabel='獨立計算粒徑 (µm)',ylabel='240 s 累積沉積 (%)');ax.set_xticks(DIAMETERS,[str(x) for x in DIAMETERS]);ax.grid(alpha=.2);ax.legend()
    fig.suptitle('粒徑效應：壁面沉積',fontsize=19);fig.text(.08,.025,'各粒徑為獨立相同源強試驗。黏附率=1，沒有凝結、再懸浮及慣性撞擊；沉積率受近壁模型影響。',fontsize=10)
    save(fig,'08_size_deposition_1to1')

HTML=r'''<!doctype html><html lang="zh-Hant"><meta charset="utf-8"><title>廚房油滴 OpenFOAM 3D</title>
<style>body{font-family:system-ui,"Noto Sans CJK TC",sans-serif;margin:0;color:#193044;background:#f3f6f8}header{padding:20px 28px;background:#123448;color:white}h1{font-size:23px;margin:0 0 6px}header p{margin:0;font-size:13px;opacity:.9}.controls{padding:14px 24px;display:flex;gap:15px;flex-wrap:wrap;align-items:center;background:white;border-bottom:1px solid #d7e1e6}label{font-size:13px;display:flex;gap:6px;align-items:center}select,input{padding:6px;border:1px solid #a5bac5;border-radius:5px}input[type=number]{width:75px}.plots{display:grid;grid-template-columns:1fr 1fr;gap:12px;padding:12px}.card{background:white;border-radius:8px;overflow:hidden}.chart{height:620px}.metrics{font-size:13px;padding:12px;border-top:1px solid #dce5eb;line-height:1.8}.notes{margin:5px 24px 25px;padding:15px;background:white;font-size:13px;line-height:1.8}.badge{background:#e9f5f0;padding:2px 8px;border-radius:5px}@media(max-width:900px){.plots{grid-template-columns:1fr}.chart{height:500px}}</style>
<header><h1>中式炒菜油滴與熱浮力｜靠左、置中對照</h1><p>OpenFOAM 2406 · 瞬態熱流與五種粒徑 · 可離線旋轉、縮放、查閱單元資料</p></header>
<div class="controls"><label>門況<select id="door"><option value="both_open">兩門全開</option><option value="both_closed">兩門全關</option><option value="left_open">左開右關</option><option value="right_open">右開左關</option></select></label>
<label>場量<select id="field"><option value="oil">油滴質量濃度</option><option value="temperature">空氣溫度</option><option value="velocity">風速</option><option value="pressure">相對縮減壓力</option></select></label>
<label>粒徑<select id="size"><option value="mix">參考質量分布</option><option value="0">0.1 µm</option><option value="1">0.3 µm</option><option value="2">1 µm</option><option value="3">3 µm</option><option value="4">10 µm</option></select></label>
<label>兩爐總排放<input id="rate" type="number" min="0.01" max="1000" step="0.1" value="14.4">mg/min</label><label>分布GSD<input id="gsd" type="number" min="1.1" max="3" step="0.1" value="1.5"></label>
<label>時間<input id="time" type="range" min="0" max="9" value="5" step="1"><span id="timeLabel"></span></label><button id="reset">重設視角</button></div>
<div class="plots"><div class="card"><div id="corner" class="chart"></div><div id="cornerMetrics" class="metrics"></div></div><div class="card"><div id="center" class="chart"></div><div id="centerMetrics" class="metrics"></div></div></div>
<div class="notes"><span class="badge">參考情境，非住家實測</span> 預設每爐頭油滴排放7.2 mg/min，質量中位徑2.5 µm，GSD1.5為假設；修改排放量會按稀薄顆粒模型線性縮放。選單粒徑時表示將相同總排放集中到該粒徑的敏感性情境。<br>圖中資料來自保存的單元中心場，抽樣顯示；每圖標示實際保存時刻。120 s停止加熱與排放，排風持續至240 s。µg/m³代表總油滴質量，不能當成PM2.5或超細粒子數。<br>計算涵蓋重力、布朗擴散、熱泳與黏附沉積；尚無火焰、鍋具／板材導熱、油滴凝結／蒸發／凝併及室外風壓。5mm門縫兩層基準網格；數值敏感性與限制見中文報告。壓力不包含完整風機、煙管或濾網。<br>走道區為廚房內y&gt;0.6m、z=0.9–1.8m，不是門外走廊末端。全室滯留與局部濃度衡量不同；沉積也不等於表面清潔。</div>
<script>__PLOTLY__</script><script>const DATA=__DATA__;
const TICKS=[0,10,30,60,100,120,130,150,180,240],D=[.1,.3,1,3,10];let cameras={};
function erf(x){let s=Math.sign(x),a=Math.abs(x),t=1/(1+.3275911*a);return s*(1-(((((1.061405429*t-1.453152027)*t)+1.421413741)*t-.284496736)*t+.254829592)*t*Math.exp(-a*a))}
function massWeights(){let g=Math.max(1.01,+document.getElementById('gsd').value),e=[0,Math.sqrt(.03),Math.sqrt(.3),Math.sqrt(3),Math.sqrt(30),Infinity];let c=e.map((v,i)=>i===0?0:i===5?1:.5*(1+erf(Math.log(v/2.5)/Math.log(g)/Math.SQRT2)));return c.slice(1).map((v,i)=>Math.max(0,v-c[i]));}
const F=[[0,1,2],[0,2,3],[4,6,5],[4,7,6],[0,4,5],[0,5,1],[1,5,6],[1,6,2],[2,6,7],[2,7,3],[3,7,4],[3,4,0]];
function box(a,b,color,opacity=.25){let p=[[a[0],a[1],a[2]],[b[0],a[1],a[2]],[b[0],b[1],a[2]],[a[0],b[1],a[2]],[a[0],a[1],b[2]],[b[0],a[1],b[2]],[b[0],b[1],b[2]],[a[0],b[1],b[2]]];return {type:'mesh3d',x:p.map(v=>v[0]),y:p.map(v=>v[1]),z:p.map(v=>v[2]),i:F.map(v=>v[0]),j:F.map(v=>v[1]),k:F.map(v=>v[2]),color,opacity,hoverinfo:'skip',showlegend:false};}
function geometry(lay,door){let hx=lay==='corner'?.7:1.55,tip=.1/Math.sqrt(2);let out=[box([0,0,0],[.7,.6,2.2],'#697b85',.45),box([.7,0,0],[3.3,.6,.9],'#b2bfc5',.3),box([.7,0,1.6],[3.3,.3,2.2],'#b6c2c7',.2),box([hx,0,1.595],[hx+.9,.3,1.61],'#3587a1',.6),box([0,0,-.008],[3.3,1.65,0],'#b8c7ce',.12)];out.push({type:'mesh3d',x:[hx,hx+.9,hx+.9,hx],y:[.3,.3,.3+tip,.3+tip],z:[1.6,1.6,1.6-tip,1.6-tip],i:[0,0],j:[1,2],k:[2,3],color:'#417e8b',opacity:.65,hoverinfo:'skip',showlegend:false});out.push({type:'scatter3d',mode:'markers',x:[hx+.25,hx+.65],y:[.3,.3],z:[.9,.9],marker:{size:6,color:'#df6f36',symbol:'circle'},text:['爐頭1','爐頭2'],hovertemplate:'%{text}<extra></extra>',showlegend:false});for(let [x,op] of [[0,['both_open','left_open'].includes(door)],[3.3,['both_open','right_open'].includes(door)]])out.push(box([x-.002,.7,op?0:.005],[x+.002,1.6,2.2],op?'#2ba879':'#b95443',op?.06:.17));return out;}
function draw(){let door=document.getElementById('door').value,field=document.getElementById('field').value,size=document.getElementById('size').value,target=TICKS[+document.getElementById('time').value],rate=Math.max(0,+document.getElementById('rate').value)*1000/60,w=size==='mix'?massWeights():D.map((v,i)=>i===+size?1:0);document.getElementById('timeLabel').textContent=target+' s';for(let lay of ['corner','center']){let item=DATA.cases[lay+'_'+door];if(!item)continue;let j=item.targetTimes.reduce((best,v,i)=>Math.abs(v-target)<Math.abs(item.targetTimes[best]-target)?i:best,0),v=item.frames[j],xx=[],yy=[],zz=[],val=[],custom=[];let unit='',cmin=0,cmax=1,scale='Viridis';for(let i=0;i<item.xyz.length;i++){let row=v[i],value;if(field==='oil'){value=row.slice(3,8).reduce((s,c,k)=>s+c*w[k],0)*rate;if(value<1e-5)continue;}else value=row[field==='temperature'?0:field==='velocity'?1:2];xx.push(item.xyz[i][0]);yy.push(item.xyz[i][1]);zz.push(item.xyz[i][2]);val.push(field==='oil'?Math.log10(Math.max(1e-8,value)):field==='velocity'?Math.log10(Math.max(.001,value)):value);custom.push(value);}if(field==='oil'){cmin=-4;cmax=3.5;scale='Magma';unit='µg/m³'}if(field==='temperature'){cmin=25;cmax=180;scale='Inferno';unit='°C'}if(field==='velocity'){cmin=-3;cmax=1.5;unit='m/s'}if(field==='pressure'){cmin=-600;cmax=5;scale='RdBu';unit='Pa'}let bar={title:unit,len:.68};if(field==='oil'){bar.tickvals=[-4,-2,0,2,3];bar.ticktext=['0.0001','0.01','1','100','1000']}if(field==='velocity'){bar.tickvals=[-3,-2,-1,0,1];bar.ticktext=['.001','.01','.1','1','10']}let trace={type:'scatter3d',mode:'markers',x:xx,y:yy,z:zz,marker:{size:2.6,color:val,cmin,cmax,colorscale:DATA.scales[scale]||scale,opacity:.65,colorbar:bar},customdata:custom,hovertemplate:'x=%{x:.3f} m<br>y=%{y:.3f} m<br>z=%{z:.3f} m<br>數值=%{customdata:.4g} '+unit+'<extra></extra>',showlegend:false};let title=(lay==='corner'?'靠左':'置中')+'｜保存場 '+Number(item.times[j].toFixed(3))+' s';let layout={title:{text:title,font:{size:17}},margin:{l:0,r:0,b:0,t:45},paper_bgcolor:'white',uirevision:lay,scene:{xaxis:{title:'x (m)',range:[0,3.3]},yaxis:{title:'y (m)',range:[0,1.65]},zaxis:{title:'z (m)',range:[0,2.2]},aspectmode:'data',camera:cameras[lay]||{eye:{x:1.8,y:2.5,z:1.4}}}};let traces=[...geometry(lay,door),trace];if(field==='velocity'){let cx=[],cy=[],cz=[],cu=[],cv=[],cw=[];for(let i=0;i<item.xyz.length;i+=Math.ceil(item.xyz.length/130)){let row=v[i],sp=row[1];if(sp<.02||row.length<11)continue;cx.push(item.xyz[i][0]);cy.push(item.xyz[i][1]);cz.push(item.xyz[i][2]);cu.push(row[8]/sp);cv.push(row[9]/sp);cw.push(row[10]/sp)}traces.push({type:'cone',x:cx,y:cy,z:cz,u:cu,v:cv,w:cw,sizemode:'absolute',sizeref:.10,anchor:'center',colorscale:[[0,'#1b3749'],[1,'#1b3749']],showscale:false,hoverinfo:'skip',opacity:.8})}Plotly.react(lay,traces,layout,{responsive:true,displaylogo:false,toImageButtonOptions:{format:'png',width:1200,height:1200,scale:2}});let b=item.bins,mean=b.reduce((s,v,i)=>s+w[i]*v.walkway_90_120_c_over_R_s_m3,0)*rate,dep=b.reduce((s,v,i)=>s+w[i]*v.deposited_fraction_end,0)*100,hood=b.reduce((s,v,i)=>s+w[i]*v.hood_fraction_at_end,0)*100;document.getElementById(lay+'Metrics').innerHTML='烹調末段走道平均：<b>'+mean.toPrecision(4)+' µg/m³</b><br>觀測期間煙機排出 '+hood.toFixed(3)+'% · 沉積 '+dep.toFixed(3)+'%<br><small>指標採完整時序；圖形取已保存場。單元抽樣顯示，原始CSV／VTU保留完整資料。風速圖箭頭僅示方向。</small>';}}
for(let id of ['door','field','size','rate','gsd','time'])document.getElementById(id).addEventListener('input',draw);document.getElementById('reset').onclick=()=>{cameras={};for(let l of ['corner','center'])Plotly.relayout(l,{'scene.camera':{eye:{x:1.8,y:2.5,z:1.4}}})};draw();for(let l of ['corner','center']){let g=document.getElementById(l);if(g.on)g.on('plotly_relayout',e=>{if(e['scene.camera'])cameras[l]=e['scene.camera']});}
</script></html>'''

def interactive(summary,available):
    payload={'cases':{},'scenario':config(),'scales':{name:[[float(t),matplotlib.colors.to_hex(plt.get_cmap(name.lower())(t))] for t in np.linspace(0,1,17)] for name in ['Magma','Inferno']}};rng=np.random.default_rng(20261015)
    for name in available:
        a=frames(name);p=a['centres'];n=len(p);sample=rng.choice(n,min(2300,n),replace=False,p=a['volumes']/a['volumes'].sum())
        special=np.flatnonzero((p[:,2]<.006)|((p[:,1]<.65)&(p[:,2]>.9)&(p[:,2]<1.65)))
        if len(special)>1600:special=rng.choice(special,1600,replace=False)
        ids=np.unique(np.r_[sample,special]);v=a['values'][:,ids,:]
        # Five significant figures retain small concentrations for log display.
        compact=np.array([float(f'{x:.5g}') for x in v.ravel()]).reshape(v.shape)
        payload['cases'][name]={'xyz':p[ids].round(5).tolist(),'times':a['times_s'].round(8).tolist(),'targetTimes':a['target_times_s'].tolist(),'frames':compact.tolist(),'bins':sorted([b for b in summary['size_cases'] if b['case']==name],key=lambda b:b['diameter_um'])}
    html=HTML.replace('__PLOTLY__',get_plotlyjs()).replace('__DATA__',json.dumps(payload,ensure_ascii=False,separators=(',',':'),allow_nan=False))
    (OUT/'kitchen_oil_openfoam_interactive_3d.html').write_text(html)

if __name__=='__main__':
    summary=json.loads((OUT/'summary.json').read_text());available=[p.name for p in (OUT/'data').glob('*') if (p/'visualization_frames.npz').exists()]
    geometry_figure();oil_figures(available);comparison_figures(summary);interactive(summary,available)
    print('Rendered',len(available),'cases',flush=True)
