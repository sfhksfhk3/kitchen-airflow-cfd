from pathlib import Path
import json,numpy as np,pandas as pd
import plotly.offline as po
from plotly.colors import sequential
from process_ceiling import ROOT,OUT,size_weights,DS
from foam_mesh import mesh_geometry
from visualize_ceiling import main_cases,label

def main(output_path=None, case_names=None, scenario_text=None):
 data={}
 cases=case_names if case_names is not None else main_cases()+[f'{l}_closed_stick{a}' for l in ['corner','center'] for a in ['0.1','0.5']]
 for case in cases:
  if not (OUT/case/'summary.json').exists():continue
  f=pd.read_csv(OUT/case/('ceiling_complete_event_faces.csv' if (OUT/case/'ceiling_complete_event_faces.csv').exists() else 'ceiling_faces.csv'));s=json.loads((OUT/case/'summary.json').read_text());g=mesh_geometry(ROOT/'runs'/case);ids=f.original_face_id.to_numpy(dtype=int)
  used=np.unique(np.concatenate([g['faces'][i] for i in ids]));index={int(i):k for k,i in enumerate(used)};tri=[];which=[]
  for k,faceid in enumerate(ids):
   vv=[index[int(i)] for i in g['faces'][faceid]]
   for j in range(1,len(vv)-1):tri.append([vv[0],vv[j],vv[j+1]]);which.append(k)
  fields=f[[f'M_{d:g}um_unit_mg_m2_per_mg' for d in DS]].to_numpy();values={str(gsd):((fields@size_weights(gsd))[which]).tolist() for gsd in [1.2,1.5,2.]}
  facevalues={str(gsd):(fields@size_weights(gsd)).tolist() for gsd in [1.2,1.5,2.]}
  data[case]=dict(face_centres=f[['x_m','y_m','z_m']].to_numpy().tolist(),face_ids=ids.tolist(),Mfaces=facevalues,label=label(case)+(('・附著率'+case.split('stick')[1]) if 'stick' in case else '・附著率1'),points=g['points'][used].tolist(),triangles=tri,M=values,area_m2=f.area_m2.tolist(),faceM=fields.tolist(),source_seconds=120,observe_seconds=240)
 payload=json.dumps(data,ensure_ascii=False,separators=(',',':'))
 html='''<!doctype html><html lang="zh-Hant"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>天花油脂沉積・離線互動3D</title><style>
 body{margin:0;background:#f4f6f8;color:#18222e;font:16px system-ui,sans-serif}main{max-width:1200px;margin:20px auto;padding:0 20px}h1{font-size:25px}section{padding:16px;background:white;border-radius:8px;margin:12px 0}.controls{display:flex;flex-wrap:wrap;gap:16px}label{display:flex;flex-direction:column;font-size:14px;gap:6px}input,select{padding:8px;font:inherit;max-width:320px}#plot{height:65vh;min-height:450px}#stats{font-size:17px;font-weight:600}.note{font-size:14px;line-height:1.65}button{padding:8px 12px;cursor:pointer}</style>
 <main><h1>天花油脂沉積・離線互動3D</h1><p><strong style="color:#a22">驗證未收斂：後段粗／細分支相差約27倍，不能作材料積油量設計值。</strong></p><p>沿用原OpenFOAM網格與瞬態流場；顯示油滴真正離開空氣、附著在天花的質量。</p>
 <section class="controls"><label>配置與附著率<select id="case"></select></label><label>假設質量分布 GSD<select id="gsd"><option>1.2</option><option selected>1.5</option><option value="2.0">2.0</option></select></label><label>每次排放總油滴（mg）<input id="dose" type="number" value="1" min="0" step="0.1"></label><label>累積日數<input id="days" type="number" value="1" min="1" max="3650"></label><label>每日煮食次數<select id="events"><option>1</option><option>2</option></select></label><label>清潔情景<select id="clean"><option value="0">不清潔</option><option value="30">每30日清潔90%（假設）</option></select></label><label>觀看方式<button id="top">天花俯視</button><button id="angle">3D斜視</button></label></section>
 <section><div id="stats"></div><div id="plot"></div></section>
 <section class="note">顏色為 log10(mg/m²)，滑過面中心可查看該面的計算質量。每次排放預設1 mg，避免把未知排放率當成精確實測。參考情境為兩爐合計14.4 mg/min、排放120秒，即28.8 mg；僅為商用文獻外推，非這個住家的量測。每次均沿用120秒排放、排風觀測至240秒，另加末端固定流場的殘留油滴清空積分。長期模型假設重複相同事件、清潔外無自然損失；殘留油滴經輸運方程分配到抽走、天花和其他表面，沒有全數算成天花沉積。<br>沒有氣態有機物冷凝、黃變、油膜厚度或漆面老化模型。三種材料的真實附著率未校準，不能從這幅圖排出材料優劣。灰色佔用區不是零油污實測。低於色階下限仍可能有正沉積；色階下限不是檢出限。<br>這是近壁閉合敏感度研究，原網格嚴格檢查仍有未通過項目；請搭配中文報告解讀。</section></main><script>__PLOTLY__</script><script>
 const COLORS=__COLORS__;const DATA=__DATA__;const sel=document.querySelector('#case');Object.entries(DATA).forEach(([k,v])=>sel.add(new Option(v.label,k)));sel.value=DATA['corner_both_closed_parallel']?'corner_both_closed_parallel':(DATA['corner_both_closed_T33_Q30']?'corner_both_closed_T33_Q30':Object.keys(DATA)[0]);
 function units(days,events,interval){let stock=0;for(let d=1;d<=days;d++){stock+=events;if(interval&&d%interval===0)stock*=0.1}return stock}
 function wire(){let x=[],y=[],z=[];const edges=[[[0,0,0],[3.3,0,0]],[[0,1.65,0],[3.3,1.65,0]],[[0,0,2.2],[3.3,0,2.2]],[[0,1.65,2.2],[3.3,1.65,2.2]],[[0,0,0],[0,1.65,0]],[[3.3,0,0],[3.3,1.65,0]],[[0,0,2.2],[0,1.65,2.2]],[[3.3,0,2.2],[3.3,1.65,2.2]],...[[0,0],[3.3,0],[0,1.65],[3.3,1.65]].map(p=>[[...p,0],[...p,2.2]])];for(const e of edges){for(const p of e){x.push(p[0]);y.push(p[1]);z.push(p[2])}x.push(null);y.push(null);z.push(null)}return {type:'scatter3d',x,y,z,mode:'lines',line:{color:'#8d99a5',width:2},hoverinfo:'skip',showlegend:false}}
 function render(){const k=sel.value,v=DATA[k],gsd=Number(document.querySelector('#gsd').value).toFixed(1);const dose=Math.max(0,Number(document.querySelector('#dose').value)||0),days=Math.max(1,Math.min(3650,Math.floor(Number(document.querySelector('#days').value)||1))),events=Number(document.querySelector('#events').value),interval=Number(document.querySelector('#clean').value);const mult=dose*units(days,events,interval),m=(v.M[gsd]||v.M['2.0']||v.M['2']).map(a=>a*mult);const actualMax=Math.max(...m,0),vmax=Math.max(actualMax,1e-30),floor=vmax*1e-8;const vals=m.map(x=>Math.log10(Math.max(x,floor)));const tr=v.triangles;
 let counts=new Map;for(let i=0;i<tr.length;i++){const q=tr[i];counts.set(i,m[i])}const mesh={type:'mesh3d',x:v.points.map(x=>x[0]),y:v.points.map(x=>x[1]),z:v.points.map(x=>x[2]),i:tr.map(x=>x[0]),j:tr.map(x=>x[1]),k:tr.map(x=>x[2]),intensity:vals,intensitymode:'cell',customdata:m,colorscale:COLORS,cmin:Math.log10(floor),cmax:Math.log10(vmax),colorbar:{title:{text:'log10(mg/m²)'}},flatshading:true,lighting:{ambient:1,diffuse:0},hoverinfo:'skip',showscale:true};
 const fmass=(v.Mfaces[gsd]||v.Mfaces['2.0']).map(x=>x*mult);const hover={type:'scatter3d',mode:'markers',x:v.face_centres.map(x=>x[0]),y:v.face_centres.map(x=>x[1]),z:v.face_centres.map(x=>x[2]+0.0001),marker:{size:4,opacity:0.001,color:'white'},customdata:fmass.map((x,i)=>[v.face_ids[i],x]),showlegend:false,hovertemplate:'原面ID=%{customdata[0]}<br>x=%{x:.2f} m，y=%{y:.2f} m<br>沉積=%{customdata[1]:.3e} mg/m²<extra>面中心</extra>'};
 Plotly.react('plot',[mesh,wire(),hover],{margin:{l:0,r:0,t:20,b:0},scene:{aspectmode:'data',xaxis:{title:'x (m)',range:[0,3.3]},yaxis:{title:'y (m)',range:[0,1.65]},zaxis:{title:'z (m)',range:[0,2.3]},camera:{eye:{x:0,y:0,z:2.5},up:{x:0,y:1,z:0},projection:{type:'orthographic'}}},uirevision:k},{responsive:true,displaylogo:false});const total=fmass.reduce((acc,x,i)=>acc+x*v.area_m2[i],0),area=v.area_m2.reduce((a,b)=>a+b,0);const ordered=fmass.map((x,i)=>[x,v.area_m2[i]]).sort((a,b)=>a[0]-b[0]);let cumulative=0,p95=0;for(const [x,a] of ordered){cumulative+=a;p95=x;if(cumulative>=.95*area)break}document.querySelector('#stats').textContent=v.label+'｜'+days+'日 × 每日'+events+'次｜每次排放'+dose+' mg｜總沉積 '+total.toExponential(3)+' mg｜平均 '+(total/area).toExponential(3)+' mg/m²｜P95 '+p95.toExponential(3)+' mg/m²｜最大 '+actualMax.toExponential(3)+' mg/m²';}
 for(const id of ['case','gsd','dose','days','events','clean'])document.querySelector('#'+id).addEventListener('change',render);document.querySelector('#top').onclick=()=>Plotly.relayout('plot',{'scene.camera':{eye:{x:0,y:0,z:2.5},up:{x:0,y:1,z:0},projection:{type:'orthographic'}}});document.querySelector('#angle').onclick=()=>Plotly.relayout('plot',{'scene.camera':{eye:{x:1.5,y:-1.8,z:1.7},up:{x:0,y:0,z:1},projection:{type:'perspective'}}});render();</script></html>'''
 scale=[[i/(len(sequential.Inferno)-1),c] for i,c in enumerate(sequential.Inferno)]
 html=html.replace('__PLOTLY__',po.get_plotlyjs()).replace('__DATA__',payload).replace('__COLORS__',json.dumps(scale))
 if scenario_text:
  html=html.replace('驗證未收斂：後段粗／細分支相差約27倍，不能作材料積油量設計值。',scenario_text['warning'])
  html=html.replace('沿用原OpenFOAM網格與瞬態流場；顯示油滴真正離開空氣、附著在天花的質量。',scenario_text['intro'])
  html=html.replace('參考情境為兩爐合計14.4 mg/min、排放120秒，即28.8 mg；僅為商用文獻外推，非這個住家的量測。',scenario_text['reference'])
 (Path(output_path) if output_path else ROOT/'ceiling_deposition_interactive_3d.html').write_text(html)
 print('interactive cases',len(data),'bytes',len(html.encode()))
if __name__=='__main__':main()
