"""Publish only the new isothermal-ceiling/prescribed-flow scenario outputs."""
from pathlib import Path
import json,sys,re,base64,zipfile,xml.etree.ElementTree as ET
import numpy as np,pandas as pd,matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.colors import LogNorm
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import font_manager
from markdown_it import MarkdownIt
from playwright.sync_api import sync_playwright
from foam_mesh import mesh_geometry
from process_ceiling import weighted_quantile,DS,size_weights
from long_term import cycle_units
ROOT=Path(__file__).resolve().parents[1];S=ROOT/'scenarios/ceiling33_Q30';O=S/'results';F=S/'figures';F.mkdir(exist_ok=True)
FONT='/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc';font_manager.fontManager.addfont(FONT);plt.rcParams.update({'font.family':font_manager.FontProperties(fname=FONT).get_name(),'axes.unicode_minus':False,'font.size':10})
DOORS={'both_open':'兩門全開','both_closed':'兩門全關','left_open':'左開右關','right_open':'右開左關'}
def name(l,d):return f'{l}_{d}_T33_Q30'
def label(c):return ('靠左' if c.startswith('corner') else '置中')+'・'+next((v for k,v in DOORS.items() if k in c),c)
def summary(c):
 p=O/c;d=json.loads((p/'summary.json').read_text())
 if (p/'complete_event_summary.json').exists():d.update(json.loads((p/'complete_event_summary.json').read_text()))
 return d
def table(d):
 def cell(x):
  if isinstance(x,(float,np.floating)):return f'{x:.4g}'
  return str(x).replace('|','/').replace('\n',' ')
 return '| '+' | '.join(map(str,d.columns))+' |\n|'+'|'.join(['---']*len(d.columns))+'|\n'+'\n'.join('| '+' | '.join(cell(x) for x in row)+' |' for row in d.itertuples(index=False,name=None))

def save(fig,n):
 fig.text(.5,.013,'天花33°C・實際流量30m³/min；沉積閉合未校準，不能用微小模型值證明不會積油',ha='center',color='#a22',fontsize=8)
 fig.savefig(F/(n+'.png'),dpi=180);fig.savefig(F/(n+'.pdf'));plt.close(fig)
def data_tables(cases):
 ss=pd.DataFrame([summary(c) for c in cases]);ss.to_csv(O/'scenario_summary.csv',index=False)
 rows=[];sizes=[]
 allcases=cases+[name(l,'both_closed')+'_stick'+str(a) for l in ['corner','center'] for a in [.1,.5]]
 for c in allcases:
  manifest=json.loads((ROOT/'runs'/c/'replay_manifest.json').read_text());mapping=json.loads((Path(manifest['source'])/'boundary_update.json').read_text());idmap=dict(zip(mapping['new_ceiling_face_ids'],mapping['old_ceiling_face_ids']))
  for csvname in ['ceiling_faces.csv','ceiling_complete_event_faces.csv']:
   fp=O/c/csvname;dd=pd.read_csv(fp);dd['v1_face_id']=dd.original_face_id.map(idmap);assert dd.v1_face_id.notna().all();dd.to_csv(fp,index=False)
  for fp in (O/c).glob('*.vtp'):
   tree=ET.parse(fp);cd=tree.find('.//CellData');native=next(x for x in cd if x.attrib.get('Name')=='original_face_id');ids=np.fromstring(native.text,sep=' ',dtype=int)
   existing=next((x for x in cd if x.attrib.get('Name')=='v1_face_id'),None)
   if existing is None:existing=ET.SubElement(cd,'DataArray',type='Int64',Name='v1_face_id',format='ascii')
   existing.text=' '.join(str(idmap[int(i)]) for i in ids);tree.write(fp,encoding='utf-8',xml_declaration=True)
  s=summary(c);f=pd.read_csv(O/c/'ceiling_complete_event_faces.csv');alpha=float(c.split('_stick')[1]) if '_stick' in c else 1
  for gsd in [1.2,1.5,2]:
   unit=f[[f'M_{d:g}um_unit_mg_m2_per_mg' for d in DS]].to_numpy()@size_weights(gsd)
   sizes.append(dict(case=c,gsd=gsd,ceiling_fraction=float(unit@f.area_m2)))
  for material in ['gypsum_matte','gypsum_satin','coated_aluminium']:
   for days in [30,180,365]:
    for events in [1,2]:
     for interval,eff in [(0,0),(30,.5),(30,.9),(30,.99)]:
      retained,peak,before,removed=cycle_units(days,events,interval,eff)
      for dose in [1,4.4,7.2,8.8,28.8]:
       rows.append(dict(case=c,material=material,sticking_probability=alpha,days=days,events_per_day=events,emitted_mg_per_event=dose,clean_every_days=interval,clean_efficiency=eff,retained_total_mg=s['M_total_unit_mg_per_mg']*dose*retained,retained_mean_mg_m2=s['M_mean_unit_mg_m2_per_mg']*dose*retained,retained_p95_mg_m2=s['M_p95_unit_mg_m2_per_mg']*dose*retained,retained_max_mg_m2=s['M_max_unit_mg_m2_per_mg']*dose*retained,peak_mean_mg_m2=s['M_mean_unit_mg_m2_per_mg']*dose*peak,removed_total_mg=s['M_total_unit_mg_per_mg']*dose*removed))
 pd.DataFrame(rows).to_csv(O/'long_term_scenarios.csv',index=False);pd.DataFrame(sizes).to_csv(O/'size_distribution_sensitivity.csv',index=False)
 return ss,allcases

def validation(cases):
 rows=[]
 base=summary(name('corner','both_closed'))
 for suffix,scope in [('_dt005','full-history dt 0.1 to 0.05 s'),('_frames4','full-history flow-frame 2 to 4 s')]:
  p=O/(name('corner','both_closed')+suffix)/'summary.json'
  if p.exists():
   v=json.loads(p.read_text());# Compare observed deposition only: clearance cancels neither history nor residual.
   b=json.loads((O/name('corner','both_closed')/'summary.json').read_text())
   rows.append(dict(test=suffix,scope=scope,base_ceiling_fraction=b['M_total_unit_mg_per_mg'],test_ceiling_fraction=v['M_total_unit_mg_per_mg'],relative_change_percent=100*(v['M_total_unit_mg_per_mg']/b['M_total_unit_mg_per_mg']-1)))
 a=O/(name('corner','both_closed')+'_startup20')/'summary.json';b=O/(name('corner','both_closed')+'_fine_startup20')/'summary.json'
 if a.exists() and b.exists():
  x=json.loads(a.read_text());y=json.loads(b.read_text());rows.append(dict(test='startup20_mesh_and_initialization',scope='0–20 s coarse/fine startup diagnostic; fine uses one-time initial flux projection and differs in startup stepping. NOT an isolated mesh-convergence test',base_ceiling_fraction=x['M_total_unit_mg_per_mg'],test_ceiling_fraction=y['M_total_unit_mg_per_mg'],relative_change_percent=100*(y['M_total_unit_mg_per_mg']/x['M_total_unit_mg_per_mg']-1)))
 pd.DataFrame(rows).to_csv(O/'sensitivity.csv',index=False)
 worst=max(summary(c)['max_relative_mass_error'] for c in cases)
 v=dict(max_aerosol_relative_mass_error=worst,particle_mass_conservation_pass=worst<1e-7,unsteady_full_event_mesh_convergence_demonstrated=False,original_strict_mesh_failures_retained=True,empirical_wall_deposition_calibration=False,gas_phase_condensation_included=False,sensitivity=rows,interpretation='Functionally working, conservative unit-response screening; absolute deposition and material ranking NOT validated.')
 (S/'VALIDATION.json').write_text(json.dumps(v,ensure_ascii=False,indent=2));return v

def plots(cases,ss):
 for d,zh in DOORS.items():
  pair=[name(l,d) for l in ['corner','center']];faces=[pd.read_csv(O/c/'ceiling_complete_event_faces.csv') for c in pair];hi=max(float(f.M_unit_mg_m2_per_mg.max()) for f in faces);hi=max(hi,1e-30);norm=LogNorm(hi*1e-8,hi)
  fig,axs=plt.subplots(2,1,figsize=(9,9));fig.subplots_adjust(left=.1,right=.81,bottom=.09,top=.86,hspace=.5)
  for ax,c,f in zip(axs,pair,faces):
   g=mesh_geometry(ROOT/'runs'/c);ids=f.original_face_id.to_numpy(int);pc=PolyCollection([g['points'][g['faces'][i],:2] for i in ids],array=np.maximum(f.M_unit_mg_m2_per_mg,norm.vmin),norm=norm,cmap='inferno');ax.add_collection(pc);ax.set(xlim=(0,3.3),ylim=(0,1.65),xlabel='x（m）',ylabel='距後牆 y（m）',title=label(c));ax.set_aspect('equal');ax.set_facecolor('#ddd');s=summary(c)
   ax.text(.01,.99,f"總量 {s['M_total_unit_mg_per_mg']:.3e} mg／每排放1mg\n7.2mg參考事件：{s['M_total_reference_mg']:.3e} mg",transform=ax.transAxes,va='top',fontsize=8,bbox={'facecolor':'white','alpha':.85,'edgecolor':'none'})
   bx=[.95,1.35] if c.startswith('corner') else [1.8,2.2];ax.scatter(bx,[.3,.3],c='lime',marker='x');j=f.M_unit_mg_m2_per_mg.idxmax()
   if s['M_total_reference_mg']>=1e-12:ax.scatter([f.loc[j,'x_m']],[f.loc[j,'y_m']],marker='*',s=60,c='deepskyblue',edgecolors='white')
   else:ax.text(.5,.1,'數值極小，不判讀實際污漬位置',transform=ax.transAxes,ha='center',color='#a22',bbox={'facecolor':'white','alpha':.9,'edgecolor':'none'})
  fig.colorbar(pc,ax=axs,shrink=.75,label='mg/m²／排放1mg（對數；色階下限非檢出限）');fig.suptitle(f'{zh}・靠左與置中\n33°C天花，30m³/min；0–240s＋固定末端流場清空積分',fontsize=14);save(fig,'ceiling_'+d+'_1to1')
 fig,ax=plt.subplots(figsize=(9,9));fig.subplots_adjust(left=.23,bottom=.12);ax.barh([label(c) for c in cases],ss.M_total_unit_mg_per_mg.clip(lower=1e-300));ax.set_xscale('log');ax.set(xlabel='天花沉積mg／每排放1mg',title='八工況天花沉積・條件式模型響應');save(fig,'total_1to1')
 fig,ax=plt.subplots(figsize=(9,9));fig.subplots_adjust(left=.23,bottom=.12);left=np.zeros(len(cases))
 for key,zh in [('hood_fraction','抽走'),('left_fraction','左門'),('right_fraction','右門'),('other_wall_fraction','其他表面'),('ceiling_fraction','天花'),('airborne_fraction','清空後剩餘')]:
  val=ss[key].to_numpy();ax.barh([label(c) for c in cases],val,left=left,label=zh);left+=val
 ax.set_xlim(0,1.01);ax.legend(loc='upper center',ncol=3);ax.set(title='已形成油滴的質量平衡',xlabel='排放質量比例');save(fig,'mass_balance_1to1')
 with PdfPages(F/'all_figures_1to1.pdf') as pdf:
  for f in F.glob('*_1to1.png'):
   fig=plt.figure(figsize=(9,9));plt.imshow(plt.imread(f));plt.axis('off');plt.tight_layout(pad=0);pdf.savefig(fig);plt.close(fig)

def report(ss,v):
 flow=pd.read_csv(S/'flow_validation.csv');ts=ss[['case','M_total_unit_mg_per_mg','M_total_reference_mg','ceiling_percent','M_mean_reference_mg_m2','M_p95_reference_mg_m2','max_x_m','max_y_m']].copy();ts.loc[ts.M_total_reference_mg<1e-12,['max_x_m','max_y_m']]=np.nan;ts['case']=ts['case'].map(label);ts.columns=['配置','沉積mg／排放1mg','7.2mg事件沉積mg','排放佔比%','平均mg/m²','P95 mg/m²','最大面x','最大面y']
 fl=flow[['case','Q_late_m3_min','hood_static_suction_late_Pa','exceeds_fan_max_1300Pa']].copy();fl['case']=fl['case'].map(label);fl.columns=['配置','實際m³/min','模型煙機入口吸力Pa','超過1300Pa']
 sens=pd.DataFrame(v['sensitivity'])
 text=f'''# 天花33°C、煙機30m³/min：油滴沉積更新

2026-10-21。此報告為新情境，舊版模型、結果及圖表保留。可重現參數見parameters.json，材料資料見materials.json，研究依據見SOURCES_zh.md。

**功能與質量守恆已檢查；沉積定量仍未完成物理驗證。不能以很小的模型沉積值保證石膏板不會積油。** 本輪沒有氣態有機物冷凝、油膜老化、清潔難度或黃變預測。

## 新邊界與實際計算

兩個爐頭各對空氣放熱2000W，共4000W；初始空氣及補風25°C。天花表面由0秒起固定33°C，在能量方程中實際與空氣交換熱量，不是把舊結果改標籤。原網格的外邊界面重新歸組以新增ceiling，沒有改變幾何、點座標、單元拓撲或內部面；新面ID有重排，對照保存在boundary_update.json；輸出表格的original_face_id是新情境內的網格面ID，另加v1_face_id指向舊網格的面ID，不能把兩欄混用。

煙機強制排出0.5m³/s。這是指定實際流量的情境，不是由一條實測風機曲線算出的工作點；出口採flowRateOutletVelocity，與舊均壓抽吸邊界的速度分布並不相同。下表吸力為房間至煙機入口的模型靜壓需求，已用25°C環境的靜水壓分布修正煙機所在高度；沒有加外部風管、濾網等損失，也沒有估計風機總壓。

{table(fl)}

**全關門僅有兩道5mm底縫補風，若所需壓差超過1300Pa，30m³/min工況不能當作該煙機實際可達的工況。** 不能同時把最大自由風量和最大靜壓當成同一工作點；實際工況需要風機曲線與管路阻力。開門方案仍需考慮未建模的外部壓力／橫風。

熱流場用kitchenCeilingAirFoam耦合求解浮力、SST湍流和能量，原生面質量通量每2秒儲存（含暫停檢查點）；開門長算例maxCo=6、最大步長0.1秒，全關門maxCo=12、最大步長0.04秒，實際由Courant限制；早期開門pilot沿用原maxCo=3、最大步長0.05秒。尚未證明完整事件流場時間步獨立。另用ceilingReplayFoam重解五種粒徑的守恆油滴輸運。油滴負荷很低，沒有反作用於空氣動量／能量，是這個稀薄氣膠模型的假設。儲存流場後不再重複計算已成功的流場。新版流場中的Yp欄位只是相容性占位，不是油滴計算成果；真正結果在runs及results對應案例。

## 排放與三種材料

採中式住家煮食文獻中的一次有機氣膠作油污代理：1.8mg/min／頭示例，雙爐3.6mg/min、120秒共7.2mg。文獻不同菜式1.1–2.2mg/min，以及舊商用CFD參數7.2mg/min／頭，另列敏感度情境。把每道菜的實驗值套到每個爐頭並相加是工程外推；不是純油脂量測，不是概率上下界。每排放1mg的結果較適合轉用；整頓飯若時長不同，熱流場也會變，不能只把本次120秒結果稱為完整一頓飯。

五種粒徑0.1/0.3/1/3/10µm，密度850kg/m³；以2.5µm質量中位徑、假設GSD1.5混合，另測1.2與2。中位徑參考體積分布，在等密度假設下轉為質量分布。不是完整超細粒子數目分布。

三類表面都固定33°C、採共同α=0.1/0.5/1敏感度。附著率掃描只在兩個全關門方案執行，其餘六個門工況為α=1，不能外推成所有工況都已完成材料敏感度。未找到可直接代表其微米油滴附著率的資料，不會編造「啞光1、Satin0.5、鋁板0.1」這種材料排序。外層漆膜而非基材接觸油滴；防潮石膏板導熱0.19W/(m·K)、漆膜乾厚約36–43µm等產品資料已記錄，但固定表面溫度時不會由材料導熱／熱容量產生額外差異。所有表面暫按水力光滑處理，並非已量測粗糙度。

Satin廚衛產品的防霉、可清洗說明與金屬板中性清潔劑維護指引，是獨立產品資料，不能當成較低沉積率的證據。耐擦洗應看具體產品測試等級；光澤名稱不足以判斷。普通啞光參考產品也允許固化後清洗。清潔效率50/90/99%只是假設，不是三種材料的實測差異。

## 天花總量、位置與表面負荷

0–120秒排放與加熱，120秒停止；持續排風到240秒，再用凍結末端熱流場的線性清空積分分配殘留粒子。這段代數清空不是指定未來時間，並非一直維持實際室溫演變。總量、平均與P95按面積統計；最大網格面的位置只對此模型有效，近數值零的工況不判讀實際污漬位置。圖表中以參考總量1e−12mg作純粹的顯示提醒，並非量測檢出限或可接受油污標準；原始CSV仍保留完整數值。相對色階會放大極小值的差異，不能把亮色直接解讀成實際嚴重積油。

{table(ts)}

M_dep是逐面沉積通量J_dep的時間積分；天花總量為ΣM_depA。CSV輸出mg/m²與kg/(m²·s)，ParaView的PVD只含0–240秒真實時間，另用ceiling_complete_event.vtp儲存清空後靜態結果。沒有把天花附近空氣濃度冒充沉積。

## 驗證與可相信程度

最大瞬態粒子相對質量殘差：{v['max_aerosol_relative_mass_error']:.3e}。這只驗證粒子收支，不驗證真實沉積量。存檔流場插值另有局部空氣連續性殘差，必須與粒子收支區分；開始時沒有原生phi的0秒狀態暫用靜止零通量，會影響最初一段插值，不能稱為所有瞬間出口都精確0.5m³/s。各案例summary.json另報此空氣殘差及Q。原網格嚴格檢查的失敗項目及近壁不足仍存在；舊工況曾有約27倍粗／細分支差異，新設定不能自動消除問題。

{table(sens) if len(sens) else '沒有可比較的敏感度結果。'}

短窗粗／細診斷的初始存檔都是25°C、相同33°C壁面、相同流量及零油滴，但細網格保留了原專案的一次初始面通量投影，並採不同的啟動步長；投影會調整初始速度。因此這是網格與初始化混合的啟動診斷，不能把差異全歸因網格，更不能證明整個120秒加熱及清空歷程網格獨立。全事件時間步與存檔間距測試使用同一流場，與網格測試分開。壁面模型含布朗／湍流擴散阻力、重力、Talbot熱泳、低St平均慣性和Wood型湍流慣性接觸；管道經驗式外推至室內天花未校準，高St細縫附近亦超出低St漂移近似的理想範圍。

## 長期累積及實際決策

CSV含每日1/2次，30/180/365日，不清潔及每30日清潔50/90/99%；保留量及清潔前峰值均計算。主要每次排放1mg，絕對例為4.4/7.2/8.8/28.8mg。假設每次事件相同、期間完全清空、附著物除清潔外不流失；沒有老化、再懸浮、滴落、冷凝或油膜飽和。

1. 最容易積油的位置：表中最大面是條件式模型熱點，需配合俯視圖；不能把未收斂模型的最大面當作施工分區的精確位置。
2. 每次有多少：表中7.2mg參考事件與每排放1mg響應均已輸出；物理驗證不足，不能直接量化住家真實油脂量。
3. 最關鍵假設：近壁輸運／湍流慣性閉合與網格、實際粒徑質量分布、補風橫流與可達抽風量、附著率，還有完全未模擬的氣態冷凝。
4. 提高風量改善多少：這輪指定30m³/min，沒有同時做33°C下完整事件的另一風量，所以不能提供新的改善百分比。舊25°C邊界短窗增風比較不會拿來充當本輪證據；全關門更需先解決補風和可達工作點。
5. 石膏板是否值得關注：目前不能用很小的模型值排除長期清潔問題。防潮不代表防油；清潔維護以外層漆膜／接縫和具體產品說明為準。選材時，耐擦洗漆或可依廠商方式清潔的完整烤漆板，有較直接的維護資料；本CFD不能斷言哪種較少積油。
6. 可信度及實測需要：流程、單位和質量守恆可查，沉積量與材料排序仍屬探索性。量化至少需實際罩口風量／房間壓差、油滴質量粒徑分布及排放率、不同區域可拆卸試片的單位面積增重或化學油脂量、表面溫度與清潔前後移除量。這些用於日後校準；本輪不要求你先提供。

## 來源及檔案

研究與官方產品來源完整連結見[SOURCES_zh.md](SOURCES_zh.md)。原專案的沉積方程和限制見../../reports/02_physics_and_sources_zh.md。主要資料在results/scenario_summary.csv、long_term_scenarios.csv、size_distribution_sensitivity.csv、各案例逐面CSV與ParaView檔。互動3D預設每次排放1mg，支援粒徑分布與清潔情境。
'''
 (S/'analysis_zh.md').write_text(text)
 html=MarkdownIt('commonmark',{'html':True}).enable('table').render(text)
 for n in ['ceiling_both_open_1to1','ceiling_both_closed_1to1','ceiling_left_open_1to1','ceiling_right_open_1to1','total_1to1','mass_balance_1to1']:
  html+=f'<div style="page-break-before:always"><img style="width:100%" src="data:image/png;base64,{base64.b64encode((F/(n+".png")).read_bytes()).decode()}"></div>'
 html='<meta charset="utf-8"><style>@page{size:A4;margin:14mm}body{font:10px "Noto Sans CJK TC",sans-serif;line-height:1.6}h1{font-size:23px}h2{font-size:16px}table{border-collapse:collapse;font-size:8px;width:100%}td,th{border:1px solid #bbb;padding:3px;overflow-wrap:anywhere}img{max-width:100%}a{color:#167}</style>'+html
 (S/'analysis_zh.html').write_text(html)
 with sync_playwright() as p:
  b=p.chromium.launch(executable_path='/usr/bin/chromium',args=['--no-sandbox']);pg=b.new_page();pg.set_content(html,wait_until='load');pg.pdf(path=str(S/'analysis_zh.pdf'),print_background=True);b.close()

def main():
 cases=[name(l,d) for d in DOORS for l in ['corner','center']]
 from run_33C_Q30 import run_flow
 flows=[run_flow(ROOT/'flow_sensitivity'/c) for c in cases]
 pd.DataFrame(flows).to_csv(S/'flow_validation.csv',index=False)
 ss,allcases=data_tables(cases);v=validation(cases);plots(cases,ss)
 import make_interactive as mi
 mi.OUT=O;mi.main(output_path=S/'interactive_3d.html',case_names=allcases,scenario_text=dict(warning='沉積定量未校準：短窗網格測試不代表全事件收斂，不能作材料壽命設計值。',intro='新耦合流場：天花33°C、出口實際流量30m³/min；兩爐總放熱4000W。全關門工況須核對1300Pa上限。',reference='參考示例為兩爐合計3.6mg/min、120秒共7.2mg一次有機氣膠；由中式住家文獻外推，非純油脂實測。'))
 p=S/'interactive_3d.html';p.write_text(p.read_text().replace('原面ID','新情境面ID'))
 report(ss,v)
 (S/'PROGRESS.json').write_text(json.dumps(dict(status='COMPLETE_PHYSICAL_VALIDATION_LIMITED',completed=8,validation=v),ensure_ascii=False,indent=2))
 # A compact review archive; do not overwrite earlier archived results.
 dest=Path('/workspace/ceiling33_Q30_deliverables');dest.mkdir(exist_ok=True)
 with zipfile.ZipFile(dest/'ceiling33_Q30_review.zip','w',zipfile.ZIP_DEFLATED) as z:
  for p in S.rglob('*'):
   if p.is_file():z.write(p,p.relative_to(S))
  for n in ['prepare_33C_Q30.py','run_33C_Q30.py','update_33C_Q30.py','report_33C_Q30.py']:z.write(ROOT/'scripts'/n,'scripts/'+n)
 print('New scenario report, CSV, ParaView, figures and 3D generated',flush=True)
if __name__=='__main__':main()
