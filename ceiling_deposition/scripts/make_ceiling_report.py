from pathlib import Path
import json,base64,re,numpy as np,pandas as pd
from markdown_it import MarkdownIt
from playwright.sync_api import sync_playwright
from process_ceiling import ROOT,OUT
from visualize_ceiling import label,main_cases

def fmt(x):
 if pd.isna(x):return '—'
 return f'{x:.3g}'
def table(df):
 cols=list(df.columns);return '| '+' | '.join(cols)+' |\n|'+ '|'.join(['---']*len(cols))+'|\n'+'\n'.join('| '+' | '.join(str(v) for v in row)+' |' for row in df.itertuples(index=False,name=None))+'\n'
def main():
 p=ROOT/'reports';p.mkdir(exist_ok=True)
 primary=pd.read_csv(OUT/'primary_complete_event_summary.csv');sens=pd.read_csv(OUT/'numerical_and_physical_sensitivity.csv');fan=pd.read_csv(OUT/'fan_intervention_comparison.csv');v=json.loads((ROOT/'verification/final_validation.json').read_text());long=pd.read_csv(OUT/'long_term_scenarios.csv')
 tab=[];details=[]
 for _,r in primary.iterrows():
  tab.append({'配置':label(r['case']),'天花mg／排放1mg':fmt(r.M_total_unit_mg_per_mg),'參考總量mg':fmt(r.M_total_reference_mg),'排放佔比%':fmt(r.ceiling_percent),'最大面x,y(m)':f'{r.max_x_m:.2f}, {r.max_y_m:.2f}' if 'both_closed' in r['case'] else '近零，不判讀熱點'})
  if 'both_closed' in r['case']:details.append({'配置':label(r['case']),'平均mg/m²':fmt(r.M_mean_reference_mg_m2),'P95 mg/m²':fmt(r.M_p95_reference_mg_m2),'最大mg/m²':fmt(r.M_max_reference_mg_m2),'240s後份額%':fmt(r.tail_deposition_share*100)})
 stats=table(pd.DataFrame(tab));detail=table(pd.DataFrame(details))
 fanrows=[]
 for _,r in fan.iterrows():fanrows.append({'位置':'靠左' if r.layout=='corner' else '置中','實際Q增幅%':fmt(r.Q_increase_percent),'對照mg':fmt(r.deposition_baseline_reference_mg),'增風mg':fmt(r.deposition_boost_reference_mg),'天花沉積減幅%':fmt(r.deposition_reduction_percent)})
 ss=sens[['test','scope','relative_change_percent']].copy();ss.columns=['測試','適用時段／範圍','沉積變化%'];ss['沉積變化%']=ss['沉積變化%'].map(fmt)
 l=long[(long.material=='gypsum_matte')&(long.sticking_probability==1)&(long.emitted_mg_per_event==28.8)&(long.clean_efficiency.isin([0,.9]))&(long.case.str.contains('both_closed'))]
 ll=[]
 for _,r in l.iterrows():ll.append({'位置':'靠左' if r['case'].startswith('corner') else '置中','日數':int(r.days),'每日次數':int(r.events_per_day),'清潔':'無' if r.clean_every_days==0 else '30日/90%','保留總量mg':fmt(r.retained_total_mg),'平均mg/m²':fmt(r.retained_mean_mg_m2),'清潔前峰值平均mg/m²':fmt(r.peak_before_clean_mean_mg_m2)})
 worst=primary.loc[primary.M_total_unit_mg_per_mg.idxmax()]
 md=f'''# 廚房天花油脂沉積分析

計算及整理：2026-10-21。原始基線為 v1.0.0；保留原幾何、網格、煙機設定和成果，另建立獨立的 OpenFOAM 沉積擴展。

**驗證結論：天花沉積定量未通過穩健性檢查。後段粗／細網格分支相差約27倍，10→20秒存檔間距改變約38%。細網格分支的流場歷程亦有差別，不能把全部差異歸因單一網格因素；但已足以否定目前量值可作可靠設計值。**

**本次已實作並運行真正的壁面質量通量與逐面累積，並非天花附近濃度圖。但絕對沉積量未經實測校準，不能把很小的模型值解讀成「石膏板不會積油」。本輪沒有氣態油脂冷凝或漆面黃變模型。**

## 1. 原專案審查與選擇

OpenCFD OpenFOAM 2406、自訂 buoyantPimpleFoam 衍生求解器、瞬態 PIMPLE。已有理想氣體密度／溫度／浮力／SST 湍流，兩爐各向空氣放熱2000 W。五個 Eulerian 油滴粒徑0.1/0.3/1/3/10 µm，含漂移擴散及總壁面 sink；不是純 passive scalar，亦沒有 Lagrangian parcels。完整審查在 `01_project_audit_zh.md`。

原邊界為 hood、left、right、cabinetSide、walls；天花原在 walls 內。以原面ID、法向和 z=2.2 m 分出 ceiling faceSet/faceZone，**不重排原面、重建幾何或換網格**。每個佈局的暴露天花面積由原網格量得；高櫃與上櫃覆蓋區沒有對應空氣接觸面。

採 Eulerian 延伸，可重用原生保守面通量 phi 和保存的 U/T/p/k/nut。另解油滴傳輸，避免直接以舊濃度乘一個新沉積速度而破壞收支。Lagrangian 對大滴軌跡有優點，但細滴需近壁隨機湍流、布朗力、熱泳、慣性與附著模型及粒子數收斂；單純 trap 接觸不是完整細粒子沉積模型。

## 2. 模型與時間定義

每個粒徑獨立解保守 Eulerian 輸運。壁面以布朗／van-Driest 湍流阻力、正負法向重力漂移、Talbot 熱泳、低 Stokes 數平均流慣性滑移，以及 Wood 型湍流慣性接觸閉合計算移除質量。主篩查採慣性接觸另計的 parallel 閉合；另保留等效阻力閉合作對照。管道經驗式外推至天花是重要未校準假設，不能視作精確碰撞軌跡。附著率只在粒子到達表面的輸運計算後作用，沒有把「撞牆全黏」冒充全部輸運機理。

逐面 J_dep=ρYv_dep，單位 kg/(m²·s)；M_dep=ΣJ_dep Δt，原生單位 kg/m²，匯出 mg/m²時乘10⁶。**方程中隱式移除的質量與逐面累積是同一份質量**。全體天花總量為ΣM_dep A，P95按面積加權。完整方程、資料來源與限制見 `02_physics_and_sources_zh.md`。

沿用120秒雙爐排放與加熱，120秒停止，排風觀測至240秒。240秒後尚有懸浮油滴，另在固定末端熱流場下解線性殘留清空積分（resolvent）；排出、天花和其他表面仍逐項守恆。代數正則尺度10¹²秒**不是物理煮食或觀測時間**，不製作該時間的物理動畫。主結果包含這段清空近似，並分開保存240秒觀測量和長尾份額；固定末端溫度／流場仍是模型限制。

## 3. 八工況天花結果

優先看每排放1mg的響應。絕對參考值使用兩爐合計14.4 mg/min、持續2分鐘，即28.8 mg已形成油滴；沿用質量中位徑2.5 µm、假設GSD1.5。它是商用文獻外推，**不是這個住家的排放實測、PM2.5或全部氣態油脂**。排放量另列0.25–4倍情境，非住家量測信賴區間。

{stats}

全關門的面積統計（同一28.8mg參考排放）：

{detail}

**開門案例的天花沉積小到不足以作工程量化；不按其指數、倍率或最大面位置判斷材料風險。顏色即使有明顯相對差異，也不代表可見積油。** 在主篩查假設下，較大天花沉積出現在全關門；最大片區位置以CSV與下圖為準。這不表示全室暴露與天花沉積總會同方向變化。

![全關門天花俯視](../figures/01_ceiling_both_closed_1to1.png)

綠×為爐頭投影，藍星為模型最大網格面，灰色為櫃體佔用區；天花最大點不必正對爐頭，取決於到達天花的油滴供應、近壁氣流、溫差和沉積速度。單一網格面極值不等於已定位的實際污漬。

![總沉積](../figures/02_total_deposition_1to1.png)

## 4. 抽走、懸浮與其他表面

逐步輸出抽油煙機、左右門淨流出、天花、其他表面及懸浮，開口採求解矩陣的實際數值通量。殘留清空後仍剩的代數殘量沒有被強行配給天花。完整去向見 `primary_complete_event_summary.csv`、各案例 `size_mass_balance.csv` 和原生 `ceiling_budget.csv`。

![質量去向](../figures/03_mass_balance_1to1.png)

## 5. 提高抽風量的實際試驗

沒有縮放舊U冒充新CFD：從約100秒同一檢查點，重新解原風量對照及目標增20%的熱流，保留1300 Pa壓力上限，求解至140秒，再以相同舊油滴初始場比較沉積。這是**相同重啟的短時控制試驗，不是整次煮食從0秒提高風量**；初始油滴來自舊模型，因此此處的mg不可與第3節的新模型完整事件直接相加或比較。

{table(pd.DataFrame(fanrows))}

「沉積減幅」負值表示這個短時試驗增加。當絕對量極小，龐大百分比不應當作有實用精度的風險差異。較強排風會改變流線、油滴在天花附近的供應和慣性接觸；本輪不支持宣稱風量增加20%便按比例減少天花積油。要量化整頓飯的改善，需量測裝機風量／補風，再做全歷程不同工作點的耦合試驗。

## 6. 30／180／365日累積與清潔

配置檔 `config/long_term.json` 提供每日1／2次、每次1mg單位排放、7.2／28.8／115.2mg三個絕對情境、無清潔與每30日清潔情景。清潔在當日烹調後進行，分別輸出清潔後保留量和清潔前峰值，效率0.5／0.9／0.99只是可改假設。

以相同短事件反覆累積，不含自然掉落、再懸浮、氧化、油膜改變附著率或氣態冷凝。**不是幾個月會發黃、幾微米油膜或漆膜壽命預測。** 如果一頓飯的排放歷程與此120秒案例不同，應更換每次排放劑量並檢查新的流場歷程，不能無條件按烹調分鐘數外推。

![單位排放的長期累積](../figures/04_long_term_normalized_1to1.png)

以下僅為28.8mg/次、附著率1、全關門的條件式量值。所有材料在同一附著率與清潔效率下數值相同：

{table(pd.DataFrame(ll))}

## 7. 材料比較

| 選項 | CFD能比較 | 應另查的實際資料 |
|---|---|---|
| 防潮石膏板＋普通啞光乳膠漆 | 同一α範圍下的油滴附著量 | 漆膜耐擦洗、耐油污、清潔劑相容性；防潮板不等於防油或可反覆濕擦 |
| 防潮石膏板＋耐擦洗Satin廚衛漆 | 同上；不憑光澤指定較低α | 實際產品耐擦洗／耐污資料與施工、通風要求；亦有耐油污的啞光配方 |
| 烤漆鋁板 | 同上；塗層未校準α | 塗層和清潔劑相容、接縫、可拆洗性；不能假定鋁板表面的塗層永不老化 |

附著率α=0.1／0.5／1在靠左與置中兩個全關門案例實際重新求解；開門案例只保留α=1，因其量值已不足以作材料定量比較。三種材料共用這個掃描範圍，沒有足夠資料為任何一種排黏附率名次。CFD相同α的結果重合是正確呈現未知，不是漏做材料差異。可擦洗性和方便清潔不等於低沉積率。基於維護規格，應優先選擇有實際耐擦洗／耐油污文件且可方便維護的系統；此建議來自材料資料，並非CFD預測塗料老化。

一手材料資料：[Gyproc防潮板](https://www.gyproc.ie/documents/product-datasheets/pds-gyproc-ireland-moisture-resistant.pdf)、[Benjamin Moore廚衛Satin TDS](https://media.benjaminmoore.com/WebServices/prod/assets/production/datasheets/TDS_0322/W322_TDS_US.pdf)、[Dulux Kitchen Matt](https://www.dulux.co.uk/en/products/dulux-easycare-kitchen-matt/cccId_1935922)、[金屬吊頂清潔維護](https://assets.hunterdouglasarchitectural.eu/documents/ceiling/Installation%20Manuals/_default/Cleaning-and-Maintenance-Luxalon-multipanel.pdf)。它們都不提供本廚房的油滴附著率量測。

## 8. 實際驗證結果與可信程度

低成本20秒流程先通過才進入正式批次。量綱、零漂移／順逆重力極限、混合箱收支、非負濃度、逐面積分、VTK讀取和每步油滴質量收支已核對。完成瞬態案例{v['completed_transient_runs']}個，清空積分{v['completed_clearance_integrals']}個。最差瞬態相對收支誤差{fmt(v['worst_transient_relative_mass_error'])}，清空積分{fmt(v['worst_clearance_relative_mass_error'])}。這是**離散收支驗證，不是物理實驗驗證**。

{table(ss)}

時間步比較為完整0–240秒；參數與加密比較為110–140秒的匹配短排放試驗。細網格使用原專案已有約3萬單元案例；它們的氣流也由舊敏感性分支取得，沒有完成0–240秒全歷程網格獨立。嚴格checkMesh仍有2項未通過；主粗網格的天花剪力估算y+最大約85，首層法向距離最大約5.6 cm，並非y+≈1完整解析，故沉積依賴近壁閉合。最大10µm的低Stokes滑移近似在門縫高速區會受限。

原熱流的雙關門熱量重建殘差較大；重播原始存檔不能修復原熱流誤差，10秒存檔間距也丟失未保存的瞬態，尤其首個0–10秒不能還原快速的風機啟動與熱羽流建立。風機工作點仍由網路假設給定，未聯立實測完整風機／煙管／濾網曲線。**絕對積油與材料壽命判斷可信度低；程式可運行和質量守恆的可信度較高。若加密或閉合使極小沉積大幅改變，就不能將小數字當作可靠上限。**

![粒徑響應](../figures/05_size_response_1to1.png)

![閉合敏感性](../figures/06_wall_model_uncertainty_1to1.png)

## 9. 最後六個問題

1. **哪裡最容易積油？** 在本次已形成油滴的主情境中，全關門較值得注意；最大面位置見第3節座標及熱圖。置中／靠左的熱點取決於近壁輸運，不只看爐頭正上方。開門的近零值不足以定位真實污漬。
2. **單次多少？** 第3節列出每排放1mg及28.8mg參考排放的條件量，已含末端流場殘留清空近似。它不是實測一頓飯，更不包含氣態物質冷凝。
3. **哪些假設最影響？** 實際排放劑量與粒徑分布、門／補風和裝機風量、天花溫差、近壁模型／網格、慣性與附著率，及未建模的氣態凝結。敏感度數字見第8節，不能用質量守恆誤差代替這些不確定性。
4. **增風改善多少？** 約100–140秒匹配試驗的實際Q與沉積變化見第5節；不保證單調改善，不能當整頓飯的改善百分比。
5. **石膏板吊頂會否有長期清潔問題？** 本輪無法排除。極低的油滴沉積不包含蒸氣／半揮發性有機物的表面吸附凝結，也未預測漆面變色。應按頻繁炒菜的維護需求選耐擦洗、耐油污且可維修的表面，不能由此宣稱普通啞光石膏吊頂長期免清潔。
6. **還需量甚麼？** 裝機風量／煙管壓差、門口補風速度；每次炒菜的油滴總質量和質量粒徑分布；天花溫度隨時間；在圖示熱點與低區放置已稱重的材料試片，按多次相同烹調測沉積及清潔後殘留；若要納入凝結，另量氣相／粒相有機物組成和分配。以同場三種塗層試片校準α及清潔效率，才可量化長期積油差異。

## 10. 交付索引

README列出完整執行步驟。`results/` 含逐面CSV、原生ID、PVD/VTP、五粒徑及混合質量收支、長期情景與敏感度表；`figures/` 含1:1 PNG/PDF圖集；`ceiling_deposition_interactive_3d.html` 可完全離線操作。`runs/` 和 `flow_sensitivity/` 保存新案例，不覆蓋 v1.0.0 原始成果。物理模型和一手文獻完整列表在 `02_physics_and_sources_zh.md`。
'''
 quality=ROOT/'verification/coupled_interventions_quality.csv'
 if quality.exists():
  q=pd.read_csv(quality)[['case','Q_late_m3_s','reconstructed_energy_residual_RMS_W']].copy()
  q.columns=['新氣流分支','後段實際Q m³/s','重建能量收支RMS W']
  for col in q.columns[1:]:q[col]=q[col].map(fmt)
  md+='\n## 附錄A：新氣流分支的收支限制\n\n'+table(q)+'\n這是診斷式重建殘差，不等於經實測的能量誤差；原風量對照的熱量殘差仍不可忽略，未被油滴收支良好所修復。\n'
 md+='\n## 附錄B：完整物理模型與一手來源\n\n'+(p/'02_physics_and_sources_zh.md').read_text()
 (p/'ceiling_deposition_analysis_zh.md').write_text(md)
 body=MarkdownIt('commonmark',{'html':True}).enable('table').render(md)
 def embed(m):
  image=(p/m[1]).resolve();return 'src="data:image/png;base64,'+base64.b64encode(image.read_bytes()).decode()+'"' if image.exists() else m[0]
 body=re.sub(r'src="([^\"]+\.png)"',embed,body)
 css='''@page{size:A4;margin:17mm 14mm}body{font-family:"Noto Sans CJK TC",sans-serif;color:#16222d;font-size:10pt;line-height:1.55}h1{font-size:23pt;color:#143e58}h2{font-size:14pt;margin-top:18pt;break-after:avoid}table{border-collapse:collapse;width:100%;font-size:8pt;margin:10pt 0;table-layout:auto}th,td{border:1px solid #b7c5cf;padding:5px;word-break:break-word}th{background:#e5eef4}tr{break-inside:avoid}thead{display:table-header-group}img{display:block;width:92%;max-height:230mm;object-fit:contain;margin:8pt auto;break-inside:avoid}a{color:#166d95}strong{font-weight:700}code{font-size:8pt;overflow-wrap:anywhere}p{orphans:3;widows:3}'''
 page='<!doctype html><html lang="zh-Hant"><meta charset="utf-8"><style>'+css+'</style><body>'+body+'</body></html>'
 (p/'ceiling_deposition_analysis_zh.html').write_text(page)
 with sync_playwright() as pw:
  b=pw.chromium.launch(executable_path='/usr/bin/chromium',headless=True,args=['--no-sandbox']);tab=b.new_page();tab.set_content(page,wait_until='load');tab.pdf(path=str(p/'ceiling_deposition_analysis_zh.pdf'),print_background=True,prefer_css_page_size=True);b.close()
 print('Chinese PDF written',len(primary),'main cases')
if __name__=='__main__':main()
