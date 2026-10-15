"""Traditional Chinese report and numeric tables from saved simulation results."""
from pathlib import Path
import json,csv,subprocess,shutil
import numpy as np
ROOT=Path(__file__).resolve().parent
LAYOUTS=['corner','center'];CASES=['both_open','both_closed','left_open','right_open']
LABEL={'corner':'靠左角落','center':'檯面中央','both_open':'兩門全開','both_closed':'兩門全關','left_open':'左開右關','right_open':'右開左關'}

def table(head,rows):return '\n'.join(['| '+' | '.join(head)+' |','| '+' | '.join(['---']*len(head))+' |']+['| '+' | '.join(str(x) for x in row)+' |' for row in rows])
def read(layout,name):return json.loads((ROOT/layout/'results'/name).read_text())
def f(x,n=2):return '未達' if x is None else f'{x:.{n}f}'

def main():
 primary={};grid=[];sensitivity=[];temporal=[];placement=[];rows=[];numbers=[]
 for layout in LAYOUTS:
  for case in CASES:
   flow=read(layout,case+'_refined_metrics.json');tr=read(layout,case+'_refined_tracer_metrics.json');coarse=read(layout,case+'_tracer_metrics.json')
   assert flow['converged'] and flow['heat_input_W']==4000
   primary[layout+'_'+case]={'flow':flow,'tracer':tr}
   row=[LABEL[layout],LABEL[case],f(flow['fan_flow_m3_min']),f(tr['steady_equivalent_residence_time_s']),f(100*tr['retained_fraction_at_120s']),f(tr['occupied_zone_normalized_concentration_at_120s']*1000,3),f(tr['clearance_t95_s']),f(tr['clearance_occupied_zone_t95_s'])];rows.append(row)
   numbers.append({'layout':layout,'door_state':case,'estimated_flow_m3_min':flow['fan_flow_m3_min'],'mean_source_weighted_residence_s':tr['steady_equivalent_residence_time_s'],'retained_fraction_at120s':tr['retained_fraction_at_120s'],'occupied_zone_c_per_R_s_m3_at120s':tr['occupied_zone_normalized_concentration_at_120s'],'global_clearance_t95_s':tr['clearance_t95_s'],'occupied_zone_clearance_t95_s':tr['clearance_occupied_zone_t95_s']})
   grid.append({'layout':layout,'case':case,'coarse_h_m':.18,'primary_h_m':.14,'residence_relative_change_percent':(tr['steady_equivalent_residence_time_s']/coarse['steady_equivalent_residence_time_s']-1)*100,'occupied_concentration_relative_change_percent':(tr['occupied_zone_normalized_concentration_at_120s']/coarse['occupied_zone_normalized_concentration_at_120s']-1)*100,'coarse_occupied_clearance_t95_s':coarse['clearance_occupied_zone_t95_s'],'primary_occupied_clearance_t95_s':tr['clearance_occupied_zone_t95_s']})
   for Sc,sfx in [(.5,'_Sc0.5'),(.7,''),(1.,'_Sc1.0')]:
    tm=read(layout,case+'_refined'+sfx+'_tracer_metrics.json');sensitivity.append({'layout':layout,'case':case,'Sc_t':Sc,'mean_residence_s':tm['steady_equivalent_residence_time_s'],'occupied_c_per_R_s_m3':tm['occupied_zone_normalized_concentration_at_120s'],'occupied_clearance_t95_s':tm['clearance_occupied_zone_t95_s']})
   if case in ['both_open','both_closed']:
    tm=read(layout,case+'_refined_dt0125_tracer_metrics.json');temporal.append({'layout':layout,'case':case,'dt_coarse_s':.25,'dt_fine_s':.125,'retained_fraction_difference':tm['retained_fraction_at_120s']-tr['retained_fraction_at_120s'],'global_clearance_t95_difference_s':tm['clearance_t95_s']-tr['clearance_t95_s'],'occupied_clearance_t95_difference_s':tm['clearance_occupied_zone_t95_s']-tr['clearance_occupied_zone_t95_s']})
 for case in CASES:
  a,b=[primary[l+'_'+case]['tracer'] for l in LAYOUTS]
  placement.append({'case':case,'center_vs_corner_residence_change_percent':(b['steady_equivalent_residence_time_s']/a['steady_equivalent_residence_time_s']-1)*100,'center_vs_corner_occupied_concentration_ratio':b['occupied_zone_normalized_concentration_at_120s']/a['occupied_zone_normalized_concentration_at_120s'],'center_vs_corner_occupied_concentration_change_percent':(b['occupied_zone_normalized_concentration_at_120s']/a['occupied_zone_normalized_concentration_at_120s']-1)*100})
 summary={'primary':primary,'primary_rows':numbers,'placement_comparison':placement,'grid_comparison':grid,'tracer_diffusion_sensitivity':sensitivity,'time_step_comparison':temporal,'grid_independence_demonstrated':False,'experimentally_validated':False}
 (ROOT/'analysis_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2))
 for filename,data in [('comparison_summary.csv',numbers),('placement_comparison.csv',placement),('grid_comparison.csv',grid),('diffusion_sensitivity.csv',sensitivity),('time_step_comparison.csv',temporal)]:
  with (ROOT/filename).open('w',newline='',encoding='utf-8-sig') as fp:
   w=csv.DictWriter(fp,fieldnames=list(data[0]));w.writeheader();w.writerows(data)
 main_table=table(['位置','門狀態','排风 m³/min','平均滯留 τ s','120s滯留 %','走道指數 ×1000','全室t95 s','走道t95 s'],rows)
 placement_table=table(['門狀態','中央總滯留改變','中央／角落走道濃度','走道濃度改變'],[[LABEL[p['case']],f(p['center_vs_corner_residence_change_percent'])+'%',f(p['center_vs_corner_occupied_concentration_ratio'])+'倍',f(p['center_vs_corner_occupied_concentration_change_percent'])+'%'] for p in placement])
 grid_table=table(['位置／門狀態','τ加密改變','走道濃度加密改變','粗網格走道t95 s','細網格走道t95 s'],[[LABEL[g['layout']]+'／'+LABEL[g['case']],f(g['residence_relative_change_percent'])+'%',f(g['occupied_concentration_relative_change_percent'])+'%',f(g['coarse_occupied_clearance_t95_s']),f(g['primary_occupied_clearance_t95_s'])] for g in grid])
 ranges=[]
 for l in LAYOUTS:
  for c in CASES:
   group=[s for s in sensitivity if s['layout']==l and s['case']==c];ranges.append([LABEL[l]+'／'+LABEL[c],f(min(s['mean_residence_s'] for s in group))+'–'+f(max(s['mean_residence_s'] for s in group)),f(min(s['occupied_c_per_R_s_m3']*1000 for s in group),3)+'–'+f(max(s['occupied_c_per_R_s_m3']*1000 for s in group),3),f(min(s['occupied_clearance_t95_s'] for s in group))+'–'+f(max(s['occupied_clearance_t95_s'] for s in group))])
 sensitivity_table=table(['位置／門狀態','Sc假設下τ範圍 s','走道指數範圍','走道t95範圍 s'],ranges)
 time_table=table(['位置／門狀態','Δt由.25→.125 s','120s滯留比例差','全室t95差 s','走道t95差 s'],[[LABEL[t['layout']]+'／'+LABEL[t['case']],'加密',f(t['retained_fraction_difference'],10),f(t['global_clearance_t95_difference_s'],3),f(t['occupied_clearance_t95_difference_s'],3)] for t in temporal])
 vf=json.loads((ROOT/'validation.json').read_text())['checks'];pchecks=[vf[l+'_'+c+'_refined'] for l in LAYOUTS for c in CASES]
 validation_table=table(['位置／門狀態','節點／四面體','動量等相對殘差','熱方程相對殘差','RT0速度RMS改變'],[[LABEL[l]+'／'+LABEL[c],str(primary[l+'_'+c]['flow']['nodes'])+'／'+str(primary[l+'_'+c]['flow']['tetrahedra']),f"{primary[l+'_'+c]['flow']['nonlinear_relative_residual']:.2e}",f"{primary[l+'_'+c]['flow']['thermal_relative_residual']:.2e}",f(primary[l+'_'+c]['tracer']['velocity_transport']['relative_RMS_cell_velocity_correction']*100)+'%'] for l in LAYOUTS for c in CASES])
 left,right,both,shut=[primary['corner_'+c]['tracer'] for c in ['left_open','right_open','both_open','both_closed']]
 central=primary['center_both_open']['tracer'];maxheat=max(abs(v['direct_heat_balance_residual_W']) for v in pchecks);maxjump=max(v['RT0_internal_normal_flux_jump_max_m3_s'] for v in pchecks);maxtracer=max(v['tracer_transient_mass_balance_max_reference_units'] for v in pchecks)
 text=f'''# 廚房油煙滯留：靠左角落與檯面中央的有限元比較

**本報告比較兩種位置 × 四種門狀態，主要數據使用加密四面體網格。**雙燃氣爐各向空氣放熱2000 W，共4000 W；兩配置都按新三視圖重算，包含只在煙機90 cm寬度內的10 cm／45°擋板。舊版每爐600 W、煙機深60 cm的結果沒有拿來直接排序本次方案。

## 1. 主要判斷與優劣

在**無室外穿堂風、兩門外同壓、同一假設風機曲線、等量油煙源**的條件下，本次模型初步傾向**保留靠左位置，並以兩門全開作為補風工況**。它不是已驗證的實機結論。

- **移到中央沒有確定的整體滯留優勢。**開門三種工況的平均滯留只減少約1–2.4%，小於本次網格改變的幅度；部分擴散假設下左右單開的排序還會改變。全關時中央的平均滯留約少8.7%，但走道濃度增加約81%，仍有取捨。
- **走道的油煙示蹤濃度是中央方案較明顯的劣勢。**主模型中，中央／靠左的走道濃度比分別約3.11倍（雙開）、1.81倍（全關）、4.16倍（左開）、2.18倍（右開）。在Sc=0.5、0.7、1.0的已測假設下，中央較高的方向均保持；這個指標比幾%的總存量優勢更穩定。
- **靠左布局若只能單開一門，左開的使用中走道濃度較低。**右開約為左開的 {right['occupied_zone_normalized_concentration_at_120s']/left['occupied_zone_normalized_concentration_at_120s']:.2f} 倍。雙開和左單開的使用中走道均值很接近，不能把約6%的差異當成已確認性能提升。
- **兩門全開的尾端清除較快、較穩定；全關不利補風。**走道的長尾清除對混合與網格很敏感。不要用本報告的幾秒全室t95判定關火後可關煙機；全室總存量降到5%時，小量油煙仍可能留在走道。

**決策建議：**若主要目標是降低油煙向走道散佈，目前數據不足以支持僅為此目的把設備搬到中央。可先保留靠左配置，確認有效補風；如實際兩門全開形成穿堂風，需量測或加入室外壓差後重算，不能直接套用本次的雙開建議。若已選中央配置，應把補風路徑與走道濃度一併檢查。

## 2. 修正後的幾何與共同條件

| 項目 | 採用值與來源 |
| --- | --- |
| 房間 | 3.3 × 1.65 × 2.2 m，按三視圖解讀 |
| 檯面 | 高0.9 m、深0.6 m；可用長2.6 m，x=0.7–3.3 m |
| 高櫃 | x=0–0.7 m、深0.6 m，全高2.2 m |
| 上櫃／煙機占用帶 | 深0.3 m，z=1.6–2.2 m，實心外包絡近似 |
| 頂吸煙機 | 寬0.9 m、深0.3 m；底面距灶面0.7 m，使用者確認 |
| 擋板 | 前緣長0.10 m、向下45°，只在煙機0.9 m寬度內，隨煙機移動；厚度未提供，按2 mm近似 |
| 靠左煙機／兩爐 | 整體中心x=1.15 m；煙機x=0.7–1.6 m；爐頭x=0.95、1.35 m |
| 中央煙機／兩爐 | 使用者確認中心x=2.0 m；煙機x=1.55–2.45 m；爐頭x=1.8、2.2 m |
| 爐頭共同位置 | y=0.30 m、z=0.90 m；中心距0.40 m |
| 左／右門 | 寬0.9 m、高2.2 m；位於x=0、3.3 m；y=0.7–1.6 m位置仍由草圖估計 |
| 門關閉補風 | 每門底縫5 mm，沿用使用者同意的假設；全關入口鋪展至0.10 m帶 |
| 熱負荷 | 使用者更新每爐對空氣放熱2000 W，兩爐共4000 W，不是油煙排放率 |
| 空氣／牆壁 | 外部25°C；壁面以3 W/(m²·K)向25°C環境散熱，無滑移 |

兩配置的流體體積均約9.18282 m³。Gmsh以實體扣除建立貼合擋板表面的四面體網格；擋板僅沿煙機90 cm範圍，沒有橫跨整條2.6 m檯面。門片全開時移出流域，未解析開啟角度、門片阻擋、鍋具、人體、過濾網或排煙管內部。

![兩配置幾何模型](figures/geometry_comparison.png)

## 3. 風量與補風的設定

30 m³/min=0.5 m³/s與最大靜壓1300 Pa不是已知的同一工作點。沿用可檢查的二次假設曲線與孔口補風阻力：

```
p_fan = 1300 [1 − (Q/0.5)²]
Δp_makeup = ρ Q² / (2 Cd² A²), ρ=1.184 kg/m³, Cd=0.6
```

雙開／單開的補風面積為3.96／1.98 m²，估算風量接近30 m³/min。全關兩門底縫面積為0.009 m²，交點給出13.5468 m³/min、約1034.9 Pa的門縫網路壓降。煙機吸口假定為完整90×30 cm底面，以邊緣降為零的正弦法向速度歸一至Q，未取得實機有效濾網面積。未加入風管、濾網、外部風壓或其他漏風路徑；移位時也沒有另計風管改造，因此两位置使用同一條件式風量。

全關時，5 mm入口鋪展成0.10 m高帶，是室內遠場的近似。它保留補風量，沒有保留真實窄縫噴流的局部動量與混合，因而**全關的局部流場及精確清除秒數限制更大**。若完全密封，穩態質量守恆要求持續排風Q=0；本報告並非完全密封模型。

## 4. 油煙滯留如何計算

沒有油煙排放量、粒徑或食物／鍋具資料，故使用**無沉降的被動、氣體式示蹤物**作相對比較。兩爐等量釋放，總參考源強R=1參考單位/s；源位置與熱源相同，在爐面上方約7 cm的平滑區域，水平尺度0.085 m、垂直尺度0.055 m。這是未建模鍋具的源形狀假設。

示蹤物釋放0–120 s，之後停止示蹤物，**風機與4000 W熱負荷仍保持不變**，再追蹤900 s。此試驗隔離位置與補風路徑的影響，不是關火後的完整冷卻模擬。

| 指標 | 定義與解讀 |
| --- | --- |
| 歸一化濃度c/R | 單位s/m³；實際濃度需乘實際排放率R，不能當成已知mg/m³ |
| 平均滯留時間τ | 穩態室內示蹤物存量/R，即平均源加權停留時間；較低表示較少總存量 |
| 120s滯留比例 | 120s時室內存量/(R×120s)，分母為累積釋放量 |
| 走道濃度指數 | y>0.6 m且0.9<z<1.8 m的體積平均c/R×1000；不是任一人的個人暴露值 |
| 全室t95 | 停止示蹤物後，總存量降至停止時5%的時間 |
| 走道t95 | 同一試驗中走道平均濃度降至停止時5%的時間；可反映低濃度長尾 |

本模型中的門口大多是向內補風，門口油煙示蹤物淨外逸接近零。若唯一實質排出口是煙機，持續釋放的示蹤物最後排向煙機接近100%是守恆邊界的結果，**不能稱為實機捕煙率100%**。本報告因此以室內庫存、走道濃度與有限時間衰減比較優劣。

## 5. 八個工況：主網格h=0.14 m、Sc=0.7、Δt=0.25 s

{main_table}

t95欄是上述固定熱負荷的示蹤物試驗，數值秒數不是關火後的操作建議。走道均值很低時，低濃度尾端仍可能主導相對5%清除時間；讀取尾端時間必須併看源強未知與敏感性。

### 移到中央，相對靠左的改變

{placement_table}

開門工況的總存量略減，但中央的走道濃度上升約118–316%；全關時總存量下降約8.7%，走道濃度上升約81%。**兩種指標沒有同一方向的全面優勢。**左单開時主模型的中央走道清除較快，但其角落清除值在不同網格／擴散假設下大幅變動，不能視為已確認的移位收益。

![八工況空間分布](figures/smoke_eight_cases.png)

水平圖為z=1.5 m，使用共用對數色階；煙機虛線是投影，並非該高度的實體。較暗代表較小c/R，不是油煙完全不存在。3D點雲取原始FEM節點的濃度，點數與點大小不代表油煙質量；精確量取自積分與CSV。

![兩指標比較](figures/retention_comparison.png)

## 6. 尾端油煙滯留

全室存量多集中在熱源與吸入口附近，快速下降不表示所有位置同樣乾淨。主模型中，靠左布局兩門全開的走道t95為4.75 s，左單開為67.25 s，右單開為82.25 s；中央左／右單開為45.5／149.25 s。**其中單開與全關的尾端時間對混合與網格很敏感，這些精確秒數尚未得到確認。**

例如靠左左單開的走道t95，由粗網格3.75 s變成主網格67.25 s；Sc=0.5–1.0又可落在4–123.5 s。可看出低濃度死區／回流尾端的解析仍不足，不能依單一秒數作實際關機時機。雙開在已測兩網格、三個Sc下的走道尾端較短，這是相對較穩定的方向，但仍只適用無外部風壓的工況。

![全室存量衰減](figures/clearance_curves.png)

![走道尾端衰減](figures/walk_clearance_curves.png)

## 7. 有限元素與守恆輸運方法

氣流／熱浮力採三維四面體P1速度、P1壓力與P1溫度的有限元素離散，以scikit-fem執行。動量／連續性使用SUPG、PSPG、grad-div=0.5 m²/s以及門口回流能量耗散；偽時間迭代尋找穩態，不代表真实秒數。温度用保守Galerkin形式和對稱代數邊擴散，抑制非物理振盪并保持全域熱量守恆。

4000 W會產生較大的局部溫升，浮力使用理想氣體密度比的非線性修正：

```
f_buoyancy = g [1 − T₀,K/T_K] e_z
ν_eff = ν + ν_t, ν_t=0.005 m²/s
α_eff = α + ν_t/Pr_t, Pr_t=0.85
```

動量、連續性與熱容量仍採固定基準密度的不可壓縮近似；這**不是完整變密度低馬赫模型**，沒有在連續性中解析熱膨脹。局部最高空氣溫度約139–165°C只是平滑源附近的節點值，不能作為火焰、鍋面或人体温度。局部密度差可約30%，相關熱羽流需更完整模型與量測校準；總熱負荷雖已指定，源形狀、湍流混合與壁傳熱仍是假設。

原始P1穩定化速度只弱滿足不可壓縮，局部散度仍不小。示蹤物輸運前採RT0 H(div)面流量重建，以最小加權面流量修正滿足每個四面體進出平衡，維持煙機法向流量與壁面零通量。輸運用這個守恆速度，原始FEM速度另存VTU。這個重建會改變局部速度，不能把其嚴格質量平衡當成原流場已達高精度。

示蹤物的有限元素方程為：

```
∂c/∂t + ∇·(u_RT0 c) = ∇·(D_eff ∇c) + S
D_eff = D_m + ν_t/Sc_t, D_m=1e−5 m²/s
```

入口採乾淨空氣的Danckwerts總通量條件，出口為自然對流流出，壁面不吸附、無沉降。使用質量集塊後退Euler與對稱保守代數邊擴散；該附加擴散行／列和為零，守恆且抑制負濃度，但帶來網格依賴的數值混合。沒有解k–ε、k–ω、LES、油滴慣性、沉降、壁面沉積或燃燒化學。

## 8. 數值核對

粗網格每布局約2500節點／1.05萬四面體；主網格每布局約4500節點／2萬四面體，熱源附近尺寸由約9 cm加密至7 cm。h是全域目標尺度，不是所有四面體邊長。16個氣流／熱工況全部達到速度與溫度相對更新<2×10⁻⁵連續3次，另核對穩態方程殘差。

{validation_table}

從匯出四面體獨立積分體積、RT0內面法向通量連續性、外邊界風量、示蹤物流出／存量，並核對CSV、NPZ、VTU與全時序質量平衡，16組資料均通過設定容差。主工況最大直接熱量餘差為 {maxheat:.3e} W（對4000 W），最大共享面法向通量跳變為 {maxjump:.3e} m³/s，最大時序示蹤物平衡差為 {maxtracer:.3e} 參考單位。最低節點溫度為25°C、最低示蹤濃度非負。

上述核對確認離散方程與匯出一致，**沒有證明網格無關或實物驗證**。RT0速度體積加權RMS修正約17–30%，是本次流場離散尚有誤差的直接指標。

## 9. 網格、擴散與時間步長敏感性

### 兩網格比較：0.18→0.14 m

{grid_table}

總滯留幾%的移位差異與網格變化同級或更小；走道平均濃度加密後多下降，但中央高於角落的方向維持。尾端t95在部分工況變化非常大，未達網格無關。只有兩網格，未計算網格收斂指數，也不給統計置信區間。

### 示蹤物Sc=0.5、0.7、1.0

{ sensitivity_table }

改變Sc只改示蹤物擴散，氣流／熱場保持相同。中央走道濃度較高的方向在三個Sc下均保持；左右單開的總存量優勢有時變號。這是混合參數敏感性，並非已校準湍流不確定性的完整量化。尾端時間的寬範圍應保留，不能只公布一個精確清除秒數。

### Δt=0.25→0.125 s（雙開、全關）

{time_table}

120s庫存對時間步長的改變很小，全室t95略有變化；局部尾端的主要限制是空間與混合，而非單一時間步长。

## 10. 限制與實際用途

本次結果可用來比較指定模型內的相對滯留、補風路徑與油煙散佈，不能轉為實機捕煙率、實際mg/m³、健康限值或關火後必須延時多少秒。實際兩門外的壓力／溫度、鍋具與食物排放率、油煙粒徑、有效濾網面積、門縫與其他漏風路徑、風管阻力都會影響结果。

擋板已在兩配置中建模，但沒有無擋板對照組，因此沒有量化「擋板改善率」。壁面無示蹤物沉積，故不涵蓋長時間附著油膜或重新揚起。全關噴流被鋪展，局部高温熱膨脹未完整求解；本報告保留这些已知限制，不把精確小數當成實機準確度。

若用於真正搬移設備的工程決策，下一步應量測使用中總風量、兩門處風向／壓差及熱羽流／油煙示蹤濃度，並用實機風機曲線、鍋具源位置校準。依目前數據，**保留靠左、優先確保補風**比只憑中央總存量的微小下降搬移更有依據。

## 11. 交付檔案

| 檔案 | 內容 |
| --- | --- |
| `interactive_comparison.html` | 離線3D，可選位置、門、物理量、截面、示蹤物時間；數據直接來自原始四面體插值 |
| `report.pdf`／`REPORT.md`／`report.html` | 本中文比較報告 |
| `comparison_summary.csv` | 八主工況的風量、平均滯留、120s存量、走道濃度與兩種清除指標 |
| `placement_comparison.csv`／`grid_comparison.csv` | 移位差異及粗／細網格比較 |
| `diffusion_sensitivity.csv`／`time_step_comparison.csv` | Sc與時間步長比較 |
| `corner/config.json`、`center/config.json` | 兩配置的幾何、共同4000W與模型參數 |
| `corner/geometry.obj`、`center/geometry.obj` | 簡化實體幾何，包含僅90cm寬度的擋板，單位m |
| `<layout>/results/*_refined_combined.vtu` | 主網格、節點速度／溫度／壓力／示蹤濃度與單元RT0速度，可供ParaView讀取 |
| `*_refined_tracer.npz`／`*_refined_tracer_nodes.csv` | 原始網格、時間快照與逐點歸一化示蹤濃度 |
| `*_refined_tracer_history.csv` | 每0.25s的庫存、排出、質量餘差、走道／爐前濃度 |
| `probes.csv` | 隨布局移動的示例呼吸點和爐前探針；不是實際使用者站位 |
| `analysis_summary.json`／`validation.json` | 比較和独立核對完整數據 |
| `build_mesh.py`、`<layout>/solve.py`、`conservative_flux.py`、`tracer.py` | 可重現幾何、有限元素、RT0重建及油煙輸運程式 |

NPZ中濃度單位為s/m³／總源強；VTU的`velocity_FEM_m_s`是節點FEM速度，`velocity_RT0_transport_m_s`是逐單元輸運速度，二者不同。全室存量以源強等效秒表示；程式運算耗時不是物理示蹤時間。壓力以室內數值補風邊界作參考，不能直接等同門縫網路約1035 Pa壓降。

## 12. 重現計算

Python 3.12，見`requirements.txt`。Gmsh需要系統可用的GLU庫，中文圖表需要CJK字型。產生HTML報告用Pandoc，PDF另用Chromium／Playwright；已匯出的HTML與PDF可直接離線讀取。

```bash
python -m pip install -r requirements.txt
for layout in corner center
do
  python build_mesh.py --layout "$layout" --h .18
  python build_mesh.py --layout "$layout" --h .14
  python "$layout/solve.py" --case all --h .18
  python tracer.py --layout "$layout" --case all --dt .25
  for door in both_open both_closed left_open right_open
  do
    python "$layout/solve.py" --case "$door" --h .14 --suffix _refined --warm "$layout/results/$door.npz"
  done
  python run_tracer_verification.py --layout "$layout"
done
python validate_compare.py
python visualize_compare.py
python make_comparison_report.py
```

重算會覆寫同名結果。求解器每10次迭代保存檢查點，可以相同配置和`--warm <checkpoint.npz>`續算；配置改動後不能接用舊檢查點。舊版計算保留於原工作專案，本次成果另存於`kitchen_compare`。
'''
 # Uniform Traditional Chinese, including hand-authored occasional simplified characters.
 for a,b in {'排风':'排風','两':'兩','左单':'左單','真实':'真實','温度':'溫度','人體温度':'人體溫度','人体':'人體','并':'並','步长':'步長','结果':'結果','这些':'這些','独立':'獨立','时间':'時間','补':'補','带':'帶','热':'熱','烟':'煙'}.items():text=text.replace(a,b)
 (ROOT/'REPORT.md').write_text(text)
 css='''body{font-family:"Noto Sans CJK TC",sans-serif;max-width:1100px;margin:30px auto;padding:0 28px;line-height:1.7;color:#24384a;font-size:14px}h1{font-size:28px;color:#15344b}h2{font-size:21px;margin-top:32px;border-bottom:1px solid #d2dfe6;padding-bottom:6px}h3{font-size:17px}table{border-collapse:collapse;width:100%;font-size:12px;margin:14px 0}td,th{border:1px solid #d4dfe7;padding:7px;vertical-align:top;text-align:left}th{background:#eaf1f6}tr:nth-child(even){background:#f8fafc}img{max-width:100%}pre{white-space:pre-wrap;background:#f1f6f9;padding:14px;overflow-wrap:anywhere}code{font-size:.9em;overflow-wrap:anywhere}strong{color:#173f58}@page{size:A4;margin:16mm}@media print{body{max-width:none;margin:0;padding:0;font-size:9.5pt}table{font-size:7.5pt}tr,img{break-inside:avoid}h2,h3{break-after:avoid}pre{font-size:8pt}}'''
 (ROOT/'report.css').write_text(css)
 subprocess.run(['pandoc',str(ROOT/'REPORT.md'),'--standalone','--embed-resources','--css','report.css','--metadata','lang=zh-Hant','--metadata','pagetitle=廚房油煙滯留：角落與中央的有限元比較','-o',str(ROOT/'report.html')],cwd=ROOT,check=True)
 from playwright.sync_api import sync_playwright
 with sync_playwright() as p:
  browser=p.chromium.launch(executable_path=shutil.which('chromium'),headless=True,args=['--no-sandbox']);page=browser.new_page();page.set_content((ROOT/'report.html').read_text(),wait_until='load');page.pdf(path=str(ROOT/'report.pdf'),format='A4',print_background=True,prefer_css_page_size=True);browser.close()
 print('Chinese report and comparison tables complete',flush=True)
if __name__=='__main__':main()
