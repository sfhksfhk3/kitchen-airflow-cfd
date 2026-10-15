"""Build the Chinese analysis from completed results; embed all scientific plots."""
from pathlib import Path
import json,csv,re,base64,html,datetime
from zoneinfo import ZoneInfo
from markdown_it import MarkdownIt
from playwright.sync_api import sync_playwright
from pypdf import PdfWriter
import numpy as np
from emission_scenario import config,weights
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'deliverables';FIG=OUT/'figures_1to1'
DOORS=['both_open','both_closed','left_open','right_open'];DL={'both_open':'兩門全開','both_closed':'兩門全關','left_open':'左開右關','right_open':'右開左關'};LL={'corner':'靠左','center':'置中'}
def n(x,s=4):return '—' if x is None else f'{x:.{s}g}'
def table(headers,rows):return '| '+' | '.join(headers)+' |\n|'+'|'.join(['---']*len(headers))+'|\n'+'\n'.join('| '+' | '.join(map(str,r))+' |' for r in rows)+'\n\n'
def pic(name,caption):return f'\n![{caption}](figures_1to1/{name}.png)\n\n*{caption}*\n'

def main():
    data=json.loads((OUT/'summary.json').read_text());v=json.loads((ROOT/'verification/final_validation.json').read_text());phys=json.loads((ROOT/'verification/particle_physics_checks.json').read_text())
    if len(data['cases'])!=8 or not all(c['complete_saved_fields'] and c['complete_history'] for c in v['cases']):raise RuntimeError('Eight completed cases with endpoint fields are required')
    by={c['case']:c for c in data['cases']};date=datetime.datetime.now(ZoneInfo('Asia/Shanghai')).strftime('%Y-%m-%d')
    md=f'''# 廚房油滴與熱浮力：OpenFOAM 八工況比較

{date}｜條件式數值分析｜兩爐對空氣顯熱合計 4000 W

本次沿用已確認的廚房尺寸，以 OpenCFD OpenFOAM 2406 的有限體積法，計算靠左及置中位置、各四種門況。模型包含溫度—密度—氣流耦合、k–ω SST，以及分粒徑油滴的沉降、擴散、熱泳和壁面沉積。舊版有限元資料保留，這份報告的表格取自新的 OpenFOAM 計算。

以下濃度採文獻參考排放情境，不能視為此住家的實測濃度。未經校準的風機工作點、近壁湍流、油滴源分布及材料熱邊界仍會影響結果；軟體轉換本身不構成實驗驗證。

本輪關火後的熱流也隨時間改變；先前固定熱流場上的示蹤衰減和這次設定不同。兩輪數據的差值不能全部歸因於更換求解器或加入油滴。

## 1. 先看八個工況

走道區指廚房內檯面前方 y>0.6 m、z=0.9–1.8 m，不是門外走廊末端。走道濃度採烹調末段90–120秒的平均；滯留積分採完整觀測到240秒的空氣存量積分除以總釋放量，尚有殘留時是完整平均滯留時間的下限。

'''
    rows=[]
    for lay in ['corner','center']:
        for d in DOORS:
            c=by[f'{lay}_{d}'];rows.append([LL[lay],DL[d],n(c['Q_90_120_mean_m3_s']*60),n(c['walkway_90_120_oil_ug_m3']),n(c['retention_integral_observed_s']),n(c['hood_fraction_at_end']*100),n(c['deposited_fraction_end']*100)])
    md+=table(['位置','門況','抽風 m³/min','走道油滴 µg/m³','觀測滯留積分 s','經煙機離室 %','沉積 %'],rows)
    md+='\n濃度使用兩爐總油滴排放14.4 mg/min、質量中位徑2.5 µm、GSD1.5的參考情境。煙機排出比例表示顆粒經吸口離開計算室內，不是濾網淨化效率或未逸散的直接捕集率。\n'
    md+=pic('05_inventory_exposure_comparison_1to1','全室滯留與局部走道濃度分別比較；兩者衡量不同。')
    md+='\n## 2. 靠左與置中各有甚麼差別\n\n'
    ratios=[];paradoxes=[]
    for d in DOORS:
        a=by['corner_'+d];b=by['center_'+d]
        ratios.append([DL[d],n(b['walkway_90_120_oil_ug_m3']/max(a['walkway_90_120_oil_ug_m3'],1e-25))+' 倍',n(100*(b['retention_integral_observed_s']/a['retention_integral_observed_s']-1))+'%',n(b['walkway_90_120_oil_ug_m3']-a['walkway_90_120_oil_ug_m3'])])
        if b['retention_integral_observed_s']<a['retention_integral_observed_s'] and b['walkway_90_120_oil_ug_m3']>a['walkway_90_120_oil_ug_m3']:paradoxes.append(DL[d])
    md+=table(['門況','置中／靠左走道濃度','置中滯留積分相對變化','走道絕對差 µg/m³'],ratios)
    if paradoxes:md+='\n在'+ '、'.join(paradoxes)+'中，置中的觀測滯留積分較短，走道平均濃度卻較高。全室指標取決於所有顆粒留在房內的時間；局部走道濃度取決於有多少顆粒進入該區。少量顆粒偏向走道，已可提高局部濃度，並不要求全室總量同時增加。\n'
    md+='\n靠左的高櫃會限制源區左側的氣流空間；置中改變了兩側補風到煙機的路徑。因此，位置效果需要連同門況判斷。當走道濃度已接近零，大倍率也可能對應很小的絕對差，應同時讀取上表的絕對量。\n\n本輪保留每一門況的數值排序，但後段網格／時間步檢查只覆蓋部分工況，不能把小幅優勢定成普遍結論，也不能僅據這份模擬決定搬動灶台。補風路徑與實際運轉風量應優先確認。\n'
    md+='\n## 3. 關門為何明顯改變氣流\n\n兩扇門全關時，合計5mm底縫面積只有0.009 m²；按估計抽風工作點0.226 m³/s，等密度換算平均補風速度約25 m/s。完整熱流計算因排氣升溫膨脹，實際門縫速度略有差異。全開時補風面積大得多，入口平均速度低。高速近地噴流與熱羽流、櫃體和房間回流相互作用，不能只用「換氣次數」判斷灶台附近的流向。\n\n5mm開口採實際幾何，基準高度兩個單元；加密檢查採四個。這仍不足以宣稱門縫近壁流動已完全收斂。其他門框、窗戶及建築漏風在此未建模。\n'
    for lay in ['corner','center']:
        ordered=sorted((by[f'{lay}_{d}'] for d in DOORS),key=lambda c:c['walkway_90_120_oil_ug_m3'])
        low,high=ordered[0],ordered[-1]
        md+=f'\n在{LL[lay]}位置，本輪基準計算的走道末段濃度最低為{DL[low["doors"]]}（{n(low["walkway_90_120_oil_ug_m3"])} µg/m³），最高為{DL[high["doors"]]}（{n(high["walkway_90_120_oil_ug_m3"])} µg/m³）。這是本模型與排風假設下的比較，並非實測操作保證。\n'
    for d in DOORS:md+=pic(f'0{DOORS.index(d)+1}_oil_{d}_1to1',DL[d]+'：距地1.5m水平切面；圖中標示實際保存時刻，綠色叉是爐頭位置的平面投影。')
    md+='\n## 4. 熄火後衰減、空氣殘留與沉積\n\n120秒時停止對空氣放熱與油滴排放，煙機繼續運轉至240秒。沒有鍋具和板材的固體餘熱，因此這是理想化關火過程。下表的95%下降是相對停止時的各自存量或濃度；沒有到達的項目只記錄觀測窗內未到達。\n\n'
    rows=[]
    for lay in ['corner','center']:
        for d in DOORS:
            c=by[f'{lay}_{d}'];fmt=lambda v:n(v)+' s' if v is not None else '>120 s（未達）'
            rows.append([LL[lay],DL[d],fmt(c['airborne_t95_after_stop_s']),fmt(c['walkway_t95_after_stop_s']),n(100*c['remaining_fraction_end']),n(100*c['door_escape_fraction_end'])])
    md+=table(['位置','門況','全室95%下降','走道95%下降','240s仍在空氣 %','經兩門離室 %'],rows)
    md+='\n經兩門離室只是離開廚房計算域，並不代表整戶污染已消失；本模型沒有門外住宅空間。壁面沉積會降低空氣濃度，同時留下表面油污。精確延遲關機秒數仍受網格、混合、烹調源和實際補風影響，表中時間不應直接當成家用煙機的操作延時建議。\n'
    md+=pic('06_clearance_four_panels_1to1','四幅合併：靠左／置中 × 全室空氣存量／走道濃度，保留絕對量級。')
    md+=pic('07_normalized_clearance_1to1','各曲線各自起點歸一化；相對比例的交叉不等於絕對濃度同時交叉。')
    md+='\n## 5. 從氣體示蹤轉為油滴模型\n\n熱流方程使用理想氣體密度ρ(p,T)、溫度相關黏度、動量及能量耦合。油滴則為稀薄球粒子，以五個獨立粒徑的質量分率Y分別求解：\n\n`∂(ρY)/∂t + ∇·[(ρu + ρv漂移)Y] = ∇·[ρ(DB + νt/Sc_t)∇Y] + S − 壁面沉積`\n\n漂移包含Stokes–Cunningham沉降與Talbot熱泳。DB為Stokes–Einstein布朗擴散。壁面通量由含van Driest阻尼的擴散阻力積分，加向壁的沉降／熱泳通量形成，以正的隱式損失項求解；黏附率設為1。這不是完整油膜模型，沒有凝結、蒸發、成核、凝併、再懸浮或慣性撞擊。\n\n'
    md+=table(['粒徑 µm','25°C布朗擴散 m²/s','沉降速度 m/s','靜止空氣下降1m約需 s'],[[n(p['diameter_um']),n(p['brownian_D_m2_s']),n(p['settling_velocity_m_s']),n(p['settling_1m_time_s'])] for p in phys['particle_properties_25C']])
    md+=f'\n參考排放量相對0.5m³/s空氣流量的質量流率比約為{phys["bulk_particle_to_air_mass_flow_ratio_at_0p5m3_s"]:.2g}，採油滴不回饋氣流的稀薄近似；溫度、密度與氣流仍耦合。油煙有顆粒相，不代表本案例必須假定很強的顆粒動量回饋。\n'
    md+=pic('08_size_deposition_1to1','各粒徑相同參考源強的沉積比較。黏附率1只表示已到達壁面者全黏住，不能作完整油污沉積的上限保證。')
    md+='\n## 6. 文獻參數、源強及其不確定性\n\n五個單粒徑各以兩爐合計1µg/s獨立計算，原始濃度以`C/R`（s/m³）交付。它們沒有被假定為五等分的實際油煙。為便於讀圖，再以總質量分布線性組合：\n\n'
    md+=table(['項目','基準／採用方式'],[['對空氣顯熱','使用者確認：每頭2000W，共4000W'],['油滴密度','850kg/m³，文獻代表性近似'],['參考總排放','每頭7.2mg/min，兩頭14.4mg/min；由商用廚房CFD的0.12mg/s油滴源轉用'],['質量中位徑','2.5µm，參考六道中式菜餚研究'],['GSD','本次假設1.5，另以1.2和2.0重新組合；不是實測分布寬度'],['湍流擴散','Sc_t=0.7，未實測校準'],['油滴導熱','0.17W/(m·K)，僅用於熱泳係數']])
    md+='\n[Chen等2018年住宅中式烹調研究](https://pubmed.ncbi.nlm.nih.gov/29302961/)報告的PM2.5平均排放2.056±3.034mg/min，涵蓋不同烹調方法，可作另一個量級參考；本報告沒有把它直接當成所有油滴總質量。參考源強7.2mg/min/頭來自[2022年商用中式廚房CFD](https://www.mdpi.com/2076-3417/12/17/8896)；中位徑參考[Chen等2023年研究](https://doi.org/10.1007/s11356-023-30380-4)。這些跨場景轉用及GSD假設會影響絕對濃度。\n\n相同熱流下，排放量變成0.25倍或4倍，稀薄線性模型的濃度和沉積質量也等比例縮放。相對位置比較須保持兩方案源強一致。質量分布寬度變更可用已有五個粒徑響應重組，不需要重跑氣流。\n\n此粒徑範圍未涵蓋全部小於0.1µm的粒子。不能用這組質量圖推算完整超細粒子數、VOC或健康風險。完整來源與採用方式見`references/SOURCES_zh.md`。\n'
    md+='\n## 7. 數值核對與仍未解決的限制\n\n'
    md+=table(['工況','顆粒質量相對誤差上限','氣體質量誤差 kg/s','烹調段熱量殘差RMS W（占4000W）'],[[c['case'],n(c['particle_mass_error_max_relative']),n(c['air_mass_balance_max_kg_s_after10']),n(c['energy_rms_after10_cooking_W'])+' ('+n(c['energy_rms_after10_cooking_W']/40)+'%)'] for c in v['cases']])
    md+='\n顆粒質量核對和熱流核對必須分開解讀：顆粒收支很小的誤差不保證熱流場已收斂。表中的熱量殘差以逐步體積積分和邊界通量重建，包含離散／迭代與診斷近似的差異；雙關工況的熱量殘差明顯較大，因此熱浮力和局部濃度預測仍有數值限制。本輪交付可檢查的比較及敏感性結果，不宣稱完成高精度熱流收斂驗證。\n'
    md+=f'\n近壁阻力積分以72個對數區間計算，對獨立自適應積分的最大沉積速度誤差約{100*phys["wall_quadrature_max_abs_relative_error_72"]:.2f}%。這只核對積分算法，並未驗證近壁沉積的物理模型。\n'
    meshes=json.loads((ROOT/'verification/mesh_inventory.json').read_text())
    md+='\n嚴格網格檢查仍有未通過項目；以下保留實際數量，不能將基本幾何可運算解讀為全部網格品質合格。\n\n'
    md+=table(['網格','單元數','小行列式單元數','最小面插值權重','未通過嚴格項目數'],[[m['case'],int(m['cells']),int(m['small_determinant_cells']),n(m['minimum_face_interpolation_weight']),m['strict_failed_check_count']] for m in meshes])
    ss=v['sensitivity']
    if ss:
        groups={}
        for s in ss:groups.setdefault(s['test'],[]).append(s)
        md+=table(['敏感性案例','五粒徑存量最大變化','走道濃度最大變化','110–120s熱量殘差RMS：基準→檢查 W'],[[name,n(100*max(abs(x['stock_relative_change']) for x in g))+'%',n(100*max(abs(x['walkway_relative_change']) for x in g))+'%',n(g[0]['parent_energy_RMS_110_120_W'])+' → '+n(g[0]['test_energy_RMS_110_120_W'])] for name,g in groups.items()])
        gridrows=[];gridmixed={}
        for name,g in groups.items():
            if g[0]['kind']!='grid':continue
            g=sorted(g,key=lambda x:x['diameter_um']);ww=weights()
            base=sum(w*x['parent_walkway_mean_R'] for w,x in zip(ww,g))*240
            test=sum(w*x['test_walkway_mean_R'] for w,x in zip(ww,g))*240
            gridmixed[name]=(base,test)
            gridrows.append([name,n(base),n(test),n(100*(test/base-1))+'%'])
        if gridrows:
            md+='\n相同110–120秒窗口，參考質量分布的走道濃度如下。這個窗口和主比較表的90–120秒不同，不能直接逐項相減。\n\n'
            md+=table(['雙關網格檢查','基準 µg/m³','加密 µg/m³','相對變化'],gridrows)
            a=gridmixed.get('corner_both_closed_grid');b=gridmixed.get('center_both_closed_grid')
            if a and b:
                if (a[0]-b[0])*(a[1]-b[1])<0:
                    md+='\n**雙關工況的靠左／置中走道濃度排序在加密後反轉，這一位置優劣未通過本輪網格敏感性檢查。不可據基準表格宣稱其中一種位置已被證明較好。**\n'
                else:
                    md+='\n雙關工況在這個後段窗口的相對排序未反轉，但絕對濃度的變化仍應作為數值不確定性解讀，不能稱作已網格獨立。\n'
    md+='\n時間步與加密網格檢查從相同保存狀態的烹調後段續算，觀察110–120秒及關火過渡。網格轉移採正權重插值，逐粒徑修正使總質量一致，再先校正質量通量；加密高度四個門縫單元，近地層逐步增厚。這是後段數值敏感性，不是獨立完整歷程的網格收斂試驗，也不是統計信賴區間。局部濃度對這些設定較敏感時，位置排序應視為暫定。\n\n網格基本檢查通過；嚴格幾何檢查對門縫附近的扁長單元仍有單元行列式、插值權重或梯度支撐警告，部分類別涉及多個近地層，不能稱作所有嚴格檢查全通過。實際數量與位置留在checkMesh記錄。入口湍流k=0.001m²/s²、ω=1s⁻¹為未校準假設；高速門縫、壁面函數、未建模漏風及外部風壓仍限制結果。\n\n入口靜溫設25°C，未對門縫高速流作總溫修正，約有0.3°C量級差異。牆面採向25°C環境的3W/(m²K)簡化散熱邊界，沒有顆粒板／不鏽鋼導熱或熱容量，亦無火焰／鍋具直接輻射；本輪空氣溫度不能用來判定板材的實際最高溫度或耐熱安全性。\n'
    md+='\n## 8. 幾何、風機及原始資料\n'
    md+=pic('00_geometry_comparison_1to1','尺寸與位置均沿用確認的三視圖；煙機前緣擋板採零厚度擋流面。')
    md+='\n風機只給最大風量30m³/min、最大靜壓1300Pa，因此仍採先前的二次曲線與孔口係數0.6估算工作流量。OpenFOAM以均勻吸口負壓控制該體積流量，讓局部吸入分布自行形成。完整風機、煙管和濾網未聯立求解；室內控制負壓不能稱作廠商靜壓的完整預測。\n\n初期指定吸口逐點速度的試算在開口邊緣出現非物理高風速，未採用。最終模型改用自由出口速度與均勻負壓；高速度集中在符合流量／面積量級的門縫。被拒絕的試算仍保留在工作區的verification資料夾，未納入本報告表格。\n\n交付資料包含：每一粒徑的逐時空氣存量與排出／沉積、參考混合情境時序、單元中心CSV、原生OpenFOAM檢查點、VTU三維場、離線互動3D、數值檢查及GPL求解器原始碼。續算造成的保存時刻偏移有在圖中明示實際時間，末段場不冒稱為精確120秒快照；定量表格讀取完整時序。\n\n重現及續算方法見`README_zh.md`。原生場為有限體積單元中心資料，不是有限元節點值。\n'
    md+='\n需特別注意：以上是給定估計排風量的條件式室內比較。二次曲線／孔口網路與實際三維門縫的阻力並未在每一時間步重新聯立，因此不能把兩者各自的壓降當作同一完整風機工作點；30m³/min也不等於裝機後必達的風量。若實際煙管、止回閥及濾網令風量下降，濃度與位置優劣都可能改變。\n\n熱量與油滴源採同一歸一化高斯體積分布：各爐水平標準差8.5cm、垂直標準差5.5cm，中心高於灶面7cm，截於灶面上方0–25cm及檯面深度內。每爐的體積積分分別保證2000W，但分布形狀未以實際鍋具量測校準，亦沒有炒鍋動作的瞬間噴發。\n'
    # Save markdown and embed images for browser/PDF portability.
    (OUT/'openfoam_oil_analysis_zh.md').write_text(md)
    body=MarkdownIt('commonmark',{'html':True}).enable('table').render(md)
    def embed(m):
        src=m[1];p=OUT/src
        if not p.exists():return m[0]
        return 'src="data:image/png;base64,'+base64.b64encode(p.read_bytes()).decode()+'"'
    body=re.sub(r'src="([^"]+)"',embed,body)
    css='''@page{size:A4;margin:16mm 15mm 17mm}body{font:10pt/1.7 "Noto Sans CJK TC","Noto Sans CJK JP",sans-serif;color:#20313c}h1{font-size:22pt;line-height:1.4}h2{font-size:15pt;border-bottom:1px solid #b9cbd4;padding-top:10pt;break-after:avoid}p{margin:7pt 0}table{border-collapse:collapse;width:100%;font-size:8pt;margin:10pt 0}td,th{border:1px solid #ccd8de;padding:4pt;vertical-align:top}th{background:#e9f0f4}tr{break-inside:avoid}img{max-width:100%;display:block;margin:14pt auto 4pt}p:has(>img){break-before:page;break-inside:avoid}h2:has(+p>img){break-before:page}h2+p:has(>img){break-before:auto}em{font-size:9pt;color:#526974}code{font-family:monospace;font-size:8pt;white-space:normal;overflow-wrap:anywhere}a{color:#1b6c8f}'''
    page='<!doctype html><html lang="zh-Hant"><meta charset="utf-8"><style>'+css+'</style><body>'+body+'</body></html>'
    (OUT/'openfoam_oil_analysis_zh.html').write_text(page)
    with sync_playwright() as p:
        b=p.chromium.launch(executable_path='/usr/bin/chromium',headless=True,args=['--no-sandbox']);tab=b.new_page();tab.set_content(page,wait_until='load');tab.pdf(path=str(OUT/'openfoam_oil_analysis_zh.pdf'),print_background=True,prefer_css_page_size=True);b.close()
    pdfs=sorted(FIG.glob('*_1to1.pdf'));writer=PdfWriter()
    for p in pdfs:writer.append(str(p))
    writer.write(str(OUT/'all_figures_1to1.pdf'))
    (OUT/'figures_gallery.html').write_text('<!doctype html><html lang="zh-Hant"><meta charset="utf-8"><title>1:1 圖集</title><body style="font-family:system-ui;max-width:1000px;margin:auto"><h1>OpenFOAM 1:1 圖集</h1>'+''.join(f'<h2>{p.stem}</h2><a href="figures_1to1/{p.stem}.pdf">PDF</a> · <a href="figures_1to1/{p.stem}.svg">SVG</a><img style="width:100%" src="figures_1to1/{p.stem}.png">' for p in pdfs)+'</body></html>')
    print('Report and square figure collection complete',flush=True)

if __name__=='__main__':main()
