"""Report derived from the completed 4000 W corner / closed-door FEM case."""
from pathlib import Path
import json, base64
import numpy as np
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'deliverables'
SOURCE=ROOT.parent/'corner/results'

def main():
    summary=json.loads((OUT/'wall_summary.json').read_text())
    wall=summary['4000W_refined']
    grid=summary['grid_comparison_h018_to_h014']
    records=[]
    for label, name, key in [('粗網格0.18 m','both_closed','4000W_coarse'),('細網格0.14 m','both_closed_refined','4000W_refined')]:
        m=json.loads((SOURCE/f'{name}_metrics.json').read_text())
        records.append((label,m,summary[key]))
    fine=records[-1][1]
    numerical={
        'both_runs_converged':all(m['converged'] for _,m,_ in records),
        'thermal_relative_residual_below_1e-4':all(m['thermal_relative_residual']<1e-4 for _,m,_ in records),
        'global_flux_relative_error_below_1e-8':all(m['net_flux_relative_error']<1e-8 for _,m,_ in records),
        'global_energy_relative_error_below_1e-6':all(abs(m['discrete_heat_balance_residual_W'])/4000<1e-6 for _,m,_ in records),
        'surface_integral_matches_robin_identity':abs(wall['total_net_heat_into_wall_W']-3*wall['area_m2']*(wall['area_average_temperature_C']-25))<1e-10,
        'face_csv_heat_sum_matches_report':abs(np.genfromtxt(OUT/'4000W_refined_wall_faces.csv',delimiter=',',names=True)['net_heat_into_wall_W'].sum()-wall['total_net_heat_into_wall_W'])<1e-10,
    }
    numerical={k:bool(v) for k,v in numerical.items()}
    if not all(numerical.values()):raise RuntimeError('Numerical verification failed')
    (OUT/'numerical_verification.json').write_text(json.dumps(dict(numerical_checks=numerical,grid=grid,physical_limitations=dict(solid_conduction_not_solved=True,direct_radiation_not_solved=True,actual_solid_temperature_not_validated=True,fixed_reference_density=True)),ensure_ascii=False,indent=2))
    rows='\n'.join(f"| {label} | {m['tetrahedra']:,} | {w['total_net_heat_into_wall_W']:.2f} | {w['area_average_temperature_C']:.2f} | {w['maximum_interface_temperature_C']:.2f} |" for label,m,w in records)
    pos=wall['maximum_temperature_position_m']
    report=f'''# 靠左雙爐、兩門全關：高櫃側面受熱估算

最新工況為每個爐頭對空氣放熱 **2000 W**，兩頭同時使用，共 **4000 W**。研究面為爐頭旁高櫃的右側垂直面（**x=0.7 m**）。此工況與已完成的有限元計算一致，因此直接沿用已存流場與溫度場，新增側面積分及可視化；9000 W中途計算已停止，其結果未用於本報告。

**以下是空氣有限元模型的簡化散熱邊界結果。沒有求解櫃板固體導熱，也沒有加入火焰／鍋具直接輻射；板材、厚度及背面散熱條件尚未提供。模型邊界溫度不能當作已驗證的實際板溫，功率也不是包括直接輻射的全部受熱量。**

## 主要結果

| 指標 | 0.14 m細網格結果 |
|---|---:|
| 空氣接觸的側面面積 | {wall['area_m2']:.6f} m² |
| 該面模型淨受熱功率 | **{wall['total_net_heat_into_wall_W']:.2f} W** |
| 面積平均邊界溫度 | **{wall['area_average_temperature_C']:.2f} °C** |
| 局部最高邊界溫度 | **{wall['maximum_interface_temperature_C']:.2f} °C** |
| 局部最低邊界溫度 | {wall['minimum_interface_temperature_C']:.2f} °C |
| 面積平均淨熱通量 | {wall['area_average_heat_flux_W_m2']:.2f} W/m² |
| 最大節點淨熱通量 | {wall['maximum_nodal_heat_flux_W_m2']:.2f} W/m² |
| 最高溫座標(x,y,z) | ({pos[0]:.3f}, {pos[1]:.3f}, {pos[2]:.3f}) m |

W是整個面每秒得到的淨熱量；W/m²是局部熱量密度。正值代表空氣域向該面／散熱參考環境傳熱。最高溫為單一邊界節點的局部熱點，不能代表整塊板。

計算面包含 y=0–0.6 m、z=0.9–1.6 m，以及 y=0.3–0.6 m、z=1.6–2.2 m。灶台以下、上櫃遮擋部分不與空氣直接接觸；名義暴露面積0.600 m²，扣除擋板接合处後為{wall['area_m2']:.6f} m²。

![模型邊界溫度分布](wall_temperature_4000W.png)

## 工況與積分方法

房間3.3×1.65×2.2 m，灶面高0.9 m。兩爐中心x=0.95和1.35 m、y=0.30 m；左爐中心距所選側面0.25 m。頂吸煙機寬0.9 m、深0.3 m，底面高1.6 m，距灶面0.7 m；其90 cm寬度內有10 cm、向下45°前擋板。

兩門全關，但每門底保留5 mm縫。沿用假設風機曲線p=1300[1−(Q/0.5)²]、孔口係數0.6，運行風量{fine['fan_flow_m3_min']:.3f} m³/min。標稱30 m³/min與最大靜壓1300 Pa不是同一運行點。門縫在室內網格中鋪展成10 cm高補風帶，未解析真實5 mm門縫射流，也未加入風管／濾網阻力。

三維四面體P1/P1有限元，採SUPG/PSPG、grad-div、非線性理想氣體浮力修正，以及保守熱量弱式和代數邊擴散。空氣方程仍為固定參考密度不可壓縮模型，有效渦黏度0.005 m²/s、湍流Pr=0.85。每爐2000 W以近灶面高斯體積熱源施加，沒有火焰或鍋具幾何。

所有固體邊界的散熱參考環境為25°C，等效散熱係數為3 W/(m²·K)。在所選表面的每個三角形上積分：

`局部淨熱通量 q'' = 3 × (T_boundary − 25)`

`總功率 Q = Σ A_triangle × 3 × (T_triangle_mean − 25)`

線性三角形面積積分是精確的。這個3 W/(m²·K)是原模型的等效邊界參數，不是櫃板材料導熱係數。輸出為穩態估算，不包括板材儲熱或升溫時間。

## 數值驗證與適用範圍

| 網格 | 四面體數 | 側面功率W | 平均邊界溫度°C | 最高邊界溫度°C |
|---|---:|---:|---:|---:|
{rows}

粗到細網格，側面功率相對細網格值改變 **{grid['total_net_heat_into_wall_W']['relative_change_percent']:+.2f}%**，平均溫升改變 **{grid['temperature_rise_average_relative_change_percent']:+.2f}%**，最高邊界溫度改變 **{grid['maximum_interface_temperature_C']['relative_change_percent']:+.2f}%**。兩級網格尚不足以證明網格無關；這些差值不是包含物理模型誤差的完整誤差界限。可將約40–44 W、平均約47–50°C、局部最高約128–133°C作為這兩級網格的條件結果範圍。

兩級已存計算均收斂。細網格非線性相對殘差{fine['nonlinear_relative_residual']:.3e}、熱量相對殘差{fine['thermal_relative_residual']:.3e}、全域體積流量相對誤差{fine['net_flux_relative_error']:.3e}。表面CSV逐三角形功率總和與報告一致，也滿足Q=3A(T平均−25)恆等式。

4000 W熱量收支：煙機帶走{fine['advective_heat_flux_W']['hood']:.4f} W，**全部**固體邊界合計散熱{fine['wall_heat_loss_W']:.4f} W，門口熱交換反力{fine['door_thermal_reaction_W']:.6f} W，離散熱收支誤差{fine['discrete_heat_balance_residual_W']:.3e} W。高櫃側面的{wall['total_net_heat_into_wall_W']:.2f} W只是全部固體邊界散熱的一部分。

全域守恆不等於每個單元強形式不可壓縮；P1/P1速度場散度／速度梯度L2比值為{fine['divergence_gradient_norm_ratio']:.3f}。固定湍流參數、數值擴散、假設風機曲線和門縫處理均影響熱量分配。局部最高空氣溫度{fine['max_nodal_temperature_C']:.1f}°C，密度變化已顯著；浮力項有非線性修正，但輸運仍用固定參考密度。因此本結果是模型條件估算，沒有實測驗證，不能用來確定實際板材最高溫或熱耐受能力。

若要取得實際櫃板溫度，需側板材質、厚度及背面空間／散熱條件，並加入固體導熱耦合。若要求包含火焰／鍋具直接輻射的全部受熱功率，還需相關輻射条件或可校準的量測。

## 資料與可視化

- `wall_summary.csv`：兩級網格的側面功率與溫度摘要。
- `4000W_refined_wall_nodes.csv`：表面座標、邊界溫度、局部熱通量。
- `4000W_refined_wall_faces.csv`：三角形面積、平均溫度、熱通量及功率。
- `4000W_refined_wall.vtu`：ParaView可開啟的3D表面資料。
- `wall_interactive_3d.html`：離線旋轉、縮放、查讀及兩級網格切換。
- 資料包另含原始4000 W空氣場、收斂記錄、模型設定及重現程式。
'''
    diagnostic_path=OUT/'hotspot_diagnostic.json'
    if diagnostic_path.exists():
        diag=json.loads(diagnostic_path.read_text())
        appendix=f'''## 熱點位置的補充核查

左爐的中心是(x,y,z)=(0.95,0.30,0.90) m，在x=0.70 m側面上的投影為(y,z)=(0.30,0.90) m。25 cm間距沿圖外x方向，不是圖中的y方向。原圖青色叉位於(y,z)=(0,1.04) m，標的是原模型側面最大值，不是爐頭位置。

既有細網格在爐頭旁x=0.7–1.1 m、y=0.2–0.4 m、z=0.9–1.3 m區域的平均速度y分量約−0.223 m/s，表示有向後牆輸運；熱點不必與爐头投影重合。然而流動輸運與数值誤差的影響尚不能完全分開。

另用既有RT0局部質量守恆速度重建場做固定流場熱方程敏感性檢查。沒有重新耦合浮力與流動，因此它不是替代原結果的新完整解。側面最高值改為{diag['RT0_wall_max_C']:.1f}°C，位置(x,y,z)=({diag['RT0_wall_max_xyz_m'][0]:.3f},{diag['RT0_wall_max_xyz_m'][1]:.3f},{diag['RT0_wall_max_xyz_m'][2]:.3f}) m；側面功率改為{diag['RT0_wall_heat_W']:.1f} W，平均邊界溫度{diag['RT0_wall_average_boundary_C']:.1f}°C。

這顯示原127.6°C熱點位置與44.5 W受熱量都對速度場數值處理敏感。兩級網格收斂與全域熱量守恆不足以確認它是真實熱點；原表格保留為原模型數據，不能當成已驗證的實際板溫或受熱量。113.9°C同樣是數值敏感性結果，不能作實際板溫結論。

![爐頭投影與熱點位置核查](wall_hotspot_explanation.png)

'''
        report=report.replace('## 資料與可視化',appendix+'## 資料與可視化')
    (OUT/'wall_heat_report_zh.md').write_text(report)
    from markdown_it import MarkdownIt
    body=MarkdownIt('commonmark').enable('table').render(report)
    plot='data:image/png;base64,'+base64.b64encode((OUT/'wall_temperature_4000W.png').read_bytes()).decode()
    body=body.replace('src="wall_temperature_4000W.png"',f'src="{plot}"')
    extra_image=OUT/'wall_hotspot_explanation.png'
    if extra_image.exists():
        embedded='data:image/png;base64,'+base64.b64encode(extra_image.read_bytes()).decode()
        body=body.replace('src="wall_hotspot_explanation.png"',f'src="{embedded}"')
    css='body{font-family:"Noto Sans CJK TC",sans-serif;color:#182535;font-size:11px;line-height:1.65;margin:30px auto;max-width:950px}h1{font-size:22px}h2{font-size:16px;margin-top:24px;break-after:avoid}table{border-collapse:collapse;width:100%;font-size:10px}td,th{border:1px solid #c9d0d7;padding:6px}th{background:#edf2f7}img{max-width:100%;max-height:650px;display:block;margin:auto}code{white-space:pre-wrap}tr{break-inside:avoid}@page{size:A4;margin:15mm}'
    document=f'<!DOCTYPE html><html lang="zh-Hant"><meta charset="utf-8"><title>高櫃側面受熱估算</title><style>{css}</style><body>{body}</body></html>'
    (OUT/'wall_heat_report_zh.html').write_text(document)
    from playwright.sync_api import sync_playwright
    with sync_playwright() as pw:
        browser=pw.chromium.launch(executable_path='/usr/bin/chromium',args=['--no-sandbox'])
        page=browser.new_page();page.set_content(document,wait_until='load');page.evaluate('document.fonts.ready')
        page.pdf(path=str(OUT/'wall_heat_report_zh.pdf'),format='A4',print_background=True)
        browser.close()
    print(json.dumps(dict(wall=wall,verification=numerical),ensure_ascii=False,indent=2))

if __name__=='__main__':main()
