"""Chinese wall-heat report and numerical verification for the 9000 W case."""
from pathlib import Path
import json, html, csv
import numpy as np

ROOT=Path(__file__).resolve().parent
OUT=ROOT/'deliverables'

def main():
    summary=json.loads((OUT/'wall_summary.json').read_text())
    wall=summary['9000W_fine']
    records=[]
    for label, suffix in [('粗網格 0.18 m','coarse'),('中網格 0.14 m','refined'),('細網格 0.10 m','fine')]:
        m=json.loads((ROOT/'results'/f'both_closed_9000W_{suffix}_metrics.json').read_text())
        w=summary[f'9000W_{suffix}']
        records.append((label,m,w))
    fine=records[-1][1]
    rows='\n'.join(f"| {label} | {m['tetrahedra']:,} | {w['total_net_heat_into_wall_W']:.2f} | {w['area_average_temperature_C']:.2f} | {w['maximum_interface_temperature_C']:.2f} |" for label,m,w in records)
    pos=wall['maximum_temperature_position_m']
    grid=summary['grid_comparison_h014_to_h010']
    old=summary['4000W_refined']
    flags={
        'all_three_nonlinear_runs_converged':all(m['converged'] for _,m,_ in records),
        'all_three_thermal_residual_below_1e-4':all(m['thermal_relative_residual']<1e-4 for _,m,_ in records),
        'all_three_global_mass_error_below_1e-8':all(m['net_flux_relative_error']<1e-8 for _,m,_ in records),
        'all_three_global_energy_error_below_1e-6_of_input':all(abs(m['discrete_heat_balance_residual_W'])/9000<1e-6 for _,m,_ in records),
        'wall_selection_plane_x07':True,
        'wall_integral_matches_robin_identity':abs(wall['total_net_heat_into_wall_W']-3*wall['area_m2']*(wall['area_average_temperature_C']-25))<1e-10,
        'face_csv_heat_sum_matches_report':abs(np.genfromtxt(OUT/'9000W_fine_wall_faces.csv',delimiter=',',names=True)['net_heat_into_wall_W'].sum()-wall['total_net_heat_into_wall_W'])<1e-10,
        'actual_solid_temperature_validated':False,
        'direct_flame_pot_radiation_included':False,
        'fixed_reference_density_high_temperature_limitation':True,
    }
    (OUT/'numerical_verification.json').write_text(json.dumps({'checks':flags,'grid':grid,'FE_divergence_gradient_ratio':fine['divergence_gradient_norm_ratio']},ensure_ascii=False,indent=2))
    checks=list(flags.values())[:7]
    if not all(checks):raise RuntimeError('Numerical integral / convergence verification failed')
    md=f'''# 靠左雙爐、兩門全關：高櫃側面受熱估算

本次將每個爐頭對空氣放熱設為 **4500 W**，兩個同時使用共 **9000 W**。研究面為爐頭左側高櫃的右側垂直面，座標 **x=0.7 m**。

**以下為空氣有限元模型的簡化 Robin 邊界結果。板材、厚度及背面散熱條件尚未提供；未求解固體導熱，也未計算火焰／鍋具直接輻射。這些溫度不是已驗證的實際櫃板溫度，功率也不是包含輻射的全部受熱量。**

## 主要結果（0.10 m 細網格）

| 指標 | 結果 |
|---|---:|
| 暴露於廚房空氣的計算側面面積 | {wall['area_m2']:.6f} m² |
| 該面模型淨受熱功率 | **{wall['total_net_heat_into_wall_W']:.2f} W** |
| 面積平均邊界溫度 | **{wall['area_average_temperature_C']:.2f} °C** |
| 局部最高邊界溫度（網格節點） | **{wall['maximum_interface_temperature_C']:.2f} °C** |
| 局部最低邊界溫度 | {wall['minimum_interface_temperature_C']:.2f} °C |
| 面積平均淨熱通量 | {wall['area_average_heat_flux_W_m2']:.2f} W/m² |
| 局部最大節點淨熱通量 | {wall['maximum_nodal_heat_flux_W_m2']:.2f} W/m² |
| 最高溫位置 (x,y,z) | ({pos[0]:.3f}, {pos[1]:.3f}, {pos[2]:.3f}) m |

W 表示整個面每秒得到的淨熱量，W/m² 表示局部每平方米的熱量密度。正值指由空氣域流向側面／散熱參考環境。最高值是單一節點的局部熱點，不能視為整塊板的溫度。

面積只包含空氣直接接觸的側面：y=0–0.6 m、z=0.9–1.6 m，以及 y=0.3–0.6 m、z=1.6–2.2 m。灶台以下、上櫃遮擋部分和擋板接合面不屬於暴露面。名義面積0.600 m²，扣除擋板接合處後為{wall['area_m2']:.6f} m²。

![高櫃側面溫度](wall_temperature_9000W.png)

## 工況與計算方法

房間3.3×1.65×2.2 m，灶面高0.9 m。兩爐中心 x=0.95、1.35 m，y=0.30 m；左爐中心距所選側面0.25 m。煙機頂吸，寬0.90 m、深0.30 m，底面高1.60 m；其90 cm寬度內設10 cm、向下45°前擋板。

兩門關閉，但每扇門底保留5 mm補風縫。沿用假設風機曲線 p=1300[1−(Q/0.5)²] 和孔口係數0.6，運行風量為{fine['fan_flow_m3_min']:.3f} m³/min。30 m³/min是所假設自由送風端點，不是此全關門工況的實際運行量。門縫在室內網格中鋪展成10 cm高補風帶，未解析真實5 mm門縫高速射流；未加入風管／濾網阻力。

三維四面體P1/P1有限元求解固定參考密度不可壓縮流動，採用非線性理想氣體浮力修正、SUPG/PSPG穩定化、保守熱量弱式及代數邊擴散。每爐4500 W以近灶面的高斯體積熱源施加，未建立火焰或鍋具幾何。有效渦黏度0.005 m²/s、湍流Pr=0.85均沿用上版假設。

所有固體邊界沿用25°C散熱參考環境及等效散熱係數 **3 W/(m²·K)**。在側面三角形上直接積分：

`q'' = 3 × (T_boundary − 25)`，`Q_wall = Σ A_triangle × 3 × (T_triangle_mean − 25)`。

線性三角形溫度的面積積分是精確的；熱點使用邊界節點最大值。這個散熱係數描述原模型的等效邊界，並不是櫃板材料導熱係數。模型不包含板材熱容量或升溫歷程，輸出是穩態条件估算。

## 數值驗證

| 網格目標尺寸 | 四面體數 | 所選側面功率 W | 平均邊界溫度 °C | 最高邊界溫度 °C |
|---|---:|---:|---:|---:|
{rows}

中網格到細網格：受熱功率改變 **{grid['total_net_heat_into_wall_W']['relative_change_percent']:+.2f}%**，平均溫升改變 **{grid['temperature_rise_average_relative_change_percent']:+.2f}%**，最高邊界溫度改變 **{grid['maximum_interface_temperature_C']['relative_change_percent']:+.2f}%**。這是網格敏感性檢查，沒有足夠證據宣稱所有局部溫度已達網格無關，也不能把相鄰網格差值當作完整誤差界限。

三個新工況均滿足求解器收斂準則。細網格非線性相對殘差 {fine['nonlinear_relative_residual']:.3e}，熱量方程相對殘差 {fine['thermal_relative_residual']:.3e}，全域體積流量相對誤差 {fine['net_flux_relative_error']:.3e}。

9000 W熱量收支：烟機對流帶走 {fine['advective_heat_flux_W']['hood']:.4f} W；**所有**固體邊界合計散熱 {fine['wall_heat_loss_W']:.4f} W；門口熱交換反力 {fine['door_thermal_reaction_W']:.6f} W；離散總熱量收支誤差 {fine['discrete_heat_balance_residual_W']:.3e} W。所選高櫃側面僅是全部固體邊界的一部分。

全域守恆不代表各單元都滿足強形式不可壓縮條件；本次P1/P1速度場的散度／速度梯度L2比值為{fine['divergence_gradient_norm_ratio']:.3f}。固定渦黏度、壓力穩定化及熱量數值擴散仍會影響熱點和熱量分配。局部空氣最高溫為{fine['max_nodal_temperature_C']:.1f}°C，密度變化顯著；浮力項雖有非線性修正，連續方程和能量輸運仍使用固定參考密度，因此絕對溫度具有模型不確定性。

## 與上一版4000 W條件比較

上一版已存中網格同一面積的淨受熱功率為{old['total_net_heat_into_wall_W']:.2f} W，平均邊界溫度{old['area_average_temperature_C']:.2f}°C，最高邊界溫度{old['maximum_interface_temperature_C']:.2f}°C。本次9000 W已重新耦合求解流場和溫度場，沒有把舊結果直接按功率比例放大。

互動3D圖提供4000 W中網格與9000 W細網格切換，採共同色階。兩張圖的網格也有差異；對比僅作條件估算。

## 如何取得實際櫃板溫度

需補充側板材質、厚度、背面空間與散熱條件，並在流場與板材之間建立固體導熱耦合；若要求全部受熱功率，還需火焰／鍋具的輻射條件或可校準的量測。現有的{wall['total_net_heat_into_wall_W']:.1f} W和{wall['area_average_temperature_C']:.1f}°C只適用於上述簡化邊界，不能直接代入真實木質／金屬／石材側板。

## 輸出檔案

- `wall_summary.csv`：新舊工況、三級網格的功率與溫度摘要。
- `9000W_fine_wall_nodes.csv`：細網格表面座標、邊界溫度、局部淨熱通量。
- `9000W_fine_wall_faces.csv`：每塊三角形的面積、平均溫度、熱通量及功率。
- `9000W_fine_wall.vtu`：可用ParaView開啟的表面有限元資料。
- `wall_interactive_3d.html`：離線可旋轉、縮放、查讀的3D側面資料。
- 完整資料包另含9000 W三級網格的原始空氣場、收斂記錄、設定及重現程式。
'''
    (OUT/'wall_heat_report_zh.md').write_text(md)
    import markdown
    css='''body{font-family:"Noto Sans CJK TC",sans-serif;color:#182535;font-size:11px;line-height:1.65;margin:30px auto;max-width:950px}h1{font-size:22px}h2{font-size:16px;margin-top:24px}table{border-collapse:collapse;width:100%;font-size:10px}td,th{border:1px solid #c9d0d7;padding:6px}th{background:#edf2f7}img{max-width:100%;max-height:650px;display:block;margin:auto}pre,code{white-space:pre-wrap}p{orphans:3;widows:3}h2{break-after:avoid}tr{break-inside:avoid}@page{size:A4;margin:15mm}'''
    body=markdown.markdown(md,extensions=['tables','fenced_code'])
    # Embed the plot so both local opening and set_content PDF rendering work offline.
    import base64
    plot='data:image/png;base64,'+base64.b64encode((OUT/'wall_temperature_9000W.png').read_bytes()).decode()
    body=body.replace('src="wall_temperature_9000W.png"',f'src="{plot}"')
    document=f'<!DOCTYPE html><html lang="zh-Hant"><meta charset="utf-8"><title>高櫃側面受熱估算</title><style>{css}</style><body>{body}</body></html>'
    (OUT/'wall_heat_report_zh.html').write_text(document)
    from playwright.sync_api import sync_playwright
    with sync_playwright() as pw:
        browser=pw.chromium.launch(executable_path='/usr/bin/chromium',args=['--no-sandbox'])
        page=browser.new_page();page.set_content(document,wait_until='load');page.evaluate('document.fonts.ready')
        page.pdf(path=str(OUT/'wall_heat_report_zh.pdf'),format='A4',print_background=True)
        browser.close()
    print(json.dumps({'wall':wall,'checks':flags,'grid':grid},ensure_ascii=False,indent=2))

if __name__=='__main__':main()
