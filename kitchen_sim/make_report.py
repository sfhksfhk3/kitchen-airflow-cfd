"""Build the Chinese report from the actual solver metrics; never invent results."""
from pathlib import Path
import json
import shutil
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parent
CFG = json.loads((ROOT / 'config.json').read_text())


def table(headers, rows):
    return '\n'.join(['| ' + ' | '.join(headers) + ' |',
                      '| ' + ' | '.join(['---'] * len(headers)) + ' |'] +
                     ['| ' + ' | '.join(map(str, row)) + ' |' for row in rows])


def num(x, digits=3):
    return f'{x:.{digits}f}'


def main():
    metrics = json.loads((ROOT / 'results/metrics.json').read_text())
    cases = CFG['cases']
    missing = [c for c in cases if c not in metrics]
    if missing:
        raise SystemExit('Report waits for completed metric files: ' + ', '.join(missing))
    rows = [metrics[c] for c in cases]
    if any(m.get('heat_formulation')!='conservative_galerkin_edge_diffusion' for m in rows):raise SystemExit('Final conservative thermal results required')
    labels = CFG['case_labels']
    h = rows[0]['mesh_target_m']
    extras = [(k, v) for k, v in metrics.items() if k not in cases and 'initial' not in k]
    result_table = table(
        ['情境', '風量 m³/min', '補風網路壓差 Pa', '走道平均速度 m/s', '走道平均溫度 °C', '最高節點溫度 °C'],
        [[labels[c], num(metrics[c]['fan_flow_m3_min'], 2), num(metrics[c]['network_makeup_pressure_drop_Pa'], 3),
          num(metrics[c]['occupied_zone_mean_speed_m_s']), num(metrics[c]['occupied_zone_mean_temperature_C'], 2),
          num(metrics[c]['max_nodal_temperature_C'], 2)] for c in cases])
    door_table = table(
        ['情境', '左門淨流入 m³/min', '右門淨流入 m³/min', '煙機流出 m³/min', '全邊界淨不平衡 / 抽風量'],
        [[labels[c], num(-metrics[c]['boundary_flux_m3_s']['left']['net_outward'] * 60),
          num(-metrics[c]['boundary_flux_m3_s']['right']['net_outward'] * 60),
          num(metrics[c]['boundary_flux_m3_s']['hood']['net_outward'] * 60),
          f"{metrics[c]['net_flux_relative_error'] * 100:.5f}%"] for c in cases])
    numerical_table = table(
        ['情境', '節點 / 四面體', '迭代次數', '更新準則達標', '速度相對更新', '溫度相對更新', '穩態動量等相對殘差'],
        [[labels[c], f"{metrics[c]['nodes']:,} / {metrics[c]['tetrahedra']:,}", metrics[c]['iterations'],
          '是' if metrics[c]['converged'] else '否', f"{metrics[c]['last_relative_velocity_update']:.2e}",
          f"{metrics[c]['last_relative_temperature_update']:.2e}",
          f"{metrics[c]['nonlinear_relative_residual']:.2e}"] for c in cases])
    heat_table = table(
        ['情境', '煙機顯熱流出 W', '牆壁散熱 W', '門口熱邊界反力 W', '離散熱平衡餘差 W'],
        [[labels[c], num(metrics[c]['advective_heat_flux_W']['hood'], 2), num(metrics[c]['wall_heat_loss_W'], 2),
          num(metrics[c]['door_thermal_reaction_W'], 2),
          num(metrics[c]['discrete_heat_balance_residual_W'], 4)] for c in cases])
    all_converged = all(m['converged'] for m in rows)
    status = ('四種基準情境均達到程式設定的迭代更新準則。' if all_converged else
              '部分基準情境未達到程式設定的迭代更新準則；相應數值是停止迭代時的近似場，不能當作已收斂的穩態結果。')
    base = metrics['both_open']
    refinement_rows = []
    closure_rows = []
    for key, m in extras:
        if m.get('case') != 'both_open':
            continue
        item = [key, num(m['mesh_target_m']), num(m['nu_t_m2_s'], 4),
                num(m['occupied_zone_mean_speed_m_s']), num(m['occupied_zone_mean_temperature_C'], 2),
                num(m['max_nodal_temperature_C'], 2), '是' if m['converged'] else '否',
                f"{m['nonlinear_relative_residual']:.2e}"]
        if abs(m['nu_t_m2_s'] - base['nu_t_m2_s']) < 1e-10 and abs(m['mesh_target_m'] - base['mesh_target_m']) > 1e-10:
            refinement_rows.append(item)
        elif abs(m['nu_t_m2_s'] - base['nu_t_m2_s']) > 1e-10:
            closure_rows.append(item)
    validation_headers = ['資料集', 'h m', 'νt m²/s', '走道均速 m/s', '走道均溫 °C', '最高節點溫度 °C', '更新準則達標', '穩態相對殘差']
    base_item = ['both_open（基準）', num(h), num(base['nu_t_m2_s'], 4), num(base['occupied_zone_mean_speed_m_s']),
                 num(base['occupied_zone_mean_temperature_C'], 2), num(base['max_nodal_temperature_C'], 2),
                 '是' if base['converged'] else '否', f"{base['nonlinear_relative_residual']:.2e}"]
    refinement_text = (table(validation_headers, [base_item] + refinement_rows) if refinement_rows else
                       '本次未取得可比較的加密網格結果，不能聲稱已完成網格獨立性分析。')
    closure_text = (table(validation_headers, [base_item] + closure_rows) if closure_rows else
                    '本次未取得固定渦黏性變動的比較結果；此模型參數造成的不確定性未量化。')
    relative_rows = []
    for key, m in extras:
        if m.get('case') == 'both_open' and (abs(m['mesh_target_m'] - h) > 1e-10 or abs(m['nu_t_m2_s'] - base['nu_t_m2_s']) > 1e-10):
            relative_rows.append([key,
                f"{(m['occupied_zone_mean_speed_m_s'] / base['occupied_zone_mean_speed_m_s'] - 1) * 100:+.2f}%",
                f"{m['occupied_zone_mean_temperature_C'] - base['occupied_zone_mean_temperature_C']:+.3f}",
                f"{m['max_nodal_temperature_C'] - base['max_nodal_temperature_C']:+.3f}"])
    differences = table(['資料集對基準差異', '走道均速相對差', '走道均溫差 °C', '最高節點溫度差 °C'], relative_rows) if relative_rows else ''
    mesh_note=''
    if 'both_open_refined' in metrics:
        refined=metrics['both_open_refined']
        speed_change=(refined['occupied_zone_mean_speed_m_s']/base['occupied_zone_mean_speed_m_s']-1)*100
        peak_change=refined['max_nodal_temperature_C']-base['max_nodal_temperature_C']
        mesh_note=f'本次加密後走道平均風速改變 {speed_change:+.2f}%，最高節點溫度改變 {peak_change:+.2f} °C。局部峰值尚未達網格無關；基準工況間幾度的溫差不宜用來排序捕煙效果。'

    image_candidates = [('model_3d.png', '三維幾何與計算流線'),
                        ('four_cases_velocity.png', '四種工況速度截面'), ('four_cases_temperature.png', '四種工況溫度截面'),
                        ('figures/convergence.png', '迭代收斂曲線'), ('figures/mesh_comparison.png', '兩門全開的網格比較')]
    figures = '\n\n'.join(f'![{caption}]({path})' for path, caption in image_candidates if (ROOT / path).exists())
    closed = metrics['both_closed']
    max_flux_error = max(m['net_flux_relative_error'] for m in rows)
    max_div_ratio = max(m['divergence_gradient_norm_ratio'] for m in rows)
    run_commands = [f'python solve.py --case all --h {h:g} --maxiter {max(250, max(m["iterations"] for m in rows))}']
    for key, m in extras:
        if m.get('case') not in cases:
            continue
        suffix = key[len(m['case']):]
        run_commands.append(f'python solve.py --case {m["case"]} --h {m["mesh_target_m"]:g} --nu-t {m["nu_t_m2_s"]:g} --suffix {suffix} --maxiter {max(250, m["iterations"])}')
    run_commands += ['python validate_exports.py', 'python visualize.py', 'python make_mesh_comparison.py', 'python make_report.py']
    command_text = '\n'.join(run_commands)
    text = f"""# 廚房抽油煙機氣流：三維有限元素初步分析

**適用範圍：**依手繪尺寸建立的條件式工程模型，用於比較兩門全開、全關、左开右關及右開左關。計算包含三維速度、室內相對壓力與空氣溫度；未模擬油滴、燃燒反應、油煙濃度或污染物捕集率。報告與圖中的數據來自 `solve.py` 的有限元素計算，風量則由另列的假設風機／補風阻力模型給定。

**主要結果：**在本文的二次風機曲線、補風流量係數 0.6、無額外風管／濾網阻力的假設下，一扇或兩扇門全開時的風量約 30.00 m³/min；兩門全關、每門底縫 5 mm 時估計為 **{closed['fan_flow_m3_min']:.2f} m³/min**，約給定 30 m³/min 的 {closed['fan_flow_m3_s'] / .5 * 100:.1f}%。全關結果對實際門縫與風機曲線非常敏感。這個數字不是以 1300 Pa 和 30 m³/min 同時運轉所得，也不是量測值。

{status}數值收斂、網格敏感性與實物驗證是不同的檢查；本次沒有實物量測驗證。

## 1. 幾何、座標與已知條件

x 軸沿廚房 3.3 m 長度由左向右；y 軸由後方牆面向走道；z 軸自地面向上。以下尺寸以公尺為單位。

| 項目 | 模型採用值 | 來源／說明 |
| --- | --- | --- |
| 室內長 × 深 × 高 | 3.30 × 1.65 × 2.20 | 長、深由草圖；2.20 m 作為室高解讀 |
| 室內毛體積 | 11.979 m³ | 含櫥櫃占用空間 |
| 左端高櫃 | x = 0–0.60、y = 0–0.60、z = 0–2.20 | 左端櫃寬 0.60 m 為估計 |
| 下櫃／檯面 | x = 0.60–3.30、y = 0–0.60、z = 0–0.90 | 檯深 0.60 m、檯高 0.90 m |
| 上櫃／煙機占用區 | x = 0.60–3.30、y = 0–0.60、z = 1.60–2.20 | 簡化為實心區域，不建內部風道 |
| 煙機吸入口 | x = 0.60–1.50、y = 0–0.60、z = 1.60 | 頂吸式、距檯面 0.70 m 已確認；0.90 m 寬，深 0.60 m 為假設 |
| 雙爐中心 | (0.85, 0.30, 0.90)、(1.25, 0.30, 0.90) | 中心距 0.40 m；左爐離高櫃邊 0.25 m |
| 左／右门 | x = 0／3.30；y = 0.70–1.60；高 2.20 | 兩門寬 0.90 m、高 2.20 m 為使用者確認；沿 y 位置為估計 |
| 關門門底縫 | 每扇 0.90 × 0.005 m | 實際未知，依使用者同意暫採每門 5 mm |
| 網格流體體積 | {base['fluid_volume_m3']:.3f} m³ | 扣除模型中實心櫃體後的積分體積 |

門片在全開時移出計算域，沒有模擬門片開啟角度或阻擋。牆厚、把手、鍋具、人體、洗碗槽、濾網細節、煙機內部、排煙管與室外環境均未建模。幾何以 `config.json` 與 `solve.py` 為準；調整尺寸時須同步調整網格／邊界定義，不能只修改設定檔就假設程式會自動重建所有邊界。

## 2. 風機與補風：30 m³/min、1300 Pa 如何使用

30 m³/min = 0.50 m³/s 是給定風量，1300 Pa 是給定最大靜壓；尚未確認 30 m³/min 是否為實機的自由出風端點。風量與最大靜壓通常對應不同工作點，**不能直接同時指定為同一工況的邊界條件**。本次為了形成可計算的示例，暫把它們視作風機曲線的自由出風／零風量兩端點：

```
Δp_fan(Q) = 1300 × [1 − (Q / 0.50)²]                  [Pa]
Δp_makeup(Q) = ρ / 2 × [Q / (Cd × A)]²               [Pa]
ρ = 1.184 kg/m³；Cd = 0.60
Q = sqrt(1300 / [1300 / 0.50² + ρ / (2 Cd² A²)])     [m³/s]
```

此模型只將門口阻力與風機曲線配對。未計風管、彎頭、止回閥、濾網、公共煙道及室外風壓，因此開門情境的約 30 m³/min 接近假設的自由出風量，不能視為現場保證。FEM 不重新求取整台風機工作點；它採用上述 Q，計算室內流場。

| 情境 | 阻力網路採用的物理補風面積 A |
| --- | --- |
| 兩門全開 | 2 × 0.90 × 2.20 = 3.960 m² |
| 左開右關、右開左關 | 0.90 × 2.20 = 1.980 m²；忽略另一扇關門底縫 |
| 兩門全關 | 2 × 0.90 × 0.005 = 0.009 m² |

全關时，在假設孔口係數下，即使把全部 1300 Pa 用於門縫，理論補風上限也僅約 **0.253 m³/s = 15.18 m³/min**；在二次曲線交點則為 {closed['fan_flow_m3_min']:.2f} m³/min、門縫壓差 {closed['network_makeup_pressure_drop_Pa']:.1f} Pa。若硬要求 30 m³/min 穿過這兩條門縫，僅門縫壓降就約 5075 Pa，已高於給定最大靜壓。

全關工況的實際門縫面積平均速度由網路估為約 {closed['fan_flow_m3_s'] / .009:.1f} m/s，**這不是室內 FEM 已解析的局部速度**。5 mm 細縫沒有做細網格解析：為計算遠離門縫的室內流動，FEM 將每側入口鋪展為 0.10 m 高的等效帶，但風量仍以真正 5 mm 縫面積計算。鋪展會改變局部噴流動量與混合，全關流場只能作遠場估算，不能從圖上判斷門縫噴流速度、噪音或精確舒適度。

若兩門及其他開口完全密封，穩態不可壓縮質量守恆要求 Q = 0。這只是完全密封的物理極限說明；本次四種計算不包含「完全密封仍維持抽風」的虛構工況，也未計啟動瞬間的壓力變化。

## 3. 物理與數值模型

模型求解三維穩態不可壓縮動量方程與熱輸運方程，使用 Boussinesq 浮力近似：

```
∇·u = 0
(u·∇)u = −∇p / ρ + ν_eff ∇²u + g β (T − T₀) e_z
∇·(u θ) = α_eff ∇²θ + q''' / (ρ cp)；θ = T − T₀
ν_eff = ν + νt；α_eff = α + νt / Prt
```

三维四面體有限元素使用 `scikit-fem` 的 P1 速度／P1 壓力／P1 溫度離散。動量／連續性方程採 SUPG／PSPG 穩定化與 grad-div 項，grad-div 係數為 0.5 m²/s。自然應力門口增加 0.5 × max(−u·n, 0) 的速度邊界耗散，以抑制回流引起的能量不穩定。溫度使用保守 Galerkin 弱式：對流項為 −∫θ u·∇v，加上邊界 ∫(u·n)θv，因此全域熱平衡不需要以局部速度散度作事後修正。另加對稱代數邊擴散，使非對角項非正、抑制振盪，並實際檢查最低節點溫度；附加擴散矩陣的行／列和為零，會增加數值混合、平滑溫度梯度。這是數值穩定化，不是額外的物理湍流模型。

採耦合固定點迭代、鬆弛與動量偽時間推进以尋找穩態；迭代步數不代表真實流動經過的秒數。偽時間項在固定點抵消，最終另以穩態離散算子檢查殘差。四個基準工況承接先前已收斂流場，表內迭代次數是最後保守熱輸運修正後的迭代；冷啟動重算所需次數可能較多。基準網格的分段座標目標間距 h = {h:g} m；h 不是所有四面體最大邊長。門底等效帶額外包含 z = 0.05 m 的節點層，但仍未解析真實 5 mm 縫。

| 物理量 | 採用值／處理 |
| --- | --- |
| 環境溫度 T₀ | 25 °C |
| 空氣密度 ρ、比熱 cp | 1.184 kg/m³、1006 J/(kg·K) |
| 重力 g、熱膨脹係數 β | 9.81 m/s²、0.003354 K⁻¹ |
| 分子運動黏度 ν、熱擴散率 α | 1.56 × 10⁻⁵ m²/s、2.20 × 10⁻⁵ m²/s |
| 固定渦黏度 νt | 0.005 m²/s，為未校準的混合強度假設 |
| 渦流普朗特數 Prt | 0.85 |
| 雙爐顯熱 | 使用者確認兩個燃氣爐同時使用，功率未知；每爐向空氣釋放 600 W，共 1200 W 為假設，並非已知燃氣輸入功率 |
| 熱源形狀 | 爐中心上方的平滑體積熱源，水平尺度 0.085 m、垂直尺度 0.055 m，分別積分歸一至 600 W |
| 固體表面 | 無滑移；以 3 W/(m²·K) 的 Robin 邊界向 25 °C 環境散熱 |
| 補風開口 | 自然應力條件加回流能量穩定化；溫度設 25 °C，包含發生局部流出的開口節點 |
| 煙機出口 | 指定向上流速，正弦分布在邊緣降至零；積分後恰為指定風量 Q；溫度使用自然擴散邊界 |

固定渦黏度是簡化的有效混合模型，未求解 k–ε、k–ω、LES 等湍流模型；沒有火焰、鍋具傳熱、蒸發與輻射交換。最高節點溫度是經平滑後熱源附近的空氣數值，不能解讀為火焰、鍋面或人體暴露溫度。高溫升處也會降低 Boussinesq 近似的準確性。

## 4. 四種情境的計算數據

以下平均值以有限元素積分權重計算。「走道區」定義為 y > 0.60 m 且 0.90 < z < 1.80 m 的流體區域，不是全房間平均，也不代表任一人體位置的局部風速。

{result_table}

煙機總面積假定為 0.90 × 0.60 = 0.54 m²；30 m³/min 對應的面積平均吸入速度約 0.926 m/s。模型使用非均勻吸入口分布，因此節點峰值速度可能明顯高於這個面平均值。

### 門口流量與質量守恆

{door_table}

正值表示門口淨補風，负值表示門口淨流出。`metrics.json` 另外提供總流入、總流出，供識別同一門口不同區域同時流入與流出的現象。上述全域流量差很小也不等於單元內處處滿足 ∇·u = 0；本計算採弱不可壓縮穩定化離散，應併看局部散度指標。

本模型中兩門全開的補風分配為左門約 {(-metrics['both_open']['boundary_flux_m3_s']['left']['net_outward']*60):.2f} m³/min、右門約 {(-metrics['both_open']['boundary_flux_m3_s']['right']['net_outward']*60):.2f} m³/min，與煙機靠左的非對稱幾何相符。只開右門時的走道區均速約 {metrics['right_open']['occupied_zone_mean_speed_m_s']:.3f} m/s，高於只開左門的 {metrics['left_open']['occupied_zone_mean_speed_m_s']:.3f} m/s，反映補風穿越較長的室內路徑。這是流場路徑比較，沒有據此推定捕煙率或哪側開門最佳。

## 5. 收斂、守恆與敏感性

收斂旗標的設定是速度與溫度相對更新量都小於 2 × 10⁻⁵，並連續達標 3 次。另行報告重組穩態算子的殘差，避免只靠變化量很小判斷求解品質。速度／溫度更新量與方程殘差的分母不同，數值不宜直接互比。

{numerical_table}

四種基準情境最大相對全域流量不平衡為 {max_flux_error * 100:.5f}%；最大 ‖∇·u‖₂ / ‖∇u‖₂ 為 {max_div_ratio:.4f}。後者是局部弱不可壓縮誤差的整體指標，不是體積百分比。四種情境的最低節點溫度最小值為 {min(m['min_nodal_temperature_C'] for m in rows):.6f} °C，可與 25 °C 入口條件比較以檢查非物理低溫振盪。局部散度指標仍偏大，細部流速與熱羽流需要更高階離散及進一步加密確認；極小的全域風量差不能視為局部流場已達高精度。線性求解器殘差與熱方程殘差請見每種情境的 `_metrics.json`；完整更新歷史見 `_convergence.csv`。

### 離散热平衡

{heat_table}

離散檢查使用「開口顯熱流出 + 牆散熱 + 門口熱邊界反力 − 1200 W」。溫度方程採保守形式，沒有事後扣除速度散度造成的熱差額。门口温度被固定为环境温度，热边界反力包括相应扩散／稳定化交换。將局部流出門口也固定為 25 °C 可能額外冷卻流場，因此門口熱反力包含這項建模影響，不能直接當作現場真實門口傳熱量。此表是對離散方程的核對，**不是油煙捕集率**。

### 網格加密檢查（兩門全開）

{refinement_text}

{mesh_note}

只有一組基準／加密比較時，不能宣稱完成網格獨立性或給出網格收斂指數。應將平均量與熱源峰值分開解讀；平滑局部熱源仍可能對网格解析度敏感。若任一比較工況未達穩態準則，差異同時包含迭代誤差。

### 固定渦黏度敏感性（兩門全開）

{closure_text}

{differences}

改變 νt 同時改變動量與熱量的有效擴散；這是模型假設敏感性，不能據此得到统计置信区间或已校準的湍流預測。仅对两门全开检查，不能替代另外三种情境的网格与模型检验。

## 6. 壓力、流線與圖像如何閱讀

FEM 匯出的 `pressure_room_reference_Pa` 以**等效補風開口的自然應力條件**作參考，表示已建模室內流域的相對機械壓力；未含完整靜水壓背景。`network_makeup_pressure_drop_Pa` 是獨立的門／縫阻力網路壓降。尤其全關情境，圖上室內相對壓力不會自動包含約千帕的門縫壓降；兩者不能混用，也不能把兩幅不同參考的數值直接相減。

流線由穩態速度場積分得到，只展示示蹤軌跡的幾何方向。它們不是實際油滴軌跡，不帶經過時間，不表示煙氣質量，也不能依「畫出的線有多少進入煙機」換算捕集率。色階以各圖標示的物理單位為準；比較情境時應採相同色階。圖上門底等效帶不可誤認為真實 10 cm 門縫。

{figures}

## 7. 可作出的判斷與尚不能作出的判斷

在指定補風面積與假設风机曲线下，兩門全關會大幅降低可用風量；至少一扇門全開時，門口阻力相對很小。一扇門全開與兩門全開在本模型中都接近同一總抽風量，差別主要來自補風路徑和室內流場分布。哪側開門較佳，需要結合三維流場、實際煙機型式、使用者站位及熱羽流外逸檢查，不能僅依房間平均風速排序。

這組模型不能證明實際捕煙率、除味率、燃燒安全、噪音、每一處舒適度或實機的最大能力。若要作設備／施工決策，最有用的補充資料是實機整條風量—靜壓曲線與風檔、風管／濾網阻力、實際煙機吸入口與櫥櫃尺寸、爐火向室內的顯熱，以及門縫／其他漏風通道尺寸。之後再建未鋪展門縫、外部補風域、瞬態湍流與污染物輸運，並以現場風量、壓力與溫度量測校準。

## 8. 檔案與數據字典

| 檔案 | 內容／用途 |
| --- | --- |
| `config.json` | 幾何、物性、風機與假設參數；求解器讀取幾何參數；更動尺寸時視覺化外形也應同步更新 |
| `solve.py` | 真正的三維有限元素求解、原始場資料與守恆指標匯出 |
| `results/metrics.json` | 四個基準情境與完成的網格比較彙總；最完整數據來源 |
| `results/<case>_metrics.json` | 單一情境摘要与数值质量指标 |
| `results/<case>_nodes.csv` | 所有網格節點的座標、三分量速度、壓力與溫度 |
| `results/<case>.vtu` | 四面體網格與三維節點場，可在 ParaView 中旋轉、切片與產生流線 |
| `results/<case>.npz` | NumPy 壓縮原始網格與場；供本附檔可視化使用 |
| `results/<case>_convergence.csv` | 每次迭代的更新量、線性殘差、溫度極值及執行時間 |
| `interactive_3d.html`、`visualize.py` | 離線互動 3D、原始四面體插值切片與流線；可切換四種情境 |
| `geometry.obj`、`probes.csv` | 櫥櫃／煙機幾何和六個固定位置的探針數據 |
| `validate_exports.py`、`results/validation.json` | 從匯出網格獨立重算體積、流量、速度梯度並核對 CSV／VTU／NPZ |
| `make_report.py`、`REPORT.md`、`report.html`、`report.pdf` | 依實際彙總數據產生的本報告；HTML 可由瀏覽器列印，PDF 為已排版版本 |

`<case>` 為 `both_open`、`both_closed`、`left_open` 或 `right_open`；其他後綴是敏感性試算。

| 欄位／變數 | 單位與定義 |
| --- | --- |
| `x_m, y_m, z_m` / `points` | m；節點位置 |
| `ux_m_s, uy_m_s, uz_m_s` / `velocity` | m/s；x、y、z 方向速度；正 z 向上 |
| `pressure_room_reference_Pa` / `pressure` | Pa；FEM 室內自然應力參考壓力，與門縫網路壓降分開 |
| `temperature_C` / `temperature` | °C；節點空氣溫度 |
| `cells` | 每個四面體的 4 個節點索引（0 起算） |
| `boundary_hood, boundary_left, boundary_right` | scikit-fem mesh facet 索引；不是節點索引，若重建需保留原網格拓撲 |
| `fan_flow_m3_s`, `fan_flow_m3_min` | m³/s、m³/min；假設風機／補風網路的工作風量 |
| `network_makeup_pressure_drop_Pa` | Pa；物理門口或 5 mm 門縫壓差估計 |
| `boundary_flux_m3_s.*.net_outward` | m³/s；積分法向流量，向計算域外為正 |
| `gross_inward`, `gross_outward` | m³/s；同一邊界的總流入／總流出絕對量，皆非負 |
| `mean_*`, `occupied_zone_mean_*` | 體積積分平均；後者只限本文定義的走道區 |
| `max_nodal_*`, `min_nodal_*` | 節點最大／最小值；可能隨網格改變 |
| `converged` | 僅表示達到已列的連續更新準則，並非實物驗證標記 |
| `nonlinear_relative_residual` | 重組穩態動量／連續性離散方程之相對残差 |
| `thermal_relative_residual` | 自由溫度節點的離散方程相對殘差 |
| `divergence_gradient_norm_ratio` | ‖∇·u‖₂ / ‖∇u‖₂；無量綱的弱不可壓縮誤差指標 |
| `*_heat_*_W`, `door_thermal_reaction_W` | W；顯熱傳輸、邊界反力與離散守恆檢查量 |
| `runtime_s`, `elapsed_s` | 電腦計算耗時；不是物理模擬時間 |

獨立匯出核對見 `results/validation.json`：以原始四面體重新計算體積、法向流量和速度梯度，並核對 CSV、VTU 與 NPZ 的一致性。這項核對不重用求解器的有限元素組裝。

## 9. 重現計算

使用 Python 3，所需核心套件為 `numpy scipy scikit-fem meshio`；生成图表需要 `matplotlib`。在本專案資料夾執行以下命令。若已有 `requirements.txt`，優先依該檔建立相同環境。生成可列印 HTML 需要系統有 Pandoc，PDF 另需 Chromium 和 Playwright；無對應工具時仍會生成 Markdown。

```bash
python -m pip install -r requirements.txt
mkdir -p results
{command_text}
```

重新求解會覆寫同名情境的結果。CSV、VTU 與 NPZ 都使用 SI 基本單位（溫度另為 °C），沒有以繪圖色階值替代物理資料。更動幾何、熱源、風機或邊界條件後，必須重跑求解與可視化；不能只改圖形標籤。
"""
    # Keep this deliverable consistently in Traditional Chinese.
    replacements = {'左开': '左開', '右门': '右門', '全關时': '全關時', '三维': '三維', '推进': '推進',
                    '负值': '負值', '離散热': '離散熱', '门口温度被固定为环境温度，热边界反力包括相应扩散／稳定化交换。':
                    '門口溫度被固定為環境溫度，熱邊界反力包括相應物理／代數邊擴散交換。',
                    '网格': '網格', '统计置信区间': '統計置信區間',
                    '仅对两门全开检查，不能替代另外三种情境的网格与模型检验。':
                    '僅對兩門全開檢查，不能替代另外三種情境的網格與模型檢驗。',
                    '风机曲线': '風機曲線', '与数值质量指标': '與數值品質指標',
                    '残差': '殘差', '生成图表': '生成圖表'}
    for a, b in replacements.items():
        text = text.replace(a, b)
    (ROOT / 'REPORT.md').write_text(text)
    css = '''body{font-family:"Noto Sans CJK TC","Noto Sans CJK SC",sans-serif;max-width:1100px;margin:36px auto;padding:0 30px;color:#243247;line-height:1.65;font-size:14px}h1{font-size:29px;color:#10243b}h2{font-size:22px;border-bottom:1px solid #ccd8e5;padding-bottom:5px;margin-top:38px}h3{font-size:17px;margin-top:26px}table{width:100%;border-collapse:collapse;margin:16px 0;font-size:12px}th,td{border:1px solid #d2dce6;padding:7px 9px;text-align:left;vertical-align:top}th{background:#eaf1f7}tr:nth-child(even){background:#f8fafc}code{font-family:monospace;font-size:.9em;overflow-wrap:anywhere}pre{background:#f0f4f8;padding:16px;white-space:pre-wrap;overflow-wrap:anywhere}img{max-width:100%}strong{color:#123453}a{color:#065e85}@media print{body{max-width:none;margin:0;padding:0;font-size:10pt}h2,h3{break-after:avoid}tr,img{break-inside:avoid}table{font-size:8pt}a{color:inherit;text-decoration:none}pre{font-size:8pt}}@page{size:A4;margin:17mm}'''
    (ROOT / 'report.css').write_text(css)
    try:
        subprocess.run(['pandoc', str(ROOT / 'REPORT.md'), '--standalone', '--embed-resources', '--metadata', 'lang=zh-Hant',
                        '--metadata', 'pagetitle=廚房抽油煙機氣流：三維有限元素初步分析',
                        '--css', 'report.css', '-o', str(ROOT / 'report.html')], check=True, cwd=ROOT)
    except FileNotFoundError:
        print('Pandoc is unavailable; REPORT.md was generated.')
    else:
        chromium = shutil.which('chromium')
        if chromium:
            from playwright.sync_api import sync_playwright
            with sync_playwright() as playwright:
                browser=playwright.chromium.launch(executable_path=chromium,headless=True,args=['--no-sandbox'])
                page=browser.new_page()
                page.set_content((ROOT/'report.html').read_text(),wait_until='load')
                page.pdf(path=str(ROOT/'report.pdf'),format='A4',print_background=True,prefer_css_page_size=True)
                browser.close()
    print('Report generated from', len(rows), 'main cases and', len(extras), 'supplemental metric datasets.')


if __name__ == '__main__':
    main()
