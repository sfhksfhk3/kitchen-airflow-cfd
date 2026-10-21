# 沉積物理、方程與資料依據

## 保守方程和逐面質量

每個粒徑分別求解：∂(ρY)/∂t + ∇·[(ρU+ρv_d)Y] − ∇·[ρ(D_B+ν_t/Sc_t)∇Y] = S − Σ_face[ρv_dep A/V]Y。

Y 為無因次質量分率；原生 phi 為 kg/s；D_B、ν_t 為 m²/s；v_dep 為 m/s。壁面不是以靠近天花濃度代替沉積：在離散方程以隱式 sink 移除的同一份質量，按面分配 J_face=ρYv_dep（kg/(m²s)），M_face←M_face+Δt J_face（kg/m²）。輸出 mg/m² 時乘10⁶。天花總沉積=Σ M_face A_face。所有牆面 sink 加總與油滴收支逐步核對；開口採矩陣的數值通量，避免另估邊界通量不守恆。

## 壁面近似

D(y)=D_B+κu_τ y[1−exp(−y⁺/26)]²/Sc_t；R=∫₀^δdy/D(y)，以72段對數積分。u_τ 使用保存 U、nut 的壁面剪應力近似 sqrt[(ν+ν_t,wall)|U_t|/δ]，不是宣稱高精度量測壁面剪力。

法向朝固體為正。漂移含重力（天花為負，不能直接忽略）、Talbot 熱泳、低 Stokes 數平均流慣性滑移 −τ_p DU/Dt。可選顯式湍泳近似預設關閉，避免與湍流慣性經驗閉合重複。

湍流慣性輸運採 Wood 型 v_imp=u_τ min[0.13,4.5×10⁻⁴(τ_pu_τ²/ν)²]。這個公式來自管道沉積，移植至廚房天花未經校準；不能因寫入方程就稱「真實慣性碰撞已驗證」。阻力對照閉合以 R_eff=1/(1/R+v_imp)，v_dep=α v_n/[1−exp(−v_n R_eff)]；v_n→0 時回到 α/R_eff，逆重力時沉積受抑制。主篩查採把 v_imp 與漂移阻力通量相加的 parallel 閉合：慣性接觸不應被當作布朗擴散而必須通過同一薄層阻力。這仍是管道經驗式的天花外推，對慣性接觸中的逆重力影響並非精確軌跡計算；兩種閉合均保留作不確定性比較，不是嚴格數學上下界。

天花溫度使用既有熱邊界值；±5 K 敏感度只變粒子熱泳驅動，不重解固體或氣流熱場。這不是石膏板和鋁板的熱傳比較。Sc_t、附著 α、熱泳、慣性、近壁積分及時序插值均有實際計算檢查。

## A / B 方法比較

採 Eulerian：五粒徑可直接延伸，沒有極小沉積事件的有限粒子抽樣噪音，成本為每步五個標量方程。可重用原 U/T/k/nut/phi，但 10 s 保存間距未包含未保存的瞬態細節。Lagrangian 需要近壁非均勻隨機湍流、布朗力/熱泳、慣性和附著模型；trap 只處理已接觸壁面的粒子，不能代替細粒子走到壁面的機理。沒有做 Lagrangian 粒子數收斂，也不把它列作本輪驗證通過項目。

## 油滴與冷凝分開

本輪只追蹤已形成、定密度球形油滴。氣態/半揮發性有機物的表面凝結、吸附、氣粒分配、化學轉化沒有建模。不能以牆溫較低便把氣態濃度全算成沉積；需組分排放、飽和蒸氣壓/活度、表面傳質及表面溫度。亦未計靜電、漆面微觀粗糙度、反彈/再懸浮、油膜改變黏附和長期老化。

## 一手來源與適用限制（2026-10-21查閱）

- [Lai & Nazaroff 2000](https://doi.org/10.1016/S0021-8502(99)00536-4)：室內光滑面布朗/湍流擴散與沉降。本次是阻力近似，不是逐式重現。
- [Sippola & Nazaroff, LBNL report, Particle Deposition from Turbulent Flow: Review of Published Research and Its Applicability to Ventilation Ducts in Commercial Buildings](https://www.osti.gov/servlets/purl/815478)，4.3.2–4.3.4、式53及71–74：Wood係數、水平面方向與管道外推限制。本次 hybrid closure 另作假設，不宣稱等同報告模型。
- [Guha 2008](https://doi.org/10.1146/annurev.fluid.40.111406.102220)：統一 Eulerian / Lagrangian 沉積機理及適用邊界。
- [Talbot et al. 1980](https://doi.org/10.1017/S0022112080001905)：熱泳係數；沿用既有油滴導熱率0.17 W/(m K)。
- [Applied Sciences 2022,12,8896](https://www.mdpi.com/2076-3417/12/17/8896)：原專案採0.12 mg/s/頭的商用油滴來源。兩頭14.4 mg/min只作參考，另列0.25–4倍；不是這個住家的量測或通用中式炒菜排放率。
- [Chen et al. 2023](https://pubmed.ncbi.nlm.nih.gov/37847366/)：六道中式菜的粒徑/捕集實驗；沿用質量中位徑2.5 µm、假設GSD1.5，另以GSD1.2/2.0重組。沒有把PM2.5量測當作全部純油脂。
- [Gyproc moisture-resistant board 資料](https://www.gyproc.ie/documents/product-datasheets/pds-gyproc-ireland-moisture-resistant.pdf)：防潮板的使用條件不等於防油或可無限濕擦。
- [Benjamin Moore Kitchen & Bath Satin W322 TDS](https://media.benjaminmoore.com/WebServices/prod/assets/production/datasheets/TDS_0322/W322_TDS_US.pdf)：特定廚衛漆產品用途；不提供油滴黏附率，不能把廚衛適用轉成 CFD α 數值。
- [Dulux Kitchen Matt](https://www.dulux.co.uk/en/products/dulux-easycare-kitchen-matt/cccId_1935922)：有耐油污的啞光配方，顯示不能只按光澤判斷性能；普通啞光需查實際產品的耐擦洗規格。
- [Hunter Douglas 金屬吊頂清潔維護](https://assets.hunterdouglasarchitectural.eu/documents/ceiling/Installation%20Manuals/_default/Cleaning-and-Maintenance-Luxalon-multipanel.pdf)：應依實際塗層和清潔劑相容性維護；這不是所有烤漆鋁板的沉積率量測。

三種表面都以相同 α=0.1/0.5/1 掃描；無資料支持任何材料的 α 排名。清潔模型同樣不預設某材料有較高效率。材料比較的可靠部分是產品獨立規格與可維護性，並非 CFD 預測漆面黃變、老化或可擦洗性。

## 停止煮食後的清空積分

240 s 仍有懸浮油滴；以固定末端熱流場解線性輸運的殘留清空 resolvent，逐面累積後續 sink，剩餘質量仍列帳。代數尺度 H=10¹² s 不是物理觀測時長，另以10⁹檢查正則敏感度。這可避免把所有剩餘量都配給天花或完全忽略它，但固定殘餘熱場仍是近似；詳見 README。
