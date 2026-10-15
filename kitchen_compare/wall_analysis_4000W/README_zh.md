# 4000 W高櫃侧面分析

最新設定：每爐對空氣放熱2000 W，雙爐4000 W；靠左布局，兩門關，每門底縫5 mm。此資料包不含中途9000 W結果。

先閱讀deliverables/wall_heat_report_zh.pdf，或直接開啟deliverables/wall_interactive_3d.html。

壓縮包保留kitchen_compare目錄結構。解壓後，在kitchen_compare中執行：

python wall_analysis_4000W/analyze_wall.py
python wall_analysis_4000W/make_wall_report.py

分析使用corner/results已完成的both_closed與both_closed_refined原始解，不會重新執行CFD。corner/solve.py與build_mesh.py是原有限元求解器和幾何網格程式，config.json包含全部設定。兩級網格均提供原始解、網格和收斂資料。

需要numpy、scipy、scikit-fem、meshio、matplotlib、plotly；產生報告另需markdown-it-py、playwright、pypdf、Chromium和Noto CJK字體。scikit-fem版本12，gmsh版本4.15（僅重新建網格需要）。

模型未求解固體導熱或直接輻射，不能把模型邊界溫度當作實際板溫。詳見中文報告。

## 熱點定位補充檢查

執行python wall_analysis_4000W/hotspot_diagnostic.py，可重現使用已存局部守恆RT0速度的固定流場熱敏感性檢查。輸入已另存於deliverables/hotspot_RT0_input.npz，不需重新做油煙示蹤物計算。此檢查沒有重新耦合浮力，不是替代原模型的新完整解；熱點位置與功率敏感性見更新報告。
