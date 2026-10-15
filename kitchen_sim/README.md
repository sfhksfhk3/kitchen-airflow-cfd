# 廚房氣流有限元模型

先開啟 `interactive_3d.html`：可離線旋轉模型，切換四個門狀態、速度／溫度／相對壓力及截面。3D 圖例的「定量切片」可點擊顯示或隱藏。

閱讀 `report.pdf` 或 `report.html` 查看假設、結果與數值檢查。`four_cases_velocity.png`、`four_cases_temperature.png` 是共用色階的四工況比較。

資料位於 `results/`：每工況提供逐點 CSV、四面體 VTU、NumPy NPZ、收斂 CSV 和指標 JSON。ParaView 可直接讀取 VTU，使用 Slice、Glyph 或 Stream Tracer 查看。`probes.csv` 匯出六個固定位置；`geometry.obj` 包含櫥櫃和煙機的簡化實體，單位均為公尺。

四個基準工況：
- `both_open`：兩門全開。
- `both_closed`：兩門全關，各門底縫假設 5 mm。
- `left_open`：左開右關。
- `right_open`：右開左關。

這是條件式初步計算。每燃氣爐向空氣釋放 600 W 顯熱為假設；風機曲線與混合強度未校準。關門補風以真正 5 mm 縫計算阻力，但室內入口鋪展至 0.10 m 帶，未解析門縫噴流。局部峰值仍受網格影響；報告提供較細網格比較，不宣稱網格無關或已驗證實機捕煙率。

## 重現

Python 3.12 環境：
```bash
python -m pip install -r requirements.txt
python solve.py --case all --h .16
python solve.py --case both_open --h .12 --suffix _refined --warm results/both_open.npz
python validate_exports.py
python visualize.py
python make_report.py
```

中文圖表需可用的 CJK 字型。產生 HTML 報告需 Pandoc；PDF 另需 Chromium 和 Playwright。已匯出的 HTML、PNG、PDF 不需要安裝 Python。

求解器每十次迭代保存 `results/<case>_checkpoint.npz`。如細網格運算中斷而檔案仍在，可沿用相同設定續算：
```bash
python solve.py --case both_open --h .12 --suffix _refined --warm results/both_open_refined_checkpoint.npz
```
修改模型參數後不可接用舊檢查點；程式會核對設定。重算會覆寫同名結果。
