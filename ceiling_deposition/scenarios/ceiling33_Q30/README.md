# 天花33°C／實際排風30m³/min更新

此情境保留原始專案和舊沉積結果，另計算熱流場及油滴沉積。不是把舊溫度結果改成33°C。

- 用戶條件：天花表面33°C；出口體積流量0.5m³/s；兩爐各放熱2000W；保留兩種爐位和四種門狀態。
- 排放及材料：`parameters.json`、`emission_scenario.json`、`materials.json`、`SOURCES_zh.md`。公開資料沒有可直接套用的三種材料油滴附著率，故保留共同敏感度，不捏造材料常數。
- 熱流場：`../../flow_sensitivity/*_T33_Q30`；原網格單元不變，只重排外部面以建立ceiling獨立壁面。新面ID不得直接索引舊版面數據。
- 油滴：`../../runs/*_T33_Q30*`；耦合空氣求解中的舊Yp欄位只作相容性占位，`ceilingReplayFoam`才是本情境正式油滴求解。
- 新結果輸出：`results`的逐面CSV、PVD/VTP、質量收支和長期敏感度；`analysis_zh.pdf`、`interactive_3d.html`、`figures`。先檢查`PROGRESS.json`，未顯示COMPLETE時不得把尚未完成的批次當成最終結果。

## 執行

需要原專案恢復到`/workspace/kitchen_openfoam`，以及已設定OpenCFD OpenFOAM2406環境的`kitchen-openfoam-2406`容器，掛載同一路徑`/workspace`。不重新生成網格。宿主Python依賴沿用原專案NumPy/SciPy/Pandas/Matplotlib/Plotly/Playwright；報告另需markdown-it-py，Chromium在`/usr/bin/chromium`。兩核心環境採兩個批次並行。

```bash
cd /workspace/kitchen-airflow-cfd
python ceiling_deposition/scripts/prepare_33C_Q30.py
# 在容器內、OpenFOAM環境已載入時：
docker exec kitchen-openfoam-2406 bash -lc 'cd /workspace/kitchen-airflow-cfd/ceiling_deposition/src/kitchenCeilingAirFoam && wmake'
# ceilingReplayFoam、ceilingTailFoam亦需先依原擴展說明編譯。
python -u ceiling_deposition/scripts/update_33C_Q30.py
```

正式批次會沿用已保存的熱流場時間目錄；完整油滴結果有RUN_COMPLETE標記時不重算。油滴尚未完成的執行不能僅靠Mdep檔案任意重啟，需依run_suite的完成判定從該單位響應初始時間再解，避免錯接累積質量。別覆蓋既有案例名稱的參數後又沿用完成標記；改參數應建新名稱。

細網格低成本啟動測試使用原專案已有細網格；主批次不會重新生成或細化網格。天花33°C不涉及固體導熱求解。長期模型不是塗料發黃／厚度／壽命預測；更不是對材料的實測認證。

## 風量與可信度

30m³/min是指定實際流量；求解器不在此情境套用1300Pa截斷，會計算所需入口靜壓。全關門5mm縫的完整240秒案例末段需要約1.88kPa，不能據此聲稱標稱最大靜壓1300Pa的風機能達到30m³/min；風管／濾網等外部阻力尚未加入。

原網格嚴格檢查仍有失敗項目，原壁面模型亦未經烹調油脂天花沉積量測校準。完整質量守恆不能代替真實沉積量驗證；不要用小數值判定石膏板肯定無積油問題。


## 2026-10-22 完成狀態

八個主工況、兩個全關門佈局各α=0.1/0.5的四組附著率試算、完整歷程油滴時間步／存檔間距測試均已完成。另有匹配初始化的10秒粗細網格診斷、20秒混合啟動診斷，以及共用80秒初始熱流場的24/30m³/min短窗測試。後者不是整餐增風改善率。

- `analysis_zh.pdf`：中文分析及六項問題的回答，含原始來源連結。
- `interactive_3d.html`：完全離線；12案例×3粒徑分布，顯示總量、平均、面積加權P95和最大值；可調單次排放量、天數、頻次及清潔情境。
- `figures/all_figures_1to1.pdf`：八張正方形圖；四種門狀態各將靠左／置中合成一張同色階圖。
- `results/scenario_summary.csv`：八工況總量與去向；`long_term_per_1mg_emitted.csv`優先提供單位排放響應，另有4320列材料／清潔／頻次／排放敏感度。
- `results/<case>/ceiling_deposition.pvd`：實際0–240秒逐面J及M；`ceiling_complete_event.vtp`是加上固定末端流場清空積分的靜態結果，沒有虛構未來物理時間。
- `EXPORT_VALIDATION.json`：12案例；靠左每案533面、置中每案550面，4.245m²暴露天花；CSV、VTK與積分收支核對。
- `PRESENTATION_VALIDATION.json`：38組UI對CSV核對、零排放測試、清潔控制及視角按鈕；離線無HTTP請求／JavaScript錯誤。

**功能可運行，不代表物理定量已通過。** 粒子相對收支殘差最大約3.8×10⁻¹¹；但匹配初始化的10秒粗細網格沉積差7.61倍，油滴步長0.1→0.05秒差6.90%，流場存檔2→4秒差48.51%。原生天花壁函數y⁺最高約1277，沒有解析沉積薄層；近壁閉合未用廚房油脂沉積實驗校準。完整熱流場時間步和網格收斂尚未證明。不可用極小模型值宣稱石膏板不會積油，更不能預測黃變時程。

單獨重新輸出／驗證（不重跑已完成CFD）：

```bash
python ceiling_deposition/scripts/report_33C_Q30.py
python ceiling_deposition/scripts/validate_exports_33C_Q30.py
python ceiling_deposition/scripts/validate_presentation_33C_Q30.py
```

`ceiling33_Q30_review.zip`保留專案相對目錄，包含本情境全部CSV／ParaView、PDF、離線3D和完整擴展程式碼。`ceiling33_Q30_native_cases.tar.gz`另存原生流場、網格、油滴解與清空案例，保留符號連結；解至原專案根目錄，再以較新的閱讀包覆蓋程式與報告。容器中仍使用相同`/workspace`路徑；移動後需依原擴展的relink_frames.py說明修正連結。閱讀CSV／VTK／HTML本身不依賴原生大檔。
