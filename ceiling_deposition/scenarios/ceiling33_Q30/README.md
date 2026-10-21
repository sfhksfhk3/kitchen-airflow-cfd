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

30m³/min是指定實際流量；求解器不在此情境套用1300Pa截斷，會計算所需入口靜壓。全關門5mm縫的初步20秒測試需要約1.8kPa，不能據此聲稱標稱最大靜壓1300Pa的風機能達到30m³/min；風管／濾網等外部阻力尚未加入。

原網格嚴格檢查仍有失敗項目，原壁面模型亦未經烹調油脂天花沉積量測校準。完整質量守恆不能代替真實沉積量驗證；不要用小數值判定石膏板肯定無積油問題。
