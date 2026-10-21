# 天花油脂沉積擴展

沿用 v1.0.0 幾何、原網格、煙機及既有熱流場。新增 `ceilingReplayFoam`（保守 Eulerian 油滴重播）和 `ceilingTailFoam`（固定末端流場的殘留清空積分）；原求解器與舊成果不修改。

**驗證未收斂：後段粗／細分支天花沉積約相差27倍；不能用此版定量判定材料長期積油。數據與可執行工具均保留供檢查及校準。**

先看 `reports/01_project_audit_zh.md`，再看 `reports/ceiling_deposition_analysis_zh.pdf`。主篩查情境為名稱含 `_parallel` 的八工況，與阻力閉合的八個對照都保留。這是未實測校準的條件研究，不能用極小數值宣稱材料不會黃變或積油。

## 重現

既有資料已從原 Release 的 `kitchen_openfoam_full_snapshot.tar.gz`（SHA256 見根目錄 docs/SHA256SUMS）恢復至 `/workspace/kitchen_openfoam`。原場只讀使用，不用重跑。其他位置可設 `KITCHEN_FLOW_ROOT`，或直接把案例路徑傳給 `scripts/prepare_replay.py`。新工作位於 `/workspace/kitchen-airflow-cfd/ceiling_deposition`；容器以相同路徑掛載 `/workspace`。

需要 OpenCFD OpenFOAM 2406（原映像 `opencfd/openfoam-default:2406`）、Python 套件見 requirements.txt。執行前確認容器 `kitchen-openfoam-2406` 可使用此工作區。使用既有容器：

```bash
docker start kitchen-openfoam-2406
docker exec kitchen-openfoam-2406 bash -lc 'export FOAM_USER_APPBIN=/workspace/kitchen-airflow-cfd/ceiling_deposition/bin; cd /workspace/kitchen-airflow-cfd/ceiling_deposition/src/ceilingReplayFoam; wmake; cd ../ceilingTailFoam; wmake'
```

先做低成本試算（不覆蓋同名已存在案例）：

```bash
python scripts/prepare_replay.py /workspace/kitchen_openfoam/cases/corner_both_closed demo --end 20 --dt 0.2 --options '{"impactionMode":"parallel"}'
docker exec kitchen-openfoam-2406 bash -lc '/workspace/kitchen-airflow-cfd/ceiling_deposition/bin/ceilingReplayFoam -case /workspace/kitchen-airflow-cfd/ceiling_deposition/runs/demo'
python scripts/validate_analytic.py
```

正式套件：`python scripts/run_suite.py`。已完成案例以 `RUN_COMPLETE.json` 和最終質量收支確認後沿用；不重跑成功案例。可中斷排程後重啟，執行中的求解器會先等候；失敗案例才自初始場重算。輸運重播不支援任意中途時間原地續算累積場，請勿把 startFrom 改成 latestTime 當作累積沉積重啟。不要同時啟動兩個 suite 排程。

`run_coupled_interventions.py` 準備及執行五個短時新氣流（內含 `prepare_flow_sensitivity.py` / `run_matched_flow_controls.py`）；`run_interventions_parallel.py` / `replay_matched_controls.py` 處理匹配的油滴輸運。提高風量比較只適用約100–140 s、相同初始檢查點；不是整次煮食所有時間均提高20%的完整測試。`verification/` 留有本次熱源關閉分支紀錄。

```bash
python scripts/clearance_integral.py
python scripts/process_ceiling.py
python scripts/build_analysis.py
python scripts/long_term.py
python scripts/visualize_ceiling.py
python scripts/make_interactive.py
python scripts/make_ceiling_report.py
```

## 讀取資料

- `results/<案例>/ceiling_deposition.pvd`：ParaView 開啟原網格面的瞬態熱圖；選 Cell Data。M 及 J 是面資料，不是把體積濃度上色。
- `results/<案例>/ceiling_complete_event.vtp`：240 s 內沉積＋固定240 s熱流場的殘留油滴清空積分。獨立靜態檔案，**沒有虛構未來物理時間**。
- `ceiling_faces.csv`：240 s 最後面密度、通量及時段平均通量；原始逐時 J（kg/m²/s）、M（kg/m²）見 `runs/<案例>/diagnostics/ceilingFaces_*.csv`。
- `ceiling_complete_event_faces.csv`：事件累積、觀測時段與長尾分項；原始面ID、座標、面積和五粒徑；M_unit 是 mg/m²／每排放1mg，M_reference 是文獻參考情境的 mg/m²。
- `primary_complete_event_summary.csv`：八工況最大面、面積加權平均/P95、總質量、質量佔比和去向。
- `numerical_and_physical_sensitivity.csv`：真正跑過的時間步、近壁閉合參數、時序間距、熱泳/慣性和後段加密對照。後段加密不等於完整網格獨立。
- `long_term_per_1mg_emitted.csv`：優先使用的長期單位排放響應；絕對情境另見 `long_term_scenarios.csv`。
- `ceiling_deposition_interactive_3d.html`：完整離線，不需外部網路；可切案例、假設附著率/粒徑分布、排放劑量、天數與清潔情境。

天花原本混在 walls；以原面標籤、法向與 z=2.2 m建立 `ceilingFaces` faceSet 和 `ceiling` faceZone。原面序、拓撲、座標與物理邊界均保留。實際暴露天花面積由原網格求得（高櫃與上櫃佔用區沒有空氣接觸面），不將它們錯當作零積油實測。

## 清空積分的時間含義

固定末端流場時，線性輸運 dM/dt=−A M 的累積停留量 F 滿足 A F=M₀。用很大的代數正則尺度 H 解 (A+I/H)F=M₀，並保留未清空質量在收支內。H=10¹² s 是**代數參數，不是模擬煮食了三萬年**；不能拿其原生目錄名作物理動畫時間。只將積分出的累積沉積合併至 complete-event 靜態熱圖。另以 H=10⁹ 檢查正則誤差。假設煙機繼續運作且末端流場可代表後續；不是再解固體冷卻或漆面老化。

## 材料和長期模型

`config/long_term.json` 可改每天次數、30/180/365日、每次排放mg、清潔間隔和效率。附著率0.1/0.5/1在兩個全關門佈局實際重解；開門只基準1。三種材料共用這個範圍；沒有實測支持黏附率排名。清潔效率亦是情境，與材料沒有自動對應。0.25–4倍排放是原商用文獻參考值的敏感度，不是住家排放的保證範圍。氣態冷凝、粗糙度、再懸浮和塗層老化未建模。

## 原生案例封裝

原生案例壓縮包另含新 runs/、flow_sensitivity/、tail_sources/ 和診斷紀錄；唯讀frame以符號連結連回舊場，避免重複打包。請先恢復原Release快照，解包至同一專案根目錄；若位置改變，執行 `python scripts/relink_frames.py --flow-root /你的路徑/kitchen_openfoam`。ParaView VTP/PVD與閱讀包本身不依賴這些連結。尾部代數積分案例有極大的數值目錄名，請只使用其 complete-event 靜態熱圖，勿當作物理時間動畫。
