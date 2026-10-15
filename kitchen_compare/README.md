# 廚房油煙滯留：角落與中央

先開啟 `interactive_comparison.html`，可離線旋轉3D、切換角落／中央、四種門狀態、濃度／風速／溫度／壓力及時間。油煙是相同源強的被動示蹤物；沒有實際排放率，濃度是c/R（s/m³），不是已知mg/m³。

`report.pdf` 或 `report.html` 是完整中文分析，含優劣、幾何、指標、數值核對與不確定性。`comparison_summary.csv` 是八個主工況的摘要。主結果用h=0.14m，粗網格h=0.18m作敏感性比較。

每個燃氣爐對空氣放熱2000 W，兩爐共4000 W。煙機90×30cm、距灶面70cm；擋板長10cm向下45°，只在煙機90cm寬內。角落整體中心x=1.15m；中央x=2.0m。兩者均依新三視圖重算，未直接沿用舊版1200W的結果。

主要判斷：中央總滯留的微小優勢未穩定，走道濃度則普遍較高。若優先降低向走道的油煙散佈，目前模型較支持保留靠左位置、在無外部穿堂風時兩門補風。精確尾端清除秒數對網格與混合假設敏感，不能作為關火後可關煙機的延時建議。

## 資料

- `corner/`、`center/`：各自的config、OBJ、求解程式及results。
- `*_refined_combined.vtu`：主結果供ParaView讀取，含FEM速度、空氣溫度、相對壓力、示蹤濃度，以及cell_data的RT0守恆輸運速度。
- `*_refined.npz`／`*_refined_nodes.csv`：原始FEM網格與節點場。
- `*_refined_tracer.npz`：示蹤濃度快照、時間、穩態濃度、RT0速度和質量權重。
- `*_refined_tracer_nodes.csv`／`*_refined_tracer_history.csv`：逐點c/R與每0.25s的存量、排出、走道濃度時序。
- `analysis_summary.json`、`validation.json`：完整比較及獨立核對。
- `grid_comparison.csv`、`diffusion_sensitivity.csv`、`time_step_comparison.csv`：網格、Sc和時間步長檢查。

ZIP內含八個主工況完整場資料、粗網格的驗證資料，以及Sc/時間步长敏感性的指標與時序；省略重複的敏感性場快照，可用附帶程式重現。舊專案保留於工作環境 `/workspace/kitchen_sim`，本次另存於 `/workspace/kitchen_compare`。

## 重現

見報告第12節。Python3.12，先`python -m pip install -r requirements.txt`。Gmsh需系統GLU庫，中文圖片需CJK字型；產生報告需Pandoc，PDF另需Chromium和Playwright。已匯出HTML、PNG、PDF可直接離線使用。

求解器每十次迭代保存檢查點；相同配置下可用`--warm <checkpoint.npz>`續算。修改幾何、熱負荷或邊界後需重新求解，不能只改標籤。

這是條件式數值比較：固定基準密度、未校準有效混合、無室外風壓、門底縫鋪展近似、無油滴沉降／壁面沉積。數值收斂與守恆不是實機驗證。
