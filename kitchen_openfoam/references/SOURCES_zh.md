# 參數來源與採用方式

使用者要求自行查閱文獻，不再追問實驗條件。下列是文獻資料及本次選定的工程近似；文獻並不構成此廚房的實測驗證。

| 參數／模型 | 來源 | 本次採用方式 |
|---|---|---|
| 油滴排放0.12 mg/s | [Simulation Study of the Capture and Purification Performance of Exhaust Fume Systems in Chinese Commercial Kitchens, 2022](https://www.mdpi.com/2076-3417/12/17/8896)，第2.2節。已取得官方公開PDF，確認為直徑2.5µm球粒子的估計排放率。文中表1的3.3×10⁻⁷ kg/s是其他參考案例，沒有混用。 | 取7.2 mg/min/頭作總油滴質量參考情境；兩頭14.4 mg/min。由商用案例轉作住宅計算，是工程假設，另列0.25–4倍排放敏感性。 |
| 六道中式菜平均體積中位徑2.5µm、眾數3.3µm | [Chen et al., Emission and capture characteristics of Chinese cooking-related fine particles, 2023](https://doi.org/10.1007/s11356-023-30380-4)；[PubMed](https://pubmed.ncbi.nlm.nih.gov/37847366/) | 以同密度球滴的質量中位徑2.5µm作參考；GSD1.5是本次假設，另比較1.2、2.0。沒有把眾數和中位徑當成同一量。 |
| 油加熱階段質量主要在1–4µm | [Gao et al., Volume-based size distribution of accumulation and coarse particles, 2013](https://www.sciencedirect.com/science/article/pii/S0360132312002740) | 支持以微米級油滴作質量比較；不推論所有中式菜式均有相同分布。 |
| 住宅中式烹調PM2.5排放2.056±3.034mg/min | [Chen et al., 2018](https://pubmed.ncbi.nlm.nih.gov/29302961/) | 作另一個排放量級參考；研究涵蓋不同烹調法。沒有把PM2.5排放直接當成完整油滴總質量排放，也沒有稱為純炒菜定值。 |
| 烹調顆粒密度約850kg/m³、數量分布可有約50nm模態 | [Cai et al., 2020, ACP 20,12721](https://acp.copernicus.org/articles/20/12721/2020/acp-20-12721-2020.html)，3.4節 | 密度850kg/m³為代表性近似。五種計算粒徑0.1–10µm未涵蓋全部<0.1µm粒子，不能預測總粒子數／超細粒子暴露。 |
| 油滴導熱係數約0.17W/(m·K) | [Temperature dependence of thermal conductivity of vegetable oils, 2017](https://www.sciencedirect.com/science/article/pii/S1290072916303970) | 代表性常數，僅用於Talbot熱泳係數；不作油滴蒸發／加熱方程。 |
| 熱泳 | [Talbot et al., 1980](https://doi.org/10.1017/S0022112080001905)；[作者機構公開稿](https://escholarship.org/uc/item/22f5r6cz) | 採含滑移修正的Talbot係數，Cs=1.17,Cm=1.14,Ct=2.18，依粒徑、氣體及油滴導熱比計算。 |
| 沉積 | [Lai & Nazaroff, 2000](https://doi.org/10.1016/S0021-8502(99)00536-4) | 參考布朗／湍流擴散及重力沉降的壁面通量思想。本次實作是明示的van Driest阻力積分近似，並非逐式重現該論文。黏附率基準1，沒有反彈、再懸浮或慣性撞擊。 |
| k–ω SST、Prt | [OpenFOAM SST模型](https://doc.openfoam.com/2306/tools/processing/models/turbulence/ras/linear-evm/rtm/kOmegaSST/)；[alphatWallFunction](https://doc.openfoam.com/2312/tools/processing/boundary-conditions/rtm/derived/thermal/alphatWallFunction/) | 使用安裝的OpenCFD2406程式碼，壁面Prt0.85；Sc_t0.7為假設，需敏感性評估。 |

雙爐對空氣放熱各2000W、幾何尺寸、擋板和5mm門縫均優先採使用者已確認設定。沒有用文獻爐具總燃氣額定功率覆蓋這些值。
