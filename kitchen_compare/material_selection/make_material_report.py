from pathlib import Path
import json,base64
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'deliverables'

def main():
 s=json.loads((OUT/'material_summary.json').read_text());a,b=s['baseline']
 ksteel=s['materials']['SUS304']['k'];rwood=.022/.14;rsteel=.0015/ksteel
 mesh_error=max(abs(group[0]['exposed_surface_t95_min']-group[-1]['exposed_surface_t95_min'])/group[-1]['exposed_surface_t95_min'] for group in [[r for r in s['mesh'] if r['covered']==c] for c in [False,True]])
 report=f'''# 22 mm顆粒板與不鏽鋼覆顆粒板：升溫條件比較

使用者確認的不鏽鋼方案是**1.5 mm鋼板緊贴或薄膠貼在22 mm顆粒板上，沒有刻意留空氣層**。本次比較裸顆粒板與覆鋼顆粒板；原廚房工況仍為雙爐每頭對空氣放熱2000 W、合計4000 W，靠左、兩門關、每門底縫5 mm。

**本次是新增的板材一維暫態導熱有限元，不是廚房從開火起的完整暫態CFD。迎熱空氣假設在t=0瞬時達到一個固定溫度，板材初始25°C；沒有把原127.6°C熱點當作已知的真實受熱條件。**

## 物性假設與來源

尚未指定產品，因此先用公開資料中的代表值，不是所選板材的實測規格。

| 材料 | 厚度 | 導熱係數 W/(m·K) | 密度 kg/m³ | 比熱 J/(kg·K) |
|---|---:|---:|---:|---:|
| 顆粒板 | 22 mm | 0.14 | 600 | 1700 |
| SUS304不鏽鋼 | 1.5 mm | 16.3 | 7930 | 500 |

顆粒板參考[UNILIN歐洲技術認證ETA 10/0019，表1引用EN12524](https://cms.esi.info/Media/documents/155730_1334915100717.pdf)。不鏽鋼參考[日本製鐵SUS304物性表](https://www.nipponsteel.com/en/product/stainless/steel_shape/)。牌號、板密度、含水率、飾面和溫度會影響實際物性。

## 加覆鋼板改變甚麼

| 指標（每平方米） | 裸22 mm顆粒板 | 1.5 mm鋼＋22 mm顆粒板 |
|---|---:|---:|
| 質量 kg/m² | 13.20 | 25.095 |
| 熱容量 kJ/(m²·K) | {a['heat_capacity_per_area_J_m2K']/1000:.2f} | {b['heat_capacity_per_area_J_m2K']/1000:.2f} |
| 材料穿板導熱熱阻 m²·K/W | {rwood:.6f} | {rwood+rsteel:.6f} |

覆鋼增加約{100*(b['heat_capacity_per_area_J_m2K']/a['heat_capacity_per_area_J_m2K']-1):.1f}%熱容量，但只增加約{100*rsteel/rwood:.3f}%材料導熱熱阻。在接觸良好、相同受熱與散熱條件下，鋼層本身幾乎不增加隔熱；升溫可因額外熱容量變慢。

## 從25°C起的条件升溫時間

示例條件：迎熱面對流係數8 W/(m²·K)，背面對25°C櫃內空氣散熱係數3 W/(m²·K)，兩層完全熱接觸。這些係數是待確認的假設，不是從原廚房流場校準得到的值。櫃內空氣維持25°C亦是固定散熱參考，沒有計入整個櫃內空氣升溫。

以「達到各面最終溫升的95%」定義接近穩態。例如最終溫升20°C，95%時仍差1°C；不是溫度等於95%的°C讀數。

| 位置 | 裸顆粒板 | 覆鋼顆粒板 |
|---|---:|---:|
| 迎熱面達95%最終溫升 | {a['exposed_surface_t95_min']:.0f}分鐘 | {b['exposed_surface_t95_min']:.0f}分鐘 |
| 顆粒板靠熱源一面達95% | {a['wood_hot_face_t95_min']:.0f}分鐘 | {b['wood_hot_face_t95_min']:.0f}分鐘 |
| 顆粒板背面達95% | {a['wood_back_face_t95_min']:.0f}分鐘 | {b['wood_back_face_t95_min']:.0f}分鐘 |

時間描述板材对已穩定熱空氣的響應。實際開火時空氣與板材同時升溫，因此不能直接將此表當作廚房開火後時間。線性、定物性且不含輻射的模型中，達相同最終溫升百分比的時間與所施加的熱空氣溫升幅度無關；絕對°C則取決於該空氣溫度。

![條件升溫與厚度溫度分布](material_transient_comparison.png)

## 溫度應如何閱讀

本模型使用歸一化溫升。若固定熱空氣溫度為Ta：

`T(t)=25+模型歸一化溫升(t)×(Ta−25)`。

示例散熱條件下，裸板迎熱面最終溫升為熱空氣溫升的{a['exposed_surface_steady_rise_fraction']:.4f}倍；覆鋼方案為{b['exposed_surface_steady_rise_fraction']:.4f}倍，兩者非常接近。覆鋼方案顆粒板靠鋼的一面為{b['wood_hot_face_steady_rise_fraction']:.4f}倍，與鋼板外面溫差很小。

例如**僅假設**熱空氣固定50°C，裸板迎熱面約{25+25*a['exposed_surface_steady_rise_fraction']:.2f}°C，覆鋼迎熱面約{25+25*b['exposed_surface_steady_rise_fraction']:.2f}°C。這是公式示例，不是廚房預測，更不表示原平均49.7°C適用於任何板材。

## 假設與數值敏感性

迎熱係數5、8、15 W/(m²·K)，背面係數0、1、3、5 W/(m²·K)的示例組合中，迎熱面95%時間範圍：裸板約{s['film_sensitivity_t95_ranges']['False'][0]:.0f}–{s['film_sensitivity_t95_ranges']['False'][1]:.0f}分鐘，覆鋼板約{s['film_sensitivity_t95_ranges']['True'][0]:.0f}–{s['film_sensitivity_t95_ranges']['True'][1]:.0f}分鐘。這是參數敏感性範圍，不是統計信賴區間，也不是所選產品的保證範圍。

接觸熱阻另測0、0.001、0.01、0.1 m²·K/W。後兩項僅為敏感性參數，不代表已知薄膠或指定空氣層。增加接觸熱阻會降低熱傳入顆粒板的程度，並改变鋼面與木面溫差；有空氣層時還要計入層內輻射，不能單憑此表選氣隙厚度。

採一維線性有限元、每單位面積的相容熱容量矩陣與導熱矩陣，兩側Robin對流邊界。時間演化由半離散矩陣的特徵模態直接求解；95%時間用求根計算，而非30秒匯出取樣插值。穩態解與解析串聯熱阻解逐項核對。顆粒板80、160、320單元的迎熱面95%時間差異小於{100*mesh_error:.4f}%，支持此一維假設下的數值穩定性。

## 選材解讀與剩餘限制

緊貼鋼板的穿板導熱熱阻很小；在本模型相同對流條件下，額外鋼層主要增加儲熱量，穩態顆粒板迎熱面並未明顯變冷。鋼板的橫向熱擴散、表面對輻射的吸收／反射可能改變真實熱點，但本一維比較尚未計算。表面飾面、膠層與背面散熱都需要在產品確定後更新。

實際側面的°C、W與開火到接近穩態的時間，仍需通過局部守恆檢查的流場、暫態熱浮力與固體導熱耦合，以及火焰／鍋具輻射條件。本次結果用於比較板材自身的熱響應，不取代原廚房模型，也沒有驗證原熱點位置。

## 輸出

`material_baseline.csv`提供示例條件主要結果；兩個`*_history.csv`為時間與各面歸一化溫升；`film_sensitivity.csv`、`contact_sensitivity.csv`、`mesh_validation.csv`提供邊界與數值驗證；`solid_transient.py`及`material_summary.json`保留程式與全部輸入。完整資料包不依賴原廚房CFD檔案。
'''
 (OUT/'material_comparison_zh.md').write_text(report)
 from markdown_it import MarkdownIt
 body=MarkdownIt('commonmark').enable('table').render(report)
 embedded='data:image/png;base64,'+base64.b64encode((OUT/'material_transient_comparison.png').read_bytes()).decode()
 body=body.replace('src="material_transient_comparison.png"',f'src="{embedded}"')
 css='body{font-family:"Noto Sans CJK TC",sans-serif;font-size:11px;line-height:1.6;color:#172a3d}h1{font-size:22px}h2{font-size:16px;break-after:avoid}table{border-collapse:collapse;width:100%;font-size:10px}th,td{border:1px solid #c3cdd8;padding:6px}th{background:#edf2f7}img{max-width:100%}tr{break-inside:avoid}code{white-space:pre-wrap}@page{size:A4;margin:16mm}'
 document=f'<!doctype html><html lang="zh-Hant"><meta charset="utf-8"><title>顆粒板覆鋼條件比較</title><style>{css}</style><body>{body}</body></html>'
 (OUT/'material_comparison_zh.html').write_text(document)
 from playwright.sync_api import sync_playwright
 with sync_playwright() as pw:
  browser=pw.chromium.launch(executable_path='/usr/bin/chromium',args=['--no-sandbox']);page=browser.new_page();page.set_content(document,wait_until='load');page.evaluate('document.fonts.ready');page.pdf(path=str(OUT/'material_comparison_zh.pdf'),format='A4',print_background=True);browser.close()
 print(dict(front_t95_minutes=[a['exposed_surface_t95_min'],b['exposed_surface_t95_min']],maximum_mesh_relative_difference=mesh_error))

if __name__=='__main__':main()
