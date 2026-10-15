# 板材条件比較

裸22 mm顆粒板 versus 1.5 mm SUS304緊貼22 mm顆粒板；初始25°C。

python solid_transient.py
python make_material_report.py

這是板材一維暫態導熱比較；迎熱空氣瞬時達定值、迎熱h=8、背面h=3 W/(m²K)、背面參考空氣25°C。不包含廚房空氣從開火起的升溫、直接輻射、橫向導熱或完整熱浮力耦合。

物性與其他假設見deliverables/material_summary.json。計算使用numpy、scipy、matplotlib（Noto CJK字體）；報告另需markdown-it-py、playwright、Chromium。PDF與CSV、圖均在deliverables。

時間是達各面95%最終溫升，不是95%的°C數值，也不是廚房實際開火後達穩時間。切勿將原廚房模型熱點直接當作此板材真實板溫。
