from pathlib import Path
import json,numpy as np,pandas as pd
from playwright.sync_api import sync_playwright
from pypdf import PdfReader
from PIL import Image
ROOT=Path(__file__).resolve().parents[1]
errors=[];requests=[]
with sync_playwright() as p:
 b=p.chromium.launch(executable_path='/usr/bin/chromium',headless=True,args=['--no-sandbox','--enable-webgl','--use-gl=angle','--use-angle=swiftshader','--enable-unsafe-swiftshader'])
 page=b.new_page(viewport={'width':1440,'height':1050});page.on('pageerror',lambda e:errors.append(str(e)));page.on('request',lambda r:requests.append(r.url) if r.url.startswith(('https://','http://')) else None)
 page.set_content((ROOT/'ceiling_deposition_interactive_3d.html').read_text(),wait_until='load');page.wait_for_function("document.querySelector('#plot').data && document.querySelector('#plot').data.length===3")
 cases=page.locator('#case option').count();assert cases==12,cases
 assert page.evaluate("Array.isArray(document.querySelector('#plot').data[0].colorscale) && document.querySelector('#plot').data[0].colorscale[0][1]==='#000004'")
 for idx in range(cases):
  page.locator('#case').select_option(index=idx)
  page.wait_for_timeout(100)
  for gsd in ['1.2','1.5','2.0']:
   page.locator('#gsd').select_option(gsd)
   size=page.evaluate("document.querySelector('#plot').data[0].intensity.length")
   assert size>500
 page.locator('#case').select_option('corner_both_closed_parallel');page.locator('#gsd').select_option('1.5');page.locator('#dose').fill('28.8');page.locator('#dose').dispatch_event('change');page.locator('#days').fill('365');page.locator('#days').dispatch_event('change');page.locator('#events').select_option('2');page.locator('#clean').select_option('30')
 # Compare displayed per-face values to the long-term CSV independently.
 value=page.evaluate("Math.max(...document.querySelector('#plot').data[2].customdata.map(x=>x[1]))")
 d=pd.read_csv(ROOT/'results/long_term_scenarios.csv');row=d[(d['case']=='corner_both_closed_parallel')&(d.material=='gypsum_matte')&(d.days==365)&(d.events_per_day==2)&(d.emitted_mg_per_event==28.8)&(d.clean_efficiency==.9)].iloc[0]
 assert np.isclose(value,row.retained_max_mg_m2,rtol=1e-12),(value,row.retained_max_mg_m2)
 page.locator('#dose').fill('0');page.locator('#dose').dispatch_event('change');assert page.evaluate("Math.max(...document.querySelector('#plot').data[2].customdata.map(x=>x[1]))")==0
 page.locator('#dose').fill('1');page.locator('#dose').dispatch_event('change');page.locator('#days').fill('1');page.locator('#days').dispatch_event('change');page.locator('#events').select_option('1');page.locator('#clean').select_option('0');page.wait_for_timeout(500);page.screenshot(path=str(ROOT/'verification/interactive_screenshot.png'),full_page=True)
 b.close()
assert not errors,errors
assert not requests,requests
pdf=PdfReader(ROOT/'reports/ceiling_deposition_analysis_zh.pdf');text=''.join(p.extract_text() for p in pdf.pages)
assert '27' in text[:3000] and '附錄B' in text and len(pdf.pages)>=10
figs=list((ROOT/'figures').glob('*_1to1.png'))
for f in figs:
 w,h=Image.open(f).size;assert w==h,(f,w,h)
l=pd.read_csv(ROOT/'results/long_term_scenarios.csv');s=pd.read_csv(ROOT/'results/primary_complete_event_summary.csv').set_index('case')
# Cleaning changes retained inventory, never creates or destroys unaccounted grease.
for c in s.index:
 q=l[l['case']==c];expected=s.loc[c,'M_total_unit_mg_per_mg']*q.emitted_mg_per_event*q.days*q.events_per_day
 assert np.allclose(q.retained_total_mg+q.removed_total_mg,expected,rtol=1e-12,atol=1e-40)
result={'offline_browser_cases':cases,'size_distributions_each_case':3,'external_HTTP_requests':requests,'javascript_errors':errors,'interactive_mass_matches_CSV':True,'zero_emission_returns_zero':True,'long_term_cleaning_mass_balance':'PASS','PDF_pages':len(pdf.pages),'report_front_page_nonconvergence_warning':True,'square_figures':len(figs)}
(ROOT/'verification/deliverable_review.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
