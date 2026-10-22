"""Check final offline controls against exported CSV, plus PDF/figure integrity."""
from pathlib import Path
import json,re
import numpy as np
import pandas as pd
from PIL import Image
from pypdf import PdfReader
from playwright.sync_api import sync_playwright
from process_ceiling import DS,size_weights
from long_term import cycle_units

ROOT=Path(__file__).resolve().parents[1]
S=ROOT/'scenarios/ceiling33_Q30'

def main():
 figures=list((S/'figures').glob('*_1to1.png'))
 assert len(figures)==8,len(figures)
 for f in figures:
  im=Image.open(f);assert im.width==im.height,(f,im.size)
 pdf=PdfReader(S/'analysis_zh.pdf');text='\n'.join(p.extract_text() for p in pdf.pages)
 assert all(x in text for x in ['天花','1300','33','近壁','Satin']),text[:200]
 square=PdfReader(S/'figures/all_figures_1to1.pdf')
 assert len(square.pages)==8
 assert all(abs(float(p.mediabox.width)-float(p.mediabox.height))<.01 for p in square.pages)
 errors=[];requests=[];checks=0
 with sync_playwright() as p:
  b=p.chromium.launch(executable_path='/usr/bin/chromium',args=['--no-sandbox'])
  page=b.new_page(viewport={'width':1400,'height':1100},device_scale_factor=1)
  page.on('pageerror',lambda e:errors.append(str(e)))
  page.on('request',lambda r:requests.append(r.url) if r.url.startswith(('http://','https://')) else None)
  page.set_content((S/'interactive_3d.html').read_text(),wait_until='load')
  page.wait_for_function("document.querySelector('#plot').data?.length===3")
  cases=page.locator('#case option').evaluate_all('(xs)=>xs.map(x=>x.value)');assert len(cases)==12,len(cases)
  def check(case,gsd,dose=1,mult=1):
   nonlocal checks
   f=pd.read_csv(S/'results'/case/'ceiling_complete_event_faces.csv')
   vals=f[[f'M_{d:g}um_unit_mg_m2_per_mg' for d in DS]].to_numpy()@size_weights(gsd)*dose*mult
   expected=float(vals@f.area_m2)
   stats=page.locator('#stats').inner_text();actual=float(re.search(r'總沉積 ([\deE+.-]+) mg',stats)[1])
   assert np.isclose(actual,expected,rtol=5.1e-4,atol=1e-300),(case,gsd,expected,stats)
   checks+=1
  for case in cases:
   page.select_option('#case',case)
   for gsd in [1.2,1.5,2.0]:
    page.select_option('#gsd',str(gsd));check(case,gsd)
  case='corner_both_closed_T33_Q30';page.select_option('#case',case);page.select_option('#gsd','1.5')
  for selector,value in [('#dose','7.2'),('#days','365')]:page.fill(selector,value);page.dispatch_event(selector,'change')
  page.select_option('#events','2');page.select_option('#clean','30')
  check(case,1.5,7.2,cycle_units(365,2,30,.9)[0])
  page.fill('#dose','0');page.dispatch_event('#dose','change');check(case,1.5,0)
  page.click('#angle');page.wait_for_function("document.querySelector('#plot').layout.scene.camera.projection.type==='perspective'")
  page.click('#top');page.wait_for_function("document.querySelector('#plot').layout.scene.camera.projection.type==='orthographic'")
  page.fill('#dose','1');page.dispatch_event('#dose','change');page.fill('#days','1');page.dispatch_event('#days','change');page.select_option('#events','1');page.select_option('#clean','0')
  page.wait_for_timeout(300)
  page.screenshot(path=str(S/'interactive_preview.png'),full_page=True)
  b.close()
 assert not errors,errors
 assert not requests,requests
 result=dict(status='PASS',ui_csv_comparisons=checks,offline_external_http_requests=requests,javascript_errors=errors,pdf_pages=len(pdf.pages),square_figures=len(figures),square_pdf_pages=len(square.pages),scope='Presentation and output consistency only; not physical validation.')
 (S/'PRESENTATION_VALIDATION.json').write_text(json.dumps(result,ensure_ascii=False,indent=2));print(json.dumps(result,ensure_ascii=False,indent=2))

if __name__=='__main__':main()
