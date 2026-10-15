"""Check completed exports and exercise the offline 3D viewer in Chromium."""
from pathlib import Path
import json,csv,math
from PIL import Image
from pypdf import PdfReader
from playwright.sync_api import sync_playwright

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'deliverables'

def main():
    summary=json.loads((OUT/'summary.json').read_text())
    assert len(summary['cases'])==8 and len(summary['size_cases'])==40
    closure={}
    for c in summary['cases']:
        total=sum(c[k] for k in ['hood_fraction_at_end','door_escape_fraction_end','deposited_fraction_end','remaining_fraction_end'])
        closure[c['case']]=abs(total-1)
        assert abs(total-1)<1e-6,(c['case'],total)
        p=OUT/'data'/c['case']
        rows=list(csv.DictReader((p/'emission_distribution_sensitivity.csv').open()))
        assert len(rows)==9
        for g in [1.2,1.5,2.]:
            values={float(r['emission_rate_factor']):float(r['walkway_90_120_oil_ug_m3']) for r in rows if float(r['assumed_mass_GSD'])==g}
            assert math.isclose(values[4.],16*values[.25],rel_tol=1e-12)
    figures=[]
    for p in sorted((OUT/'figures_1to1').glob('*.png')):
        size=Image.open(p).size;assert size==(3600,3600),(p,size)
        pdf=PdfReader(p.with_suffix('.pdf'));box=pdf.pages[0].mediabox
        assert abs(float(box.width)-float(box.height))<.01
        figures.append({'file':p.name,'pixels':size})
    assert len(figures)==9
    report=PdfReader(OUT/'openfoam_oil_analysis_zh.pdf')
    texts=[p.extract_text() for p in report.pages]
    assert '廚房' in ''.join(texts) and '4000' in ''.join(texts)
    assert all(t.strip() for t in texts),'Blank PDF page'
    assert len(PdfReader(OUT/'all_figures_1to1.pdf').pages)==9
    errors=[];network=[]
    with sync_playwright() as p:
        b=p.chromium.launch(executable_path='/usr/bin/chromium',headless=True,args=['--no-sandbox','--enable-unsafe-swiftshader','--use-angle=swiftshader'])
        tab=b.new_page(viewport={'width':1440,'height':1120})
        tab.on('pageerror',lambda e:errors.append(str(e)))
        tab.on('request',lambda req:network.append(req.url) if req.url.startswith(('http:','https:')) else None)
        tab.set_content((OUT/'kitchen_oil_openfoam_interactive_3d.html').read_text(),wait_until='load',timeout=120000)
        tab.wait_for_function("document.getElementById('corner')._fullLayout && document.getElementById('center')._fullLayout",timeout=60000)
        assert tab.evaluate('Object.keys(DATA.cases).length')==8
        for door in ['both_open','both_closed','left_open','right_open']:
            tab.select_option('#door',door)
            for field in ['oil','temperature','velocity','pressure']:
                tab.select_option('#field',field)
                tab.wait_for_timeout(250)
                result=tab.evaluate("['corner','center'].map(l=>({traces:document.getElementById(l).data.length,points:document.getElementById(l).data.find(t=>t.customdata)?.x.length||0}))")
                assert all(x['traces']>=9 and x['points']>0 for x in result),(door,field,result)
                if field=='oil':
                    assert tab.evaluate("document.getElementById('corner').data.find(t=>t.customdata).marker.colorscale[0][1]")=='#000004'
        tab.select_option('#door','both_closed');tab.select_option('#field','oil')
        for size in ['0','1','2','3','4','mix']:
            tab.select_option('#size',size);tab.wait_for_timeout(150)
        tab.fill('#gsd','2');tab.fill('#rate','57.6');tab.wait_for_timeout(300)
        assert 'NaN' not in tab.locator('#cornerMetrics').inner_text()
        tab.fill('#rate','14.4');tab.fill('#gsd','1.5')
        for tick in ['0','5','9']:
            tab.locator('#time').evaluate('(el,v)=>{el.value=v;el.dispatchEvent(new Event("input",{bubbles:true}))}',tick);tab.wait_for_timeout(300)
            assert f"{[0,10,30,60,100,120,130,150,180,240][int(tick)]} s"==tab.locator('#timeLabel').inner_text()
        tab.locator('#time').evaluate('(el)=>{el.value="5";el.dispatchEvent(new Event("input",{bubbles:true}))}');tab.wait_for_timeout(800)
        tab.screenshot(path=str(ROOT/'verification/final_viewer_oil.png'),full_page=True)
        tab.select_option('#field','velocity');tab.wait_for_timeout(800)
        tab.screenshot(path=str(ROOT/'verification/final_viewer_velocity.png'),full_page=True)
        tab.click('#reset');tab.wait_for_timeout(200)
        b.close()
    assert not errors,errors
    assert not network,network
    result={'summary_case_count':8,'independent_size_case_count':40,'mixture_mass_fraction_closure':closure,'square_figures':figures,'report_pages':len(report.pages),'pdf_chinese_text_extractable':True,'offline_viewer_cases':8,'viewer_door_field_combinations_checked':16,'viewer_sizes_checked':6,'viewer_time_ticks_checked':[0,120,240],'javascript_errors':errors,'external_network_requests':network,'note':'Automated structure and interaction checks; screenshots and selected scientific figures also need visual inspection.'}
    (ROOT/'verification/artifact_review.json').write_text(json.dumps(result,ensure_ascii=False,indent=2))
    print(json.dumps(result,ensure_ascii=False,indent=2))

if __name__=='__main__':main()
