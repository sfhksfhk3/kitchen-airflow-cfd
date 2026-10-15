"""Wait on existing solves, repair only missing terminal fields, export deliverables."""
from pathlib import Path
import time,json,subprocess,datetime
from concurrent.futures import ThreadPoolExecutor
from run_batch import run_case
from foam_io import times
ROOT=Path(__file__).resolve().parents[1]

def phase(name):
    p=ROOT/'pipeline_progress.json';d=json.loads(p.read_text()) if p.exists() else {};d[name]=datetime.datetime.now(datetime.timezone.utc).isoformat();p.write_text(json.dumps(d,indent=2));print('PHASE',name,flush=True)
def run(script):
    subprocess.run(['python',str(ROOT/'scripts'/script)],check=True)

while True:
    p=ROOT/'batch_status.json';status=json.loads(p.read_text()) if p.exists() else []
    if len(status)==8:break
    time.sleep(5)
phase('primary_main_runs_returned')
with ThreadPoolExecutor(max_workers=2) as pool:
    out=list(pool.map(run_case,sorted((ROOT/'cases').iterdir())))
if not all(v['status'] in ['complete','already_complete'] for v in out):raise RuntimeError(out)
phase('terminal_fields_complete')
expected=['center_both_closed_dt','center_both_closed_grid','center_both_open_dt','corner_both_closed_grid']
while True:
    done=True
    for name in expected:
        p=ROOT/'sensitivity'/name
        state=p/'run_status.json'
        if state.exists() and json.loads(state.read_text())['status']=='failed':raise RuntimeError(f'Sensitivity failed: {name}')
        if not p.exists() or max(times(p),default=0)<139.99:done=False
    if done:break
    time.sleep(5)
phase('sensitivity_complete')
for name in ['verify_physics.py','process_results.py','validate_results.py','visualize_results.py','make_report.py','review_artifacts.py']:run(name)
phase('reports_and_data_built')
d=json.loads((ROOT/'PROGRESS.json').read_text());d.update({'overall_status':'COMPUTATION_AND_EXPORT_COMPLETE_PENDING_ARTIFACT_REVIEW','updated_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'next_tasks':['review browser/figures/report','package with SHA256 and finalize']});(ROOT/'PROGRESS.json').write_text(json.dumps(d,ensure_ascii=False,indent=2))
