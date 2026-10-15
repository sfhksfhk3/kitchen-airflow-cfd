"""Run independent OpenFOAM cases, resume latest saved time, keep a durable log."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
import argparse,subprocess,json,time,re,datetime
ROOT=Path(__file__).resolve().parents[1]

def latest(p):
    return max([float(x.name) for x in p.iterdir() if x.is_dir() and re.fullmatch(r'\d+(?:\.\d+)?',x.name)],default=0.)

def run_case(p):
    p=Path(p);now=latest(p)
    control=(p/'system/controlDict').read_text();end=float(re.search(r'endTime\s+([\d.]+)',control)[1])
    if now>=end-1e-8:return {'case':p.name,'status':'already_complete','last_time_s':now}
    stamp=datetime.datetime.now().strftime('%Y%m%d_%H%M%S');log=p/f'run_{stamp}.log'
    print('START',p.name,'from',now,'to',end,flush=True);start=time.time()
    with log.open('w') as f:
        rc=subprocess.call(['docker','exec','-u','1000:1000','kitchen-openfoam-2406','bash','-lc',f'/workspace/kitchen_openfoam/bin/kitchenAerosolFoam -case {p}'],stdout=f,stderr=subprocess.STDOUT)
    result={'case':p.name,'status':'complete' if rc==0 and latest(p)>=end-1e-8 else 'failed','returncode':rc,'last_time_s':latest(p),'runtime_wall_s':time.time()-start,'log':str(log)}
    (p/'run_status.json').write_text(json.dumps(result,indent=2));print('DONE',json.dumps(result),flush=True)
    return result

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--workers',type=int,default=3);ap.add_argument('cases',nargs='*');a=ap.parse_args()
    cases=[ROOT/'cases'/x for x in a.cases] if a.cases else sorted((ROOT/'cases').glob('*'))
    statuses=[]
    with ThreadPoolExecutor(max_workers=a.workers) as pool:
        for f in as_completed([pool.submit(run_case,p) for p in cases]):
            statuses.append(f.result());(ROOT/'batch_status.json').write_text(json.dumps(statuses,indent=2))
