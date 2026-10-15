"""Keep startup steps small, allow implicit quasi-steady gap-jet time steps.
This change is recorded and must be accompanied by a smaller-step sensitivity run.
"""
from pathlib import Path
import re,time,json,datetime
ROOT=Path(__file__).resolve().parents[1]
done=set();events=[]
while len(done)<2:
    for layout in ['corner','center']:
        if layout in done:continue
        p=ROOT/'cases'/f'{layout}_both_closed'
        files=list((p/'diagnostics').glob('*.csv')) if (p/'diagnostics').exists() else []
        if not files:continue
        f=max(files,key=lambda q:q.stat().st_mtime)
        with f.open('rb') as stream:
            stream.seek(max(0,f.stat().st_size-5000));lines=stream.read().decode().splitlines()
        good=[s for s in lines if re.match(r'^\d',s) and len(s.split(','))==25]
        if not good:continue
        t=float(good[-1].split(',')[0])
        if t<10:continue
        c=p/'system/controlDict';s=c.read_text();s=s.replace('maxCo 3;','maxCo 12;').replace('maxDeltaT 0.05;','maxDeltaT 0.04;');c.write_text(s)
        event={'case':p.name,'physical_time_s':t,'maxCo':12,'maxDeltaT_s':.04,'recorded_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'requirement':'Validate against maxCo3 and maxDeltaT0.02 on matched saved state'}
        events.append(event);done.add(layout);(ROOT/'runtime_timestep_events.json').write_text(json.dumps(events,indent=2));print(event,flush=True)
    time.sleep(3)
