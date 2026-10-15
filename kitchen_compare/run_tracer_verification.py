from pathlib import Path
import argparse,json
from tracer import simulate
ROOT=Path(__file__).resolve().parent
if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('--layout',required=True);a.add_argument('--case',default='all');v=a.parse_args()
 for case in (['both_open','both_closed','left_open','right_open'] if v.case=='all' else [v.case]):
  simulate(v.layout,case,'_refined',dt=.25)
  for Sc in [.5,1.0]:simulate(v.layout,case,'_refined',Sc=Sc,dt=.25,suffix='_Sc'+str(Sc))
 # Worst supply case and one open case for temporal refinement.
 for case in ([c for c in ['both_open','both_closed'] if v.case in ['all',c]]):simulate(v.layout,case,'_refined',dt=.125,suffix='_dt0125')
