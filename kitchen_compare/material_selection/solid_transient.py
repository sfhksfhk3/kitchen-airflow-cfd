"""Conditional 1D through-thickness transient solid FEM comparison.

This is NOT a transient kitchen CFD or a coupled wall-air solution. A step
in surrounding air temperature is imposed at t=0. No direct radiation or
lateral spreading is included. Constant properties and Robin films assumed.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS','1')
os.environ.setdefault('OMP_NUM_THREADS','1')
from pathlib import Path
import json,csv
import numpy as np
from scipy.linalg import eigh,solve
from scipy.optimize import brentq
ROOT=Path(__file__).resolve().parent
OUT=ROOT/'deliverables'
OUT.mkdir(exist_ok=True)
MATERIALS={
 'particleboard':dict(k=0.14,rho=600.,cp=1700.,thickness_m=.022,
                     source='https://cms.esi.info/Media/documents/155730_1334915100717.pdf'),
 'SUS304':dict(k=16.3,rho=7930.,cp=500.,thickness_m=.0015,
               source='https://www.nipponsteel.com/en/product/stainless/steel_shape/')}

def model(covered=False,hfront=8.,hback=3.,contact_R=0.,wood_elements=160):
 layers=(['SUS304'] if covered else [])+['particleboard']
 pts=[0.];elements=[];contact_pair=None;interface=0
 for li,key in enumerate(layers):
  mat=MATERIALS[key];ne=12 if key=='SUS304' else wood_elements
  if li>0:
   interface=len(pts)-1
   if contact_R>0:
    pts.append(pts[-1]);contact_pair=(len(pts)-2,len(pts)-1)
   interface=len(pts)-1
  first=len(pts)-1;x0=pts[-1]
  pts.extend((x0+np.linspace(0,mat['thickness_m'],ne+1)[1:]).tolist())
  for i in range(first,first+ne):elements.append((i,i+1,key))
 n=len(pts);M=np.zeros((n,n));K=M.copy()
 for i,j,key in elements:
  mat=MATERIALS[key];dx=pts[j]-pts[i]
  M[np.ix_([i,j],[i,j])]+=mat['rho']*mat['cp']*dx/6*np.array([[2,1],[1,2]])
  K[np.ix_([i,j],[i,j])]+=mat['k']/dx*np.array([[1,-1],[-1,1]])
 if contact_pair:
  K[np.ix_(contact_pair,contact_pair)]+=np.array([[1,-1],[-1,1]])/contact_R
 K[0,0]+=hfront;K[-1,-1]+=hback
 f=np.zeros(n);f[0]=hfront
 steady=solve(K,f,assume_a='pos')
 lam,V=eigh(K,M);coeff=V.T@M@steady
 assert lam.min()>0
 def field(t):return steady-V@(coeff*np.exp(-lam*t))
 def t95(node):
  return brentq(lambda t:field(t)[node]/steady[node]-.95,0.,72*3600)
 names={0:'exposed_surface',interface:'wood_hot_face',n-1:'wood_back_face'}
 # Retain explicit metrics even when bare-board surface and wood face coincide.
 result=dict(covered=covered,hfront_W_m2K=hfront,hback_W_m2K=hback,
             contact_R_m2K_W=contact_R,wood_elements=wood_elements,
             exposed_surface_t95_min=t95(0)/60,wood_hot_face_t95_min=t95(interface)/60,
             wood_back_face_t95_min=t95(n-1)/60,
             exposed_surface_steady_rise_fraction=float(steady[0]),
             wood_hot_face_steady_rise_fraction=float(steady[interface]),
             wood_back_face_steady_rise_fraction=float(steady[-1]),
             front_back_steady_conductance_W_m2K=float(hback*steady[-1]),
             heat_capacity_per_area_J_m2K=float(M.sum()))
 # Independently check the exact steady thermal resistance solution.
 Rsolid=sum(MATERIALS[k]['thickness_m']/MATERIALS[k]['k'] for k in layers)+contact_R
 q=0 if hback==0 else 1/(1/hfront+Rsolid+1/hback)
 exact_front=1-q/hfront;exact_back=1 if hback==0 else q/hback
 assert abs(steady[0]-exact_front)<1e-8 and abs(steady[-1]-exact_back)<1e-8
 times=np.r_[0,np.arange(30,6*3600+1,30)]
 history=np.array([[t,*field(t)[[0,interface,n-1]]] for t in times])
 return result,history,np.array(pts),steady

def main():
 baseline=[]
 for covered,label in [(False,'22mm_particleboard'),(True,'1p5mm_SUS304_on_22mm_particleboard')]:
  result,history,x,steady=model(covered)
  result['name']=label;baseline.append(result)
  np.savetxt(OUT/f'{label}_history.csv',history,delimiter=',',header='time_s,exposed_surface_normalized_rise,wood_hot_face_normalized_rise,wood_back_face_normalized_rise',comments='')
  np.savetxt(OUT/f'{label}_steady_profile.csv',np.c_[x,steady],delimiter=',',header='distance_from_exposed_face_m,normalized_temperature_rise',comments='')
 sensitivity=[]
 for covered in [False,True]:
  for hf in [5.,8.,15.]:
   for hb in [0.,1.,3.,5.]:
    r,*_=model(covered,hfront=hf,hback=hb,wood_elements=100);sensitivity.append(r)
 contact=[]
 for resistance in [0.,.001,.01,.1]:
  r,*_=model(True,contact_R=resistance);contact.append(r)
 mesh=[]
 for covered in [False,True]:
  for ne in [80,160,320]:
   r,*_=model(covered,wood_elements=ne);mesh.append(r)
 for filename,rows in [('material_baseline.csv',baseline),('film_sensitivity.csv',sensitivity),('contact_sensitivity.csv',contact),('mesh_validation.csv',mesh)]:
  keys=list(dict.fromkeys(k for row in rows for k in row))
  with (OUT/filename).open('w') as f:
   w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)
 summary=dict(materials=MATERIALS,assumptions=dict(initial_temperature_C=25,
       imposed_air_temperature='step to constant warm-air value; normalized, magnitude unspecified',
       ambient_back_C=25,baseline_hfront_W_m2K=8,baseline_hback_W_m2K=3,
       perfect_contact=True,no_direct_radiation=True,no_lateral_spreading=True,
       no_air_startup_or_buoyancy_recoupling=True,scope='conditional solid response, not actual kitchen startup prediction'),
       baseline=baseline,mesh=mesh,contact_sensitivity=contact,
       film_sensitivity_t95_ranges={str(c):[min(r['exposed_surface_t95_min'] for r in sensitivity if r['covered']==c),max(r['exposed_surface_t95_min'] for r in sensitivity if r['covered']==c)] for c in [False,True]})
 (OUT/'material_summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2))
 import matplotlib
 matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 plt.rcParams['font.family']='Noto Sans CJK JP'
 fig,axes=plt.subplots(1,2,figsize=(12,5))
 for covered,label,color in [(False,'22 mm顆粒板','#ad6200'),(True,'1.5 mm不鏽鋼覆22 mm顆粒板','#006eaa')]:
  r,H,x,S=model(covered)
  axes[0].plot(H[:,0]/60,H[:,1],label=label,color=color)
  axes[0].plot(H[:,0]/60,H[:,3],ls='--',color=color,alpha=.6)
  axes[1].plot(x*1000,S,label=label,color=color)
 axes[0].set(xlabel='從板材受熱起的時間（分鐘）',ylabel='溫升／熱空氣相對25°C的溫升',xlim=(0,240),ylim=(0,1),title='實線：迎熱面；虛線：顆粒板背面')
 axes[1].set(xlabel='從迎熱面起的深度（mm）',ylabel='穩態溫升／熱空氣溫升',ylim=(0,1),title='緊貼條件：兩層界面沒有溫度跳躍')
 for ax in axes:ax.grid(alpha=.2);ax.legend(fontsize=9)
 fig.suptitle('板材一維暫態有限元：條件比較，並非廚房從開火起的預測\n初始25°C；熱空氣瞬時達定值；迎熱面 h=8、背面 h=3 W/(m²·K)',fontsize=12)
 fig.tight_layout();fig.savefig(OUT/'material_transient_comparison.png',dpi=160);plt.close(fig)
 print(json.dumps(summary,ensure_ascii=False,indent=2))

if __name__=='__main__':main()
