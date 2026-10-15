"""Frozen-flow thermal sensitivity; not a newly coupled CFD solution."""
from pathlib import Path
import sys,json,numpy as np
from skfem import MeshTet,Basis,FacetBasis,ElementTetP1,asm
from scipy.sparse.linalg import spsolve
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parent/'corner'))
import solve

def main():
 p=ROOT.parent/'corner/results/both_closed_refined.npz';d=np.load(p)
 m=MeshTet(d['points'].T,d['cells'].T);bs=Basis(m,ElementTetP1(),intorder=3);f=solve.faces(m,'both_closed');coords=bs.global_coordinates();source=np.zeros_like(coords[0]);P=solve.P
 for bx,by,bz in solve.CFG['burners']:
  s=np.exp(-((coords[0]-bx)**2+(coords[1]-by)**2)/(2*.085**2)-((coords[2]-bz-.07)**2)/(2*.055**2));s*=((coords[2]>bz)&(coords[2]<bz+.25)&(coords[1]<solve.CD));source+=s*2000/np.sum(s*bs.dx)
 fbw=FacetBasis(m,ElementTetP1(),facets=f['walls'],intorder=3);fbt=FacetBasis(m,ElementTetP1(),facets=np.r_[f['hood'],f['left'],f['right']],intorder=3)
 wmat=asm(solve.mass,fbw)*3/(P['rho']*P['cp']);rhs=asm(solve.positive_heat_load,bs,source=source/(P['rho']*P['cp']));dn=np.unique(np.r_[m.facets[:,f['left']].ravel(),m.facets[:,f['right']].ravel()]);free=np.setdiff1d(np.arange(bs.N),dn)
 cached=ROOT/'deliverables/hotspot_RT0_input.npz'
 tracer=np.load(cached if cached.exists() else ROOT.parent/'corner/results/both_closed_refined_tracer.npz');urt=tracer['velocity_RT0_cells'];vel=np.repeat(urt.T[:,:,None],bs.dx.shape[1],axis=2);vn=np.sum(urt[fbt.tind].T[:,:,None]*fbt.normals,axis=0)
 A=solve.monotone_heat(asm(solve.heat_galerkin,bs,vel=vel,nu=P['alpha_molecular']+P['nu_t_assumed']/P['Pr_t'])+asm(solve.heat_boundary,fbt,vn=vn)+wmat)
 theta=np.zeros(bs.N);theta[free]=spsolve(A[free][:,free].tocsc(),rhs[free]);T=theta+25
 wallfac=m.boundary_facets();v=m.p[:,m.facets[:,wallfac]];c=v.mean(axis=1);selected=wallfac[(np.abs(c[0]-.7)<1e-8)&(c[1]<.60000001)];nodes=np.unique(m.facets[:,selected]);imax=nodes[T[nodes].argmax()];maxorig=nodes[d['temperature'][nodes].argmax()]
 fbside=FacetBasis(m,ElementTetP1(),facets=selected,intorder=3);area=float(fbside.dx.sum());power=float(3*np.sum(fbside.interpolate(theta)*fbside.dx))
 result={'purpose':'frozen-flow thermal diagnostic only; RT0 velocity from existing passive-tracer calculation; buoyancy NOT recoupled', 'original_wall_max_C':float(d['temperature'][maxorig]),'original_wall_max_xyz_m':m.p[:,maxorig].tolist(),'RT0_wall_max_C':float(T[imax]),'RT0_wall_max_xyz_m':m.p[:,imax].tolist(),'RT0_air_max_C':float(T.max()),'RT0_wall_heat_W':power,'RT0_wall_average_boundary_C':25+power/(3*area),'RT0_thermal_residual':float(np.linalg.norm((A@theta-rhs)[free])/np.linalg.norm(rhs[free]))}
 (ROOT/'deliverables/hotspot_diagnostic.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
 # Persist just the fixed-flow diagnostic input, avoiding a large tracer dependency.
 np.savez_compressed(ROOT/'deliverables/hotspot_RT0_input.npz',velocity_RT0_cells=urt)
 import matplotlib
 matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 import matplotlib.tri as mtri
 plt.rcParams['font.family']='Noto Sans CJK JP'
 _,local=np.unique(m.facets[:,selected].T,return_inverse=True);local=local.reshape(-1,3)
 triang=mtri.Triangulation(m.p[1,nodes],m.p[2,nodes],local)
 fig,axes=plt.subplots(1,2,figsize=(10,10),sharey=True)
 for ax,values,peak,title in zip(axes,[d['temperature'],T],[maxorig,imax],['原始已存模型','局部守恆修正後的熱場敏感性檢查']):
  plot=ax.tripcolor(triang,values[nodes],shading='gouraud',cmap='inferno',vmin=25,vmax=133)
  ax.plot(.3,.9,marker='x',color='limegreen',ms=13,mew=3,clip_on=False)
  ax.axvline(.3,ls=':',color='limegreen',lw=1,alpha=.5)
  ax.plot(m.p[1,peak],m.p[2,peak],marker='x',color='cyan',ms=10,mew=2)
  ax.set_aspect('equal');ax.set_xlabel('距後牆 y（m）');ax.set_title(title,fontsize=12)
  ax.set_xlim(-.035,.635);ax.set_ylim(.85,2.23)
  ax.annotate(f'{values[peak]:.1f}°C',xy=(m.p[1,peak],m.p[2,peak]),xytext=(.04,1.42),color='deepskyblue',arrowprops=dict(arrowstyle='->',color='deepskyblue'))
 axes[0].set_ylabel('高度 z（m）')
 fig.suptitle('爐頭位置與側面熱點是兩個不同的位置\n兩圖都是 x=0.70 m 側面；綠色 × 是左爐中心投影',fontsize=14)
 fig.subplots_adjust(top=.88,bottom=.16,left=.09,right=.86,wspace=.20)
 fig.colorbar(plot,cax=fig.add_axes([.89,.20,.025,.62]),label='模型邊界溫度（°C）')
 fig.text(.09,.055,'左爐中心：(x,y,z)=(0.95,0.30,0.90) m，距側面 0.25 m（沿圖外 x 方向）\n青色 ×：各數值處理下的側面最高值；其位置與數值會改變\n右圖保持既有流場的守恆修正速度，不是重新耦合求解結果；兩圖均未計板材導熱／直接輻射',fontsize=11)
 fig.savefig(ROOT/'deliverables/wall_hotspot_explanation.png',dpi=160);plt.close(fig)

if __name__=='__main__':main()
