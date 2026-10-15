"""Independent checks on sizes, dilution, wall quadrature and histories."""
from pathlib import Path
import json,numpy as np
from scipy.integrate import quad
from emission_scenario import config,DIAMETERS
ROOT=Path(__file__).resolve().parents[1]

def properties(d_um,T=298.15):
    d=d_um*1e-6;p=101325;R=287.05;rho=p/(R*T);mu=1.458e-6*T**1.5/(T+110.4)
    lam=mu/p*np.sqrt(np.pi*R*T/2);Kn=2*lam/d;Cc=1+Kn*(1.257+.4*np.exp(-1.1/Kn))
    return rho,mu,1.380649e-23*T*Cc/(3*np.pi*mu*d),(850-rho)*9.81*d*d*Cc/(18*mu)

def resistance(D,dist,ustar,nu,n):
    y=np.r_[0,dist*10**np.linspace(-9+9/n,0,n)]
    diffusivity=D+.41*ustar*y*(-np.expm1(-y*ustar/(26*nu)))**2/.7
    return np.trapezoid(1/diffusivity,y)

def main():
    rows=[];errs=[]
    for d in DIAMETERS:
        rho,mu,D,vs=properties(d);rows.append({'diameter_um':float(d),'brownian_D_m2_s':D,'settling_velocity_m_s':vs,'settling_1m_time_s':1/vs,'response_time_s':vs/9.81})
        for T in [298.15,450.]:
            rho,mu,D,_=properties(d,T);nu=mu/rho
            for dist in [.00125,.025,.075]:
                for ustar in [.01,.1,1.]:
                    # Independent adaptive quadrature over log distance, plus
                    # the vanishingly small [0, 1e-12*dist] analytic segment.
                    tiny=dist*1e-12
                    fun=lambda s:np.exp(s)/(D+.41*ustar*np.exp(s)*(-np.expm1(-np.exp(s)*ustar/(26*nu)))**2/.7)
                    ref=tiny/D+quad(fun,np.log(tiny),np.log(dist),epsrel=1e-10)[0]
                    a=resistance(D,dist,ustar,nu,72);b=resistance(D,dist,ustar,nu,144)
                    errs.append({'diameter_um':float(d),'distance_m':dist,'ustar_m_s':ustar,'temperature_K':T,'vdep72_relative_error':ref/a-1,'vdep144_relative_error':ref/b-1})
    c=config();Rkg=c['total_rate_ug_s']*1e-9
    out={'particle_properties_25C':rows,'wall_quadrature_max_abs_relative_error_72':max(abs(x['vdep72_relative_error']) for x in errs),'wall_quadrature_max_abs_relative_error_144':max(abs(x['vdep144_relative_error']) for x in errs),'wall_quadrature_details':errs,'bulk_particle_to_air_mass_flow_ratio_at_0p5m3_s':Rkg/(.5*properties(1)[0]),'interpretation':'quadrature and constitutive scaling checks; not validation against cooking experiments'}
    (ROOT/'verification/particle_physics_checks.json').write_text(json.dumps(out,indent=2))
    print({k:v for k,v in out.items() if k not in ['wall_quadrature_details']})

if __name__=='__main__':main()
