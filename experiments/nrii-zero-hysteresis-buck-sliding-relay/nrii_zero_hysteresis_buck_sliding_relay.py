#!/usr/bin/env python3
"""Electrical switching-singularity benchmark: NRII smooth branches + semismooth relay."""
import csv, json, math, time
from pathlib import Path
import numpy as np

VIN=24.0; L=100e-6; C=100e-6; R=4.0; IREF=3.0
ORDER=40; TOL=1e-12; GAMMA=1.0


def A():
    return np.array([[0.0,-1.0/L],[1.0/C,-1.0/(R*C)]],float)

def b(u):
    return np.array([u*VIN/L,0.0],float)

def midpoint_res(y1,y0,dt,u):
    return y1-y0-dt*(A()@((y0+y1)*0.5)+b(u))

def direct_branch(y0,dt,u):
    i0,v0=y0
    M=np.array([[L/dt,.5],[-.5,C/dt+1/(2*R)]],float)
    rhs=np.array([L/dt*i0+u*VIN-.5*v0,C/dt*v0+.5*i0-v0/(2*R)],float)
    return np.linalg.solve(M,rhs)

def nrii_branch(y0,dt,u):
    M=A(); c=dt*(M@y0+b(u)); y=y0+c; tail=np.linalg.norm(c,np.inf); mx=tail
    for _ in range(2,ORDER+1):
        c=(dt*.5)*(M@c); y=y+c; tail=np.linalg.norm(c,np.inf); mx=max(mx,tail)
    r=np.linalg.norm(midpoint_res(y,y0,dt,u),np.inf)
    valid=(y[0]<IREF) if u==1 else (y[0]>IREF)
    d=direct_branch(y0,dt,u)
    return dict(u=int(u),candidate_i_A=float(y[0]),candidate_v_V=float(y[1]),residual_inf=float(r),
                tail_inf=float(tail),max_coeff_inf=float(mx),order=ORDER,domain_valid=bool(valid),
                direct_verifier_max_abs_error=float(np.max(np.abs(y-d))))

def binary_rj(x,y0,dt):
    i,v=x; i0,v0=y0; u=1 if i<=IREF else 0; vm=.5*(v0+v); im=.5*(i0+i)
    r=np.array([L/dt*(i-i0)-(u*VIN-vm),C/dt*(v-v0)-(im-vm/R)],float)
    J=np.array([[L/dt,.5],[-.5,C/dt+1/(2*R)]],float)
    return r,J,u

def strict_newton(y0,dt,n=8):
    x=y0.copy(); hist=[]; seen={}; cyc=None
    for k in range(n):
        r,J,u=binary_rj(x,y0,dt); rn=float(np.linalg.norm(r,np.inf))
        hist.append(dict(iteration=k,i_A=float(x[0]),v_V=float(x[1]),active_u=u,residual_inf=rn))
        key=(round(float(x[0]),13),round(float(x[1]),13),u)
        if key in seen and cyc is None: cyc=dict(first_seen=seen[key],repeat_at=k,period=k-seen[key])
        else: seen[key]=k
        x=x+np.linalg.solve(J,-r)
    r,_,u=binary_rj(x,y0,dt); rn=float(np.linalg.norm(r,np.inf))
    return dict(converged=rn<1e-10,cycle=cyc,history=hist,final_i_A=float(x[0]),final_v_V=float(x[1]),final_active_u=u,final_residual_inf=rn)

def gen_rj(x,y0,dt):
    i,v,d=x; i0,v0=y0; vm=.5*(v0+v); im=.5*(i0+i); s=i-IREF; q=d-GAMMA*s
    proj=min(1.0,max(0.0,q)); r=np.array([L/dt*(i-i0)-(d*VIN-vm),C/dt*(v-v0)-(im-vm/R),d-proj],float)
    dp=0.0 if (q<-1e-12 or q>1+1e-12) else 1.0
    J=np.array([[L/dt,.5,-VIN],[-.5,C/dt+1/(2*R),0.0],[GAMMA*dp,0.0,1-dp]],float)
    return r,J,q,proj

def semismooth(y0,dt,d0=1.0):
    x=np.array([y0[0],y0[1],d0],float); hist=[]
    for k in range(20):
        r,J,q,proj=gen_rj(x,y0,dt); rn=float(np.linalg.norm(r,np.inf))
        hist.append(dict(iteration=k,i_A=float(x[0]),v_V=float(x[1]),duty=float(x[2]),residual_inf=rn,projection_argument=q,projection_value=proj))
        if rn<=TOL: break
        try: dx=np.linalg.solve(J,-r)
        except np.linalg.LinAlgError: return dict(converged=False,reason='singular generalized Jacobian',history=hist)
        alpha=1.0
        for _ in range(24):
            xt=x+alpha*dx; rt,*_=gen_rj(xt,y0,dt)
            if np.linalg.norm(rt,np.inf)<=(1-1e-4*alpha)*rn+1e-15: x=xt; break
            alpha*=.5
        else: x=x+dx
    r,_,q,proj=gen_rj(x,y0,dt); rn=float(np.linalg.norm(r,np.inf)); i,v,d=map(float,x)
    mode='SLIDING' if 1e-10<d<1-1e-10 and abs(i-IREF)<=1e-9 else ('OFF' if d<=1e-10 else ('ON' if d>=1-1e-10 else 'UNRESOLVED'))
    return dict(converged=rn<=TOL,mode=mode,i_A=i,v_V=v,duty=d,residual_inf=rn,iterations_recorded=len(hist),newton_updates=max(0,len(hist)-1),history=hist)

def closed_sliding(y0,dt):
    i0,v0=y0; i1=IREF
    v1=(C/dt*v0+.5*(i0+i1)-v0/(2*R))/(C/dt+1/(2*R))
    d=(L/dt*(i1-i0)+.5*(v0+v1))/VIN
    return dict(valid=bool(0<=d<=1),i_A=float(i1),v_V=float(v1),duty=float(d))

def surface(v):
    on=(VIN-v)/L; off=-v/L
    return dict(di_on_A_per_s=on,di_off_A_per_s=off,sliding_condition=bool(on>0 and off<0))

def run_case(v0,dt):
    y0=np.array([IREF,v0],float); off=nrii_branch(y0,dt,0.0); on=nrii_branch(y0,dt,1.0); nw=strict_newton(y0,dt); ss=semismooth(y0,dt); cl=closed_sliding(y0,dt)
    if cl['valid'] and ss['converged'] and ss['mode']=='SLIDING': cl['relay_max_abs_error']=float(max(abs(ss['i_A']-cl['i_A']),abs(ss['v_V']-cl['v_V']),abs(ss['duty']-cl['duty'])))
    cert=not off['domain_valid'] and not on['domain_valid']
    return dict(initial=dict(i_A=IREF,v_V=v0,dt_s=dt),surface=surface(v0),nrii_off=off,nrii_on=on,strict_binary_newton=nw,projected_semismooth_relay=ss,
                closed_form_sliding_verifier=cl,no_binary_root_certificate=dict(certified=cert,reason='Both unique affine branch roots violate their own switching inequalities.' if cert else 'not certified'))

def hyst_freq(v,delta):
    up=(VIN-v)/L; dn=v/L
    return math.inf if delta<=0 or up<=0 or dn<=0 else 1.0/(2*delta/up+2*delta/dn)

def main():
    t=time.perf_counter(); out=Path(__file__).resolve().parent; primary=run_case(10.0,10e-6)
    dts=[1e-7,3e-7,1e-6,3e-6,1e-5,3e-5,5e-5]; ds=[]
    for dt in dts:
        c=run_case(10.0,dt); o,nw,ss=c['nrii_off'],c['strict_binary_newton'],c['projected_semismooth_relay']; on=c['nrii_on']
        ds.append(dict(dt_s=dt,off_i_A=o['candidate_i_A'],off_domain_valid=o['domain_valid'],on_i_A=on['candidate_i_A'],on_domain_valid=on['domain_valid'],
                       newton_converged=nw['converged'],newton_cycle_period=None if nw['cycle'] is None else nw['cycle']['period'],relay_converged=ss['converged'],relay_mode=ss['mode'],
                       relay_duty=ss['duty'],relay_i_A=ss['i_A'],relay_v_V=ss['v_V'],relay_residual_inf=ss['residual_inf'],relay_updates=ss['newton_updates'],nrii_max_residual_inf=max(o['residual_inf'],on['residual_inf'])))
    vs=[]
    for v0 in [-1.0,2.0,10.0,12.0,14.0,23.0,25.0]:
        c=run_case(v0,10e-6); ss=c['projected_semismooth_relay']; sf=c['surface']
        vs.append(dict(v0_V=v0,di_on_A_per_s=sf['di_on_A_per_s'],di_off_A_per_s=sf['di_off_A_per_s'],surface_sliding_condition=sf['sliding_condition'],relay_converged=ss['converged'],relay_mode=ss['mode'],relay_i_A=ss['i_A'],relay_v_V=ss['v_V'],relay_duty=ss['duty'],relay_residual_inf=ss['residual_inf']))
    hs=[dict(half_band_A=d,estimated_switching_Hz_at_v10=hyst_freq(10.0,d)) for d in [1e-1,1e-2,1e-3,1e-4,1e-5,1e-6]]
    results=dict(test_name='NRII Zero-Hysteresis Buck Sliding Relay',status='research benchmark; not production validation',parameters=dict(Vin=VIN,L=L,C=C,R=R,Iref=IREF),primary_case=primary,dt_sweep=ds,voltage_sweep=vs,hysteresis_limit=hs,elapsed_s=time.perf_counter()-t)
    (out/'nrii_zero_hysteresis_buck_sliding_relay_results.json').write_text(json.dumps(results,indent=2),encoding='utf-8')
    for name,rows in [('NRII_Zero_Hysteresis_Buck_Sliding_Relay_Sweep.csv',ds),('NRII_Zero_Hysteresis_Buck_Voltage_Gate.csv',vs),('NRII_Zero_Hysteresis_Buck_Hysteresis_Limit.csv',hs)]:
        with (out/name).open('w',newline='',encoding='utf-8') as f: w=csv.DictWriter(f,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    p=primary
    print(json.dumps(dict(test_name=results['test_name'],primary_sliding_condition=p['surface']['sliding_condition'],no_binary_root_certified=p['no_binary_root_certificate']['certified'],
                          strict_binary_newton_converged=p['strict_binary_newton']['converged'],strict_binary_newton_cycle=p['strict_binary_newton']['cycle'],relay_converged=p['projected_semismooth_relay']['converged'],
                          relay_mode=p['projected_semismooth_relay']['mode'],relay_duty=p['projected_semismooth_relay']['duty'],relay_i_A=p['projected_semismooth_relay']['i_A'],relay_v_V=p['projected_semismooth_relay']['v_V'],
                          relay_residual_inf=p['projected_semismooth_relay']['residual_inf'],relay_vs_closed_form_max_abs_error=p['closed_form_sliding_verifier'].get('relay_max_abs_error'),
                          nrii_vs_direct_max_abs_error=max(p['nrii_on']['direct_verifier_max_abs_error'],p['nrii_off']['direct_verifier_max_abs_error']),all_dt_cases_relay_converged=all(r['relay_converged'] for r in ds),
                          all_dt_cases_binary_newton_failed=all(not r['newton_converged'] for r in ds),elapsed_s=results['elapsed_s']),indent=2))
if __name__=='__main__': main()
