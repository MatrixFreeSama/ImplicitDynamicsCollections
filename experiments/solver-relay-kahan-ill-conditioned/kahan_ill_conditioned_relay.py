#!/usr/bin/env python3
import time, json, csv, os
import numpy as np
import scipy.linalg as la
import mpmath as mp

N = 128
THETA = '1.2'
MP_ORACLE_DPS = 100
TARGET_FWD = 1e-12
OUT = os.path.dirname(os.path.abspath(__file__))


def kahan_mp(n, theta_str):
    th = mp.mpf(theta_str)
    c = mp.cos(th); s = mp.sin(th)
    A = mp.matrix(n, n)
    for i in range(n):
        si = s**i
        for j in range(n):
            if j < i:
                A[i,j] = mp.mpf('0')
            elif j == i:
                A[i,j] = si
            else:
                A[i,j] = -c*si
    return A


def x_true_np(n):
    return np.array([((-1.0)**i)*(1.0+((i*37)%11)/10.0) for i in range(n)], dtype=np.float64)


def vec_norm_mp(v):
    return mp.sqrt(mp.fsum(x*x for x in v))


def matvec_mp(A, v):
    n = A.rows
    return [mp.fsum(A[i,j]*v[j] for j in range(n)) for i in range(n)]


def matTvec_mp(A, v):
    n = A.rows
    return [mp.fsum(A[i,j]*v[i] for i in range(n)) for j in range(n)]


def upper_solve_mp(A, b):
    n = A.rows
    x = [mp.mpf('0')]*n
    for i in range(n-1, -1, -1):
        s = mp.fsum(A[i,j]*x[j] for j in range(i+1, n))
        x[i] = (b[i]-s)/A[i,i]
    return x


def lower_solve_AT_mp(A, b):
    n = A.rows
    x = [mp.mpf('0')]*n
    for i in range(n):
        s = mp.fsum(A[j,i]*x[j] for j in range(i))
        x[i] = (b[i]-s)/A[i,i]
    return x


def estimate_cond2_mp(A, iters=80):
    n = A.rows
    v = [mp.mpf(1)/mp.sqrt(n)]*n
    for _ in range(iters):
        Av = matvec_mp(A, v)
        w = matTvec_mp(A, Av)
        nw = vec_norm_mp(w)
        v = [x/nw for x in w]
    sigma_max = vec_norm_mp(matvec_mp(A, v))

    v = [mp.mpf(1)/mp.sqrt(n)]*n
    for _ in range(iters):
        y = upper_solve_mp(A, v)
        w = lower_solve_AT_mp(A, y)
        nw = vec_norm_mp(w)
        v = [x/nw for x in w]
    sigma_inv_max = vec_norm_mp(upper_solve_mp(A, v))
    sigma_min = 1/sigma_inv_max
    return sigma_max, sigma_min, sigma_max/sigma_min


def backward_error(A, x, b):
    r = A @ x - b
    return np.linalg.norm(r, 2)/(np.linalg.norm(A,2)*np.linalg.norm(x,2)+np.linalg.norm(b,2))


def forward_error(x, xt):
    return np.linalg.norm(x-xt,2)/np.linalg.norm(xt,2)


def build_oracle(dps=MP_ORACLE_DPS):
    mp.mp.dps = dps
    A_mp = kahan_mp(N, THETA)
    xt = x_true_np(N)
    xt_mp = mp.matrix([mp.mpf(str(v)) for v in xt])
    b_mp = A_mp * xt_mp
    A64 = np.array([[float(A_mp[i,j]) for j in range(N)] for i in range(N)], dtype=np.float64)
    b64 = np.array([float(b_mp[i]) for i in range(N)], dtype=np.float64)
    return A_mp, xt, xt_mp, b_mp, A64, b64


def run():
    A_mp, xt, xt_mp, b_mp, A, b = build_oracle()
    sigma_max_mp, sigma_min_mp, cond_mp = estimate_cond2_mp(A_mp)

    baselines = []
    t0 = time.perf_counter(); lu, piv = la.lu_factor(A); x_lu = la.lu_solve((lu,piv), b); t_lu=(time.perf_counter()-t0)*1e3
    baselines.append(('double_LU', forward_error(x_lu,xt), backward_error(A,x_lu,b), N, '', t_lu))

    t0=time.perf_counter(); Q,R=la.qr(A,mode='economic'); x_qr=la.solve_triangular(R,Q.T@b); t_qr=(time.perf_counter()-t0)*1e3
    baselines.append(('double_QR', forward_error(x_qr,xt), backward_error(A,x_qr,b), N, '', t_qr))

    t0=time.perf_counter(); Qc,Rc,p=la.qr(A,pivoting=True,mode='economic'); y=la.solve_triangular(Rc,Qc.T@b); x_cp=np.empty_like(y); x_cp[p]=y; t_cp=(time.perf_counter()-t0)*1e3
    cp_tol = np.finfo(np.float64).eps*N*np.max(np.abs(np.diag(Rc)))
    cp_rank = int(np.sum(np.abs(np.diag(Rc)) > cp_tol))
    baselines.append(('double_CPQR', forward_error(x_cp,xt), backward_error(A,x_cp,b), cp_rank, float(np.min(np.abs(np.diag(Rc)))), t_cp))

    t0=time.perf_counter(); x_svd,_,rank_svd,svals=np.linalg.lstsq(A,b,rcond=None); t_svd=(time.perf_counter()-t0)*1e3
    baselines.append(('double_SVD_lstsq', forward_error(x_svd,xt), backward_error(A,x_svd,b), int(rank_svd), float(svals[-1]), t_svd))

    x = x_lu.copy()
    double_ref=[]
    for k in range(6):
        double_ref.append((k,forward_error(x,xt),backward_error(A,x,b)))
        r=b-A@x
        dx=la.lu_solve((lu,piv),r)
        x=x+dx

    relay=[]
    for dps in [18,20,22,24,26,30,40,60,80]:
        mp.mp.dps=dps
        Ahi=kahan_mp(N,THETA)
        xthi=mp.matrix([mp.mpf(str(v)) for v in xt])
        bhi=Ahi*xthi
        x0mp=mp.matrix([mp.mpf(float(v)) for v in x_lu])
        rhi=bhi-Ahi*x0mp
        r64=np.array([float(rhi[i]) for i in range(N)])
        t0=time.perf_counter(); dx=la.lu_solve((lu,piv),r64); xr=x_lu+dx; elapsed=(time.perf_counter()-t0)*1e3
        relay.append((dps,forward_error(xr,xt),backward_error(A,xr,b),np.linalg.norm(dx,2)/np.linalg.norm(xr,2),elapsed))

    mp_solve=[]
    for dps in [20,25,30,40,50]:
        mp.mp.dps=dps
        Ahi=kahan_mp(N,THETA)
        xthi=mp.matrix([mp.mpf(str(v)) for v in xt])
        bhi=Ahi*xthi
        t0=time.perf_counter(); xhi=mp.lu_solve(Ahi,bhi); elapsed=(time.perf_counter()-t0)*1e3
        rel=vec_norm_mp([xhi[i]-xthi[i] for i in range(N)])/vec_norm_mp([xthi[i] for i in range(N)])
        mp_solve.append((dps,mp.nstr(rel,20),elapsed))

    with open(os.path.join(OUT,'BASELINES.csv'),'w',newline='') as f:
        w=csv.writer(f); w.writerow(['method','forward_error','backward_error','reported_rank','rank_indicator','time_ms'])
        w.writerows(baselines)
    with open(os.path.join(OUT,'DOUBLE_REFINEMENT.csv'),'w',newline='') as f:
        w=csv.writer(f); w.writerow(['iteration','forward_error','backward_error']); w.writerows(double_ref)
    with open(os.path.join(OUT,'PRECISION_RELAY.csv'),'w',newline='') as f:
        w=csv.writer(f); w.writerow(['residual_dps','forward_error_after_one_correction','backward_error_after_one_correction','relative_correction_norm','double_correction_time_ms']); w.writerows(relay)
    with open(os.path.join(OUT,'MULTIPRECISION_DIRECT.csv'),'w',newline='') as f:
        w=csv.writer(f); w.writerow(['working_dps','forward_error','time_ms']); w.writerows(mp_solve)
    with open(os.path.join(OUT,'MATRIX_DIAGNOSTICS.csv'),'w',newline='') as f:
        w=csv.writer(f); w.writerow(['quantity','value'])
        w.writerow(['n',N]); w.writerow(['theta',THETA]); w.writerow(['pert',0])
        w.writerow(['sigma_max_mp_power',mp.nstr(sigma_max_mp,30)]); w.writerow(['sigma_min_mp_inverse_power',mp.nstr(sigma_min_mp,30)]); w.writerow(['cond2_mp_estimate',mp.nstr(cond_mp,30)])
        w.writerow(['cpqr_min_abs_Rdiag',float(np.min(np.abs(np.diag(Rc))))]); w.writerow(['cpqr_rank_tol',cp_tol]); w.writerow(['cpqr_reported_rank',cp_rank]); w.writerow(['svd_reported_rank',int(rank_svd)]); w.writerow(['svd_smallest_sigma_double',float(svals[-1])])
        w.writerow(['cpqr_Rdiag_to_sigma_min_mp_ratio',mp.nstr(mp.mpf(float(np.min(np.abs(np.diag(Rc)))))/sigma_min_mp,30)])

    summary = {
        'benchmark':'Kahan-128 ill-conditioned linear system', 'n':N, 'theta':THETA, 'pert':0,
        'cond2_mp_estimate': mp.nstr(cond_mp,18),
        'sigma_min_mp_estimate': mp.nstr(sigma_min_mp,18),
        'double_LU_forward_error': baselines[0][1],
        'double_LU_backward_error': baselines[0][2],
        'double_CPQR_reported_rank': cp_rank,
        'double_SVD_reported_rank': int(rank_svd),
        'minimum_residual_dps_meeting_target_after_one_correction': min((r[0] for r in relay if r[1] <= TARGET_FWD), default=None),
        'target_forward_error': TARGET_FWD,
    }
    with open(os.path.join(OUT,'SUMMARY.json'),'w') as f: json.dump(summary,f,indent=2)

if __name__ == '__main__':
    run()
