#!/usr/bin/env python3
"""
Dense constrained-dynamics singular-KKT continuation benchmark.

Purpose
-------
Construct a literature-grounded constrained multibody adaptation in which
the physical joint-space trajectory remains smooth while a dense KKT
formulation loses numerical rank at a manipulator singular configuration.

The benchmark compares:
1. a full-rank Newton/KKT continuation contract that refuses a step when the
   KKT Jacobian loses numerical rank;
2. an SVD pseudoinverse specialist control that can resolve the rank-deficient
   linear subproblem;
3. a Matrix-Free NRII implicit-midpoint continuation in the smooth reduced
   physical coordinates, which never assembles or solves the KKT matrix.

This is not a bitwise reproduction of a published manipulator example.
The literature basis is the documented failure of conventional multibody
constraint formulations near singular configurations and ill-conditioned
DAE Jacobians.

No residual-to-correction loop is used inside NRII.
"""

import csv
import json
import math
import platform
import time
from pathlib import Path

import numpy as np

SEED = 20260817

LINK_LENGTHS = np.array([1.0, 0.8, 0.6], dtype=float)
LINK_MASSES = np.array([1.2, 1.0, 0.8], dtype=float)

NQ = 24
N_AUX = NQ - 3
NCON = 10
KKT_N = NQ + NCON

DT = 0.0025
SINGULAR_STEP = 200
SINGULAR_TIME = DT * SINGULAR_STEP
N_STEPS = 400
TARGET_TIME = DT * N_STEPS
NRII_ORDER = 16

EPS = np.finfo(np.float64).eps

OUT = Path(__file__).resolve().parent
RNG = np.random.default_rng(SEED)


def end_effector_jacobian(q):
    l1, l2, l3 = LINK_LENGTHS
    q1, q2, q3 = q
    t1 = q1
    t2 = q1 + q2
    t3 = q1 + q2 + q3
    s1, c1 = math.sin(t1), math.cos(t1)
    s2, c2 = math.sin(t2), math.cos(t2)
    s3, c3 = math.sin(t3), math.cos(t3)

    J = np.zeros((2, 3), dtype=float)
    J[0, 0] = -(l1*s1 + l2*s2 + l3*s3)
    J[1, 0] =  (l1*c1 + l2*c2 + l3*c3)
    J[0, 1] = -(l2*s2 + l3*s3)
    J[1, 1] =  (l2*c2 + l3*c3)
    J[0, 2] = -l3*s3
    J[1, 2] =  l3*c3
    return J


def manipulator_inertia(q):
    l = LINK_LENGTHS
    m = LINK_MASSES
    I = m*l*l/12.0
    q1, q2, q3 = q
    theta = np.array([q1, q1+q2, q1+q2+q3], dtype=float)

    M = np.zeros((3, 3), dtype=float)
    for i in range(3):
        Jv = np.zeros((2, 3), dtype=float)
        for j in range(i+1):
            dx = 0.0
            dy = 0.0
            for k in range(j, i+1):
                coeff = l[k] * (0.5 if k == i else 1.0)
                dx += -coeff * math.sin(theta[k])
                dy +=  coeff * math.cos(theta[k])
            Jv[:, j] = (dx, dy)
        jw = np.zeros(3, dtype=float)
        jw[:i+1] = 1.0
        M += m[i] * (Jv.T @ Jv) + I[i] * np.outer(jw, jw)
    return M


_Q_raw = RNG.standard_normal((NQ, NQ))
Q, _ = np.linalg.qr(_Q_raw)

_aux_raw = RNG.standard_normal((N_AUX, NCON-2))
_aux_q, _ = np.linalg.qr(_aux_raw)
AUX_CONSTRAINT_ROWS = _aux_q[:, :NCON-2].T

AUX_MASSES = np.logspace(-2.0, 2.0, N_AUX)


def lifted_mass_and_constraints(qcore):
    M0 = np.zeros((NQ, NQ), dtype=float)
    M0[:3, :3] = manipulator_inertia(qcore)
    M0[3:, 3:] = np.diag(AUX_MASSES)

    G0 = np.zeros((NCON, NQ), dtype=float)
    G0[:2, :3] = end_effector_jacobian(qcore)
    G0[2:, 3:] = AUX_CONSTRAINT_ROWS

    M = Q.T @ M0 @ Q
    G = G0 @ Q
    KKT = np.block([
        [M, G.T],
        [G, np.zeros((NCON, NCON), dtype=float)]
    ])
    return M, G, KKT


FREQ = np.zeros(NQ, dtype=float)
FREQ[0] = 0.0
FREQ[1] = math.pi
FREQ[2] = math.pi
FREQ[3:] = np.linspace(0.7, 2.5, N_AUX)

CENTER = np.zeros(NQ, dtype=float)
CENTER[0] = 0.35


def apply_physical_operator(y):
    u = y[:NQ]
    v = y[NQ:]
    out = np.empty_like(y)
    out[:NQ] = v
    out[NQ:] = -(FREQ*FREQ) * u
    return out


def nrii_midpoint_step(y0):
    c = DT * apply_physical_operator(y0)
    y = y0 + c
    tail = float(np.linalg.norm(c, np.inf))

    for _ in range(2, NRII_ORDER + 1):
        c = 0.5 * DT * apply_physical_operator(c)
        y = y + c
        tail = float(np.linalg.norm(c, np.inf))

    residual = y - y0 - 0.5*DT * (
        apply_physical_operator(y0) + apply_physical_operator(y)
    )
    return y, float(np.linalg.norm(residual, np.inf)), tail


def exact_midpoint_step(y0):
    u0 = y0[:NQ]
    v0 = y0[NQ:]
    u1 = np.empty(NQ, dtype=float)
    v1 = np.empty(NQ, dtype=float)

    for i, w in enumerate(FREQ):
        if w == 0.0:
            u1[i] = u0[i] + DT*v0[i]
            v1[i] = v0[i]
            continue

        a = 0.5*DT
        det = 1.0 + (a*w)*(a*w)
        rhs0 = u0[i] + a*v0[i]
        rhs1 = v0[i] - a*w*w*u0[i]
        u1[i] = (rhs0 + a*rhs1) / det
        v1[i] = (-a*w*w*rhs0 + rhs1) / det

    return np.concatenate([u1, v1])


def build_initial_state():
    u = np.zeros(NQ, dtype=float)
    v = np.zeros(NQ, dtype=float)

    amp = 0.25
    w = FREQ[1]
    phi = 2.0 * math.atan(w*DT/2.0)

    u[1] = -amp * math.sin(SINGULAR_STEP * phi)
    v[1] =  amp * w * math.cos(SINGULAR_STEP * phi)

    u[2] = 0.6 * u[1]
    v[2] = 0.6 * v[1]

    aux_rng = np.random.default_rng(SEED + 1)
    u[3:] = aux_rng.normal(scale=0.02, size=N_AUX)
    v[3:] = aux_rng.normal(scale=0.02, size=N_AUX)

    return np.concatenate([u, v])


def kkt_diagnostics(qcore, a_z, relative_rhs_perturbation=0.0):
    M, G, KKT = lifted_mass_and_constraints(qcore)
    a_x_true = Q.T @ a_z
    rhs = np.concatenate([M @ a_x_true, G @ a_x_true])

    if relative_rhs_perturbation:
        rhs = rhs.copy()
        rhs[-1] += relative_rhs_perturbation * max(np.linalg.norm(rhs), 1.0)

    s = np.linalg.svd(KKT, compute_uv=False)
    sigma_max = float(s[0])
    sigma_min = float(s[-1])
    rank_tol = EPS * max(KKT.shape) * sigma_max
    rank = int(np.sum(s > rank_tol))
    cond2 = float(sigma_max / sigma_min) if sigma_min > 0.0 else math.inf

    direct_ok = True
    direct_error = math.nan
    direct_backward = math.nan
    direct_lambda_norm = math.nan
    try:
        sol = np.linalg.solve(KKT, rhs)
        direct_error = float(
            np.linalg.norm(sol[:NQ] - a_x_true) /
            max(np.linalg.norm(a_x_true), 1e-30)
        )
        direct_lambda_norm = float(np.linalg.norm(sol[NQ:]))
        denom = (
            np.linalg.norm(KKT, 2)*np.linalg.norm(sol, 2)
            + np.linalg.norm(rhs, 2)
        )
        direct_backward = float(
            np.linalg.norm(KKT @ sol - rhs, 2) / max(denom, 1e-300)
        )
    except np.linalg.LinAlgError:
        direct_ok = False

    pinv = np.linalg.pinv(KKT, rcond=1e-13)
    psol = pinv @ rhs
    pinv_error = float(
        np.linalg.norm(psol[:NQ] - a_x_true) /
        max(np.linalg.norm(a_x_true), 1e-30)
    )
    pinv_lambda_norm = float(np.linalg.norm(psol[NQ:]))
    pinv_residual_rel = float(
        np.linalg.norm(KKT @ psol - rhs, 2) /
        max(np.linalg.norm(rhs, 2), 1e-30)
    )

    return {
        "sigma_max": sigma_max,
        "sigma_min": sigma_min,
        "cond2": cond2,
        "rank": rank,
        "rank_tol": float(rank_tol),
        "full_rank": bool(rank == KKT_N),
        "direct_solve_returned": bool(direct_ok),
        "direct_acceleration_forward_error": direct_error,
        "direct_backward_error": direct_backward,
        "direct_lambda_norm": direct_lambda_norm,
        "pinv_acceleration_forward_error": pinv_error,
        "pinv_lambda_norm": pinv_lambda_norm,
        "pinv_residual_rel": pinv_residual_rel,
    }


def singularity_sweep():
    rows = []
    for q2 in [1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-10, 1e-12, 0.0]:
        qcore = np.array([0.35, q2, 0.6*q2], dtype=float)
        a_z = np.zeros(NQ, dtype=float)
        a_z[:3] = np.array([0.0, -1.0, -0.6])
        a_z[3:] = np.linspace(-0.3, 0.4, N_AUX)

        d0 = kkt_diagnostics(qcore, a_z, relative_rhs_perturbation=0.0)
        dn = kkt_diagnostics(qcore, a_z, relative_rhs_perturbation=1e-14)

        J = end_effector_jacobian(qcore)
        sj = np.linalg.svd(J, compute_uv=False)

        rows.append({
            "q2_rad": q2,
            "q3_rad": 0.6*q2,
            "task_jacobian_sigma_min": float(sj[-1]),
            "task_jacobian_rank": int(np.linalg.matrix_rank(J)),
            "kkt_sigma_min": d0["sigma_min"],
            "kkt_cond2": d0["cond2"],
            "kkt_numerical_rank": d0["rank"],
            "kkt_full_rank": d0["full_rank"],
            "direct_accel_forward_error_exact_rhs": d0["direct_acceleration_forward_error"],
            "direct_lambda_norm_exact_rhs": d0["direct_lambda_norm"],
            "direct_backward_error_exact_rhs": d0["direct_backward_error"],
            "direct_accel_forward_error_rhs_perturb_1e-14": dn["direct_acceleration_forward_error"],
            "direct_lambda_norm_rhs_perturb_1e-14": dn["direct_lambda_norm"],
            "pinv_accel_forward_error_rhs_perturb_1e-14": dn["pinv_acceleration_forward_error"],
        })
    return rows


def run_continuation():
    y_nrii = build_initial_state()
    y_oracle = y_nrii.copy()

    trace = []
    max_nrii_residual = 0.0
    max_nrii_vs_oracle = 0.0

    full_rank_baseline_alive = True
    baseline_committed_steps = 0
    baseline_stop = None

    singular_snapshot = None

    for step in range(1, N_STEPS + 1):
        y_nrii, nrii_resid, tail = nrii_midpoint_step(y_nrii)
        y_oracle = exact_midpoint_step(y_oracle)

        state_err = float(np.max(np.abs(y_nrii - y_oracle)))
        max_nrii_vs_oracle = max(max_nrii_vs_oracle, state_err)
        max_nrii_residual = max(max_nrii_residual, nrii_resid)

        z = y_nrii[:NQ] + CENTER
        v = y_nrii[NQ:]
        qcore = z[:3].copy()

        a_z = -(FREQ*FREQ) * (z - CENTER)
        diag = kkt_diagnostics(qcore, a_z, relative_rhs_perturbation=0.0)

        if full_rank_baseline_alive:
            if not diag["full_rank"]:
                full_rank_baseline_alive = False
                baseline_stop = {
                    "attempted_step": step,
                    "attempted_time_s": step*DT,
                    "last_committed_step": step-1,
                    "last_committed_time_s": (step-1)*DT,
                    "reason": "rank_deficient_KKT_Jacobian",
                    "kkt_rank": diag["rank"],
                    "kkt_dimension": KKT_N,
                    "kkt_sigma_min": diag["sigma_min"],
                    "kkt_cond2": diag["cond2"],
                }
            else:
                baseline_committed_steps = step

        if step == SINGULAR_STEP:
            singular_snapshot = {
                "step": step,
                "time_s": step*DT,
                "q_core": qcore.tolist(),
                "q2_abs": float(abs(qcore[1])),
                "q3_abs": float(abs(qcore[2])),
                "task_jacobian_singular_values": np.linalg.svd(
                    end_effector_jacobian(qcore), compute_uv=False
                ).tolist(),
                "kkt": diag,
            }

        trace.append({
            "step": step,
            "time_s": step*DT,
            "q1_rad": float(qcore[0]),
            "q2_rad": float(qcore[1]),
            "q3_rad": float(qcore[2]),
            "v2_rad_s": float(v[1]),
            "nrii_residual_inf": nrii_resid,
            "nrii_tail_inf": tail,
            "nrii_vs_midpoint_oracle_max_abs": state_err,
            "kkt_sigma_min": diag["sigma_min"],
            "kkt_cond2": diag["cond2"],
            "kkt_rank": diag["rank"],
            "kkt_full_rank": diag["full_rank"],
            "full_rank_newton_kkt_baseline_alive_after_step": full_rank_baseline_alive,
            "pinv_accel_forward_error": diag["pinv_acceleration_forward_error"],
        })

    return {
        "trace": trace,
        "summary": {
            "target_time_s": TARGET_TIME,
            "target_steps": N_STEPS,
            "singular_time_s": SINGULAR_TIME,
            "singular_step": SINGULAR_STEP,
            "full_rank_newton_kkt_baseline_reached_target": bool(
                baseline_committed_steps == N_STEPS
            ),
            "full_rank_newton_kkt_baseline_committed_steps": baseline_committed_steps,
            "full_rank_newton_kkt_baseline_stop": baseline_stop,
            "nrii_reached_target": True,
            "nrii_committed_steps": N_STEPS,
            "nrii_max_residual_inf": max_nrii_residual,
            "nrii_vs_midpoint_oracle_max_abs": max_nrii_vs_oracle,
            "singular_snapshot": singular_snapshot,
        },
    }


def write_csv(path, rows):
    if not rows:
        return
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)


def main():
    t0 = time.perf_counter()

    sweep = singularity_sweep()
    cont = run_continuation()

    parameters = [
        ("seed", SEED),
        ("physical_joint_coordinates", 3),
        ("lifted_generalized_coordinates", NQ),
        ("auxiliary_coordinates", N_AUX),
        ("constraints", NCON),
        ("dense_KKT_dimension", KKT_N),
        ("dt_s", DT),
        ("singular_step", SINGULAR_STEP),
        ("singular_time_s", SINGULAR_TIME),
        ("target_steps", N_STEPS),
        ("target_time_s", TARGET_TIME),
        ("NRII_order", NRII_ORDER),
        ("link_lengths", ";".join(map(str, LINK_LENGTHS))),
        ("link_masses", ";".join(map(str, LINK_MASSES))),
        ("aux_mass_min", float(AUX_MASSES.min())),
        ("aux_mass_max", float(AUX_MASSES.max())),
    ]

    write_csv(OUT / "SINGULARITY_SWEEP.csv", sweep)
    write_csv(OUT / "CONTINUATION_TRACE.csv", cont["trace"])

    with (OUT / "MODEL_PARAMETERS.csv").open("w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["parameter", "value"])
        w.writerows(parameters)

    results = {
        "test_name": "Dense Constrained-Dynamics Singular-KKT Continuation Benchmark",
        "status": "literature-grounded benchmark adaptation; not a bitwise reproduction of a published implementation",
        "runtime": {
            "python": platform.python_version(),
            "platform": platform.platform(),
            "numpy": np.__version__,
            "elapsed_s": time.perf_counter() - t0,
        },
        "literature_basis": [
            {
                "citation": "Ider & Amirouche, Computers & Structures 33(1), 1989, pp.129-137",
                "doi": "10.1016/0045-7949(89)90135-1",
                "supported_point": "constraint equations can become linearly dependent at instantaneous singular configurations, producing singularities in conventional multibody simulation algorithms",
            },
            {
                "citation": "Parida & Raha, Applied Mathematics and Computation 215(3), 2009, pp.1224-1243",
                "doi": "10.1016/j.amc.2009.06.063",
                "supported_point": "singular configurations and high-index constraints can produce ill-conditioned discretization Jacobians; regularization increases small singular values",
            },
        ],
        "contracts": {
            "matrix_dependent_baseline": "full-rank Newton/KKT correction contract; stop when the KKT Jacobian loses numerical rank",
            "specialist_control": "SVD pseudoinverse on the same rank-deficient KKT subproblem",
            "nrii": "Matrix-Free fixed-order NRII implicit-midpoint coefficient construction in smooth physical coordinates; residual audit only",
        },
        "continuation": cont["summary"],
        "singularity_sweep": sweep,
        "guardrails": [
            "The benchmark is a literature-grounded adaptation, not the exact three-link numerical setup from the 1989 paper.",
            "The exact rank-deficient KKT system is consistent; the physical acceleration remains defined while the Lagrange multiplier is non-unique.",
            "The SVD pseudoinverse control is included because singular KKT does not imply that every established method must fail.",
            "The NRII path bypasses this KKT formulation; it does not prove that NRII resolves every physically underdetermined singularity.",
            "No residual-driven state correction is used inside NRII.",
        ],
    }

    (OUT / "RESULTS.json").write_text(
        json.dumps(results, indent=2), encoding="utf-8"
    )

    print(json.dumps({
        "test_name": results["test_name"],
        "baseline_reached_target": cont["summary"]["full_rank_newton_kkt_baseline_reached_target"],
        "baseline_stop": cont["summary"]["full_rank_newton_kkt_baseline_stop"],
        "nrii_reached_target": cont["summary"]["nrii_reached_target"],
        "nrii_max_residual_inf": cont["summary"]["nrii_max_residual_inf"],
        "nrii_vs_midpoint_oracle_max_abs": cont["summary"]["nrii_vs_midpoint_oracle_max_abs"],
        "singular_snapshot": cont["summary"]["singular_snapshot"],
        "elapsed_s": results["runtime"]["elapsed_s"],
    }, indent=2))


if __name__ == "__main__":
    main()
