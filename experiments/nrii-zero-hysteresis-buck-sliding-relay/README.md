# NRII Zero-Hysteresis Buck Sliding Relay

## Status

Research benchmark only. This is not a production-ready converter simulator and not a universal performance claim for NRII. The purpose is to validate solver-domain separation on an electrical switching singularity.

## Model

Ideal buck converter:

\[
L\dot i = dV_{in}-v,\qquad
C\dot v = i-\frac{v}{R}.
\]

Parameters:

- \(V_{in}=24\,\mathrm{V}\)
- \(L=100\,\mu\mathrm{H}\)
- \(C=100\,\mu\mathrm{F}\)
- \(R=4\,\Omega\)
- \(I_{ref}=3\,\mathrm{A}\)

The ideal zero-hysteresis relay is

\[
d=1\quad\text{for }i<I_{ref},\qquad
d=0\quad\text{for }i>I_{ref}.
\]

At the switching surface \(s=i-I_{ref}=0\), if the ON and OFF vector fields point toward the surface, the binary model has a sliding/chattering singularity. The generalized relay law is written as

\[
0\in s+N_{[0,1]}(d),
\]

or equivalently

\[
d-\operatorname{Proj}_{[0,1]}(d-\gamma s)=0.
\]

This allows \(d\in[0,1]\) on the sliding surface.

## Discrete contract

The smooth ON and OFF subsystems use implicit midpoint. For each fixed switch state, NRII parameterizes the discrete equation by \(h\in[0,1]\):

\[
Y(h)=y_n+h\Delta t\left[A\frac{y_n+Y(h)}2+b_d\right].
\]

With \(Y(h)=\sum C_kh^k\):

\[
C_0=y_n,\qquad
C_1=\Delta t(Ay_n+b_d),\qquad
C_k=\frac{\Delta t}{2}AC_{k-1},\;k\ge2.
\]

The residual is evaluated only after candidate construction. It never corrects the NRII state.

The smooth candidate must also satisfy its switching-domain contract:

- ON candidate: \(i_{n+1}<I_{ref}\)
- OFF candidate: \(i_{n+1}>I_{ref}\)

An equation-valid but domain-invalid candidate is rejected.

## Primary case

Initial state:

\[
i_n=3\,\mathrm{A},\qquad v_n=10\,\mathrm{V},\qquad \Delta t=10\,\mu\mathrm{s}.
\]

At the switching surface:

\[
\dot i_{ON}=140000\,\mathrm{A/s},\qquad
\dot i_{OFF}=-100000\,\mathrm{A/s}.
\]

Both fields point toward the switching surface, so this is an interior sliding case.

### Smooth NRII branches

NRII order: 40.

OFF candidate:

\[
i_{n+1}=2.000000000000001\,\mathrm A,\qquad v_{n+1}=10.0\,\mathrm V.
\]

Midpoint residual infinity norm: \(8.88\times10^{-16}\). OFF requires \(i_{n+1}>3\,\mathrm A\), so the candidate is domain-invalid.

ON candidate:

\[
i_{n+1}=4.394088669950738\,\mathrm A,\qquad v_{n+1}=10.118226600985222\,\mathrm V.
\]

Midpoint residual infinity norm: \(7.77\times10^{-16}\). ON requires \(i_{n+1}<3\,\mathrm A\), so this candidate is also domain-invalid.

An independent dense 2x2 midpoint solve agrees with the NRII branch candidates to a maximum absolute error of \(8.88\times10^{-16}\).

Because each fixed-switch branch is affine and has one unique root, and both unique roots violate their own switching inequalities, the benchmark emits a `NoBinaryRootCertificate` for this one-step binary problem.

## Strict ordinary Newton control

The control applies ordinary Newton directly to the exact discontinuous binary residual. The tie at \(i=I_{ref}\) is assigned to ON only to select the first branch.

```text
k=0: i=3.0000000000, v=10.0000000000, branch=ON
k=1: i=4.3940886700, v=10.1182266010, branch=OFF
k=2: i=2.0000000000, v=10.0000000000, branch=ON
k=3: i=4.3940886700, v=10.1182266010, branch=OFF
...
```

After the first jump, the residual infinity norm remains approximately 24 rather than approaching zero. The iteration enters an exact period-2 active-set cycle.

This does **not** mean Newton methods are universally unable to solve switching systems. A smoothed relay, a complementarity formulation, a semismooth method, or another generalized nonsmooth treatment changes the mathematical contract and can resolve the event.

## Projected semismooth relay

When both smooth NRII candidates fail only the switching-domain gate, the relay passes a switching event to a projected semismooth solver for

\[
0\in s+N_{[0,1]}(d).
\]

Starting from the natural event state and \(d_0=1\), the projected semismooth solve takes one Newton update and returns

\[
i_{n+1}=3.0\,\mathrm A,
\]

\[
v_{n+1}=10.049382716049383\,\mathrm V,
\]

\[
d=0.41769547325102885.
\]

Generalized residual infinity norm: \(2.22\times10^{-15}\). An independent closed-form sliding solution agrees to \(5.55\times10^{-17}\) maximum absolute error.

## Time-step sweep

The benchmark sweeps `1e-7, 3e-7, 1e-6, 3e-6, 1e-5, 3e-5, 5e-5 s` for the same switching-surface initial state. In every tested case:

- both fixed binary NRII roots are equation-valid but switching-domain-invalid;
- strict ordinary Newton on the exact discontinuous binary residual fails to converge and enters a period-2 cycle;
- the projected semismooth relay converges;
- the relay mode is `SLIDING`.

The sweep is a validation set for this deliberately locked model, not an asymptotic complexity or performance benchmark.

## Hysteresis-to-zero limit

For a symmetric current hysteresis half-band \(\delta\), using the local slopes at \(v\approx10\,\mathrm V\),

\[
T\approx\frac{2\delta}{\dot i_{ON}}+\frac{2\delta}{|\dot i_{OFF}|},\qquad f\approx\frac1T.
\]

The generated sweep shows \(f\propto1/\delta\), so the switching frequency diverges as the idealized hysteresis width approaches zero. The generalized sliding solution replaces that unresolved infinite switching sequence with an interior duty ratio.

## Solver-relay interpretation

```text
smooth ON/OFF implicit midpoint
            |
            v
     Matrix-Free NRII
            |
            +--> equation residual audit
            |
            +--> branch-domain audit
                     |
        both branches invalid at switch
                     |
                     v
             SwitchingEventBaton
                     |
                     v
        projected semismooth solver
                     |
                     v
       SLIDING / ON / OFF generalized state
```

The intended claim is narrow:

> A solver can be valid on each smooth electrical subsystem and still be invalid at the switching singularity. Solver incompatibility at the local mathematical contract does not imply that the global relay must fail.

## Files

- `nrii_zero_hysteresis_buck_sliding_relay.py`
- `RELAY_RESULTS.json`
- `DT_SWEEP.csv`
- `VOLTAGE_GATE.csv`
- `HYSTERESIS_LIMIT.csv`

Run:

```bash
python nrii_zero_hysteresis_buck_sliding_relay.py
```

Dependencies: Python 3 and NumPy.
