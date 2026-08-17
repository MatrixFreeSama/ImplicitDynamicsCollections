# Bare-C NRII oscillating-field / moving-area flux benchmark

## Status

Literature-grounded **combined stress benchmark**. It is not a bitwise reproduction of any single paper.

The implementation is ISO C17 and uses no external numerical library. The hot path contains no global matrix assembly and no global linear solve.

Recorded build:

```bash
gcc -O3 -march=native -std=c17 nrii_oscillating_field_moving_area_flux.c -lm -o bench
./bench .
```

A separate AddressSanitizer + UndefinedBehaviorSanitizer run completed with exit code 0 and no reported sanitizer error.

## Literature basis

Two published mechanisms are deliberately combined.

1. V. Marusov, **“Measurement of a time-periodic magnetic field by rotating coil,”** arXiv:1211.3974. The paper treats rotating-coil measurement when the magnetic field is periodic in time and derives flux spectra containing combinations of field-cycle and coil-rotation frequencies.

2. S. Liu, H. Liang, B. Xiong, **“An out-of-plane electromagnetic induction based resonant MEMS magnetometer,”** *Sensors and Actuators A: Physical* 285 (2019), DOI `10.1016/j.sna.2018.11.003`. The paper states that the 4S beam is excited in a contractive-extensional mode so that the closed area of the induction coil changes with time, producing an induced electromotive force. The reported resonant scale is around 37.6 kHz.

The present benchmark combines the **time-periodic field / rotating surface** mechanism of the first paper with the **time-varying closed area** mechanism of the second. The frequencies and amplitudes below, other than retaining the 37.6 kHz scale for area oscillation, are synthetic benchmark parameters.

## Model

```text
B(t) = B0 [cos(wB t), sin(wB t), 0]
S(t) = A(t) [cos(wN t), sin(wN t), 0]
A(t) = A0 [1 + epsA cos(wA t)]
Phi(t) = B(t) dot S(t)
       = B0 A0 [1 + epsA cos(wA t)] cos((wB-wN)t)
emf(t) = -dPhi/dt
```

Locked primary parameters:

```text
field-vector frequency      fB = 40.0 kHz
area frequency              fA = 37.6 kHz
surface-normal frequency    fN = 39.1 kHz
B0                              = 1 mT
A0                              = 1e-6 m^2
relative area modulation epsA  = 0.15
dt                              = 10 ns
duration                        = 10 ms
NRII order                      = 16
steps                           = 1,000,000
```

## Frequency mixing

Let `fd = |fB-fN| = 900 Hz`. The product produces three principal flux components:

```text
carrier:         fd       =   900 Hz
lower sideband:  fA-fd    = 36.7 kHz
upper sideband:  fA+fd    = 38.5 kHz
```

Thus rapidly rotating field and surface vectors produce a slow relative carrier, while the breathing area injects two high-frequency sidebands. Faraday voltage weights those components by angular frequency, so small flux sidebands can dominate the voltage waveform.

## NRII discrete contract

Each oscillator pair `[c,s]` satisfies

```text
d/dt [c,s] = w[-s,c].
```

For implicit midpoint:

```text
y1 = y0 + dt A (y0+y1)/2
C1 = dt A y0
Ck = (dt/2) A C{k-1},  k >= 2
```

NRII uses fixed order 16. Residual is evaluated after candidate construction and never creates a correction direction.

The independent solver oracle is the exact 2x2 Cayley-form implicit-midpoint update for each oscillator block. A second continuous oracle uses the closed-form trigonometric field, area, normal, flux, and Faraday-voltage expressions.

## Primary recorded result

The C17 run completed all `1,000,000` steps.

```text
max implicit-midpoint residual inf norm       2.713e-16
max state error vs discrete midpoint oracle   9.730e-11
max state error vs continuous exact           1.323e-3
max flux error vs continuous exact            1.850e-13 Wb
max emf error vs continuous exact             3.951e-8 V
```

During 10 ms:

```text
Bx sign changes              800
area-deviation sign changes  752
flux sign changes             18
emf sign changes             681
```

The field vector and area therefore oscillate hundreds of times while their product/integral observable follows a substantially different zero-crossing pattern.

## Spectral verification

At `dt = 10 ns`, the relative amplitude errors are:

```text
900 Hz carrier flux           9.12e-7
36.7 kHz lower flux sideband  1.10e-5
38.5 kHz upper flux sideband  8.45e-6

900 Hz carrier emf            6.32e-7
36.7 kHz lower emf sideband   1.05e-5
38.5 kHz upper emf sideband   8.94e-6
```

## Area-modulation control

With `epsA=0`, the 36.7 kHz and 38.5 kHz flux sidebands collapse to numerical leakage near `1e-18 Wb`, and the emf has 17 zero crossings over 10 ms.

With `epsA=0.15`, both flux sidebands are approximately `7.5e-11 Wb`, while the emf zero-crossing count rises to 681.

This isolates the consequence of the changing area magnitude from field-vector and surface-normal rotation.

## Time-step sweep

`DT_SWEEP.csv` uses `200 ns, 100 ns, 50 ns, 20 ns, 10 ns`. The implicit-equation residual stays near binary64 roundoff while continuous-solution phase error decreases under refinement. This separates algebraic implicit-solve accuracy from temporal discretization accuracy.

## Interpretation

Supported claim:

> A Matrix-Free fixed-order NRII implicit-midpoint construction can advance a rapidly oscillating electromagnetic-field vector together with a rapidly rotating and periodically changing oriented surface area, while reproducing the analytically predicted carrier/sideband flux structure and Faraday-voltage spectrum without a global matrix solve.

Not supported:

- exact reproduction of either cited device;
- universal superiority over established electromagnetic or structural integrators;
- reduction of arbitrary moving-surface Maxwell problems to this six-state model;
- full Maxwell PDE back-reaction or finite-element structural coupling.

## Files

- `nrii_oscillating_field_moving_area_flux.c`
- `RESULTS.json`
- `MODEL_PARAMETERS.csv`
- `SPECTRUM.csv`
- `AREA_MODULATION_CONTROL.csv`
- `DT_SWEEP.csv`

The larger decimated trajectory is retained in the downloadable benchmark archive rather than committed to this branch.