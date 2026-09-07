# Physical-Unit Si PIN Device Follow-Up

**Date:** 2026-09-07  
**Implementation:** bare C, FP64, single-threaded CPU  
**Status:** device-level controlled validation study; not a production TCAD solver.

## 1. Purpose

This study moves beyond normalized semiconductor constitutive benchmarks and restores a physical-unit silicon PIN device with self-consistent electrostatics, finite-volume Scharfetter-Gummel transport, physical contacts, recombination, and an avalanche criterion.

The central question is not whether NRII can represent the equations. It is whether full-domain NRII remains computationally useful after the fastest physical relaxation modes of a real heavily doped semiconductor contact are restored.

The answer is negative for the present full-domain formulation, and the failure mechanism is physically interpretable.

## 2. Device and numerical model

The model is a 24 um one-dimensional silicon PIN diode at 300 K with p+ and n+ contacts at approximately 1e18 cm^-3. The implementation uses:

- physical electron charge and silicon permittivity;
- intrinsic carrier concentration and physical electron/hole mobilities;
- self-consistent Poisson electrostatics;
- finite-volume Scharfetter-Gummel carrier flux;
- ohmic contact states;
- SRH recombination;
- an Okuto-Crowell-style electron avalanche integral for breakdown estimation;
- a direct analytic block-tridiagonal Newton solve for the coupled one-dimensional system.

The direct Newton baseline is intentionally strong. In one dimension it has O(N) structured linear algebra and does not incur Krylov or global-reduction overhead.

## 3. DC and breakdown verification

The avalanche integral crossing converges with mesh refinement:

| Nodes | dx (um) | Estimated breakdown voltage (V) | Emax at crossing (V/cm) |
|---:|---:|---:|---:|
| 61 | 0.400 | -492.004 | 2.418e+05 |
| 121 | 0.200 | -487.502 | 2.420e+05 |
| 241 | 0.100 | -485.307 | 2.421e+05 |
| 481 | 0.050 | -484.233 | 2.421e+05 |
| 961 | 0.025 | -483.718 | 2.422e+05 |

A linear-in-dx extrapolation gives approximately **-483.12 V**, and a quadratic fit gives **-483.18 V**. The 961-node result is **-483.72 V**.

This is an ionization-integral breakdown estimate rather than a complete avalanche-multiplication IV solution. It is used here as a device-scale consistency quantity.

## 4. Full-domain NRII discovers the charge-relaxation time

Across N = 61, 121 and 241, the measured NRII coefficient-tail ratio is fitted by

q_eff ~= dt / tau_fit

with

- tau_fit = **4.942426e-15 s**
- R^2 = **0.999999805**

Using the physical contact parameters in the code, the textbook charge-relaxation time

tau_q = epsilon_si / (q * mu_n * N_contact)

is

- tau_q = **4.541682e-15 s**

The fitted NRII scale differs by approximately **8.82%**.

This gives a direct physical interpretation of the NRII radius in this device. A heavily doped contact contains a dielectric/charge-relaxation mode on a few-femtosecond scale. Full-domain coefficient construction sees that mode even when the requested device transient is much slower.

For a linear relaxation mode

d rho / dt = -rho / tau_q,

backward Euler gives a rational implicit update whose forward power-series representation has ratio dt/tau_q. The observed device data follow the same scaling.

## 5. Full-domain performance result

The full-domain NRII path is therefore constrained by approximately

dt < O(tau_q)

for power-series convergence, with rapid order growth as q_eff approaches one.

In the valid low-q region, direct block-Newton is already faster. This is not a weak-baseline artifact: the one-dimensional Newton correction uses a structured O(N) block-Thomas solve.

The result establishes a clear Newton-dominant region:

- stiff, heavily doped contact modes;
- cheap structured direct global solve;
- no expensive Krylov/global-reduction layer available for NRII to remove.

## 6. Contact-mode condensation experiment

A multiscale experiment freezes the rapidly relaxing p+/n+ contact carrier modes and applies coefficient propagation only to the slower drift region.

This substantially enlarges the observed coefficient radius. For example, around dt = 1e-13 s the condensed models can retain q_eff of order 0.3-0.4 rather than the full-domain value far above one.

However, the condensed state is not a complete solution of the original equations. Its residual in the uncondensed system remains around 1e-4 in representative runs. It must therefore be classified as an approximate constructive seed, not a converged NRII solution.

## 7. NRII-to-Newton relay

The condensed NRII state was then passed to the complete direct Newton solver.

The Newton correction restores the complete-system residual to the direct-Newton accuracy level, validating the relay mathematically. The present timing is not favorable: direct Newton from the original state remains faster than constructing the condensed NRII seed and then correcting it.

The relay therefore demonstrates information transfer, not a current speed advantage.

## 8. Revised physical preflight criterion

For semiconductor regions dominated by charge relaxation, the generic spectral preflight quantity

q_pred ~= dt * rho(J)

can be replaced by an inexpensive physical estimate

q_pred ~= dt / tau_min,

with

tau_min ~= min_x epsilon(x) / sigma(x).

This provides a material-parameter gate before any NRII coefficient is generated.

It also identifies a principled multiscale route: fast quasi-equilibrium/contact modes should be eliminated, condensed, or assigned to another solver before NRII is applied to slower depletion, heterojunction, avalanche-front, thermal, or trap dynamics.

## 9. Main conclusion

The physical-unit device study changes the semiconductor NRII interpretation.

The favorable architecture is not full-domain NRII over p+, intrinsic/drift, and n+ regions. A more defensible architecture is:

fast-contact or quasi-equilibrium specialist
+ global Poisson/electrostatic solve
+ NRII on slower local carrier/constitutive regions
+ full-system correction or relay when required.

The result also gives q_eff a direct physical interpretation in this device: it measures the requested implicit step relative to the fastest unresolved relaxation mode.

## 10. Limitations

This is still not a production TCAD validation package. Remaining items include calibrated device geometry/material decks, full avalanche multiplication feedback rather than an ionization-integral criterion, Fermi-Dirac statistics where required, high-field mobility models, contact resistance, calibrated lifetime models, production Gummel/continuation strategies, and multidimensional verification.

The strongest current result is the physically interpretable negative boundary, not a device-level speedup claim.
