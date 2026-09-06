# NRII Semiconductor Constitutive-Model Campaign

**Date:** 2026-09-06  
**Implementation:** bare C, FP64, single-threaded CPU  
**Scope:** controlled normalized 3D constitutive/transport benchmarks, not production TCAD device validation.

## 1. Campaign design

The main campaign contains **1800 validated points**:

- 5 cubic grid sizes: 10^3, 14^3, 18^3, 22^3, 26^3;
- 12 constitutive families;
- 5 constitutive-strength levels;
- 6 implicit-step levels per strength.

The tested families are:

1. Shockley-Read-Hall-type rational recombination;
2. Auger recombination;
3. radiative recombination;
4. one explicit trap state;
5. eight explicit trap states;
6. sixteen explicit trap states;
7. field-dependent mobility;
8. electrothermal mobility/heating coupling;
9. impact-ionization source;
10. heterojunction-localized interface recombination and mobility contrast;
11. combined PIN-like electrothermal + eight-trap + impact-ionization model;
12. a deliberately nonlocal global-mean negative control.

All main-campaign points satisfy

    relative NRII residual <= 1e-8
    relative NRII/Newton difference <= 1e-7

with the observed maxima substantially below these gates.

The comparison is against a matrix-free Newton-GMRES implementation using the same implicit relation and the same constitutive functions. These are controlled solver-structure tests. They do not include production Scharfetter-Gummel discretization, AMG/ILU, commercial TCAD continuation strategies, exact semiconductor material calibration, or a self-consistent elliptic Poisson solve.

## 2. Main measured ranking

| Constitutive family | Points | Median speedup | IQR | Median q_eff | Median NRII order |
|---|---:|---:|---:|---:|---:|
| field_mobility | 150 | 2.83x | 2.37-3.09x | 0.0130 | 6.0 |
| radiative | 150 | 2.73x | 2.31-2.98x | 0.0136 | 6.0 |
| heterojunction_interface | 150 | 2.64x | 2.43-2.91x | 0.0342 | 7.0 |
| impact_ionization | 150 | 2.49x | 2.30-2.83x | 0.0161 | 7.0 |
| auger | 150 | 2.14x | 1.75-2.35x | 0.0137 | 6.0 |
| trap_1 | 150 | 1.85x | 1.25-2.29x | 0.0099 | 6.0 |
| global_mean_negative_control | 150 | 1.39x | 1.25-1.49x | 0.0205 | 7.0 |
| trap_8 | 150 | 1.36x | 0.99-1.64x | 0.0081 | 5.0 |
| combined_pin_electrothermal_traps | 150 | 1.33x | 1.13-1.51x | 0.0149 | 6.0 |
| trap_16 | 150 | 1.18x | 0.87-1.46x | 0.0071 | 5.0 |
| srh | 150 | 1.13x | 0.96-1.40x | 0.0246 | 7.0 |
| electrothermal | 150 | 1.13x | 0.91-1.41x | 0.0323 | 8.0 |

The most stable favorable families in this campaign are field-dependent mobility, radiative recombination, heterojunction-localized coupling, impact ionization below its coefficient-radius boundary, and Auger recombination. Their aggregate median speed ratios are approximately:

- field-dependent mobility: **2.83x**;
- radiative recombination: **2.73x**;
- heterojunction interface: **2.64x**;
- impact ionization: **2.49x**;
- Auger recombination: **2.14x**.

These values are medians over the complete five-size parameter campaign, not selected maxima.

## 3. New structural variable: coefficient-algebra density

Locality and low q are not sufficient to predict NRII speed.

The explicit-trap series provides a controlled example. The median speedup over all five grid sizes changes from approximately

    1 trap  : 1.85x
    8 traps : 1.36x
    16 traps: 1.18x

even though q_eff becomes smaller rather than larger as the number of local trap variables increases.

The reason is that every retained NRII order must propagate the local trap coefficient algebra. The appropriate NRII cost model therefore requires a per-order constitutive cost:

    T_NRII ~ p * N_active * C_coeff

where C_coeff includes convolution, reciprocal-series, trap-state, and local multiphysics work.

The favorable-region criterion should consequently be extended from

    local + low q + expensive global correction

to

    local + low q + controlled Gamma + low/moderate coefficient-algebra density
    + expensive competing global correction.

## 4. Impact-ionization boundary sweep

A separate 72-point impact-ionization stress sweep extends the source strength and step size well beyond the main campaign.

- valid points: **65/72**;
- maximum q_eff among valid points: **0.608**;
- minimum q_eff among failed points: **0.693**.

With the current 24-order cap, valid NRII operation extends to approximately q_eff = 0.61 in the sampled directions. Accuracy failures first appear around q_eff = 0.69, and deliberately pushed cases reach q_eff > 1.

Apparent large speed ratios on failed points are rejected. For example, the q_eff > 1 cases can terminate quickly while producing an invalid solution; these are not speed wins.

The experiment supports the expected sequence

    impact strength / step size increases
        -> q_eff rises
        -> retained order rises
        -> the NRII advantage contracts
        -> the power-series path reaches its practical order/accuracy boundary.

## 5. Combined electrothermal + traps + impact model

The combined model is deliberately coefficient-heavy. Its main-campaign median speedup is **1.33x**, but in the stronger 72-point boundary sweep the median speedup over valid points falls to **0.67x**.

This is a counterexample to the statement that a local multiphysics block automatically favors NRII. The model remains local, and q_eff can remain below one, but each NRII coefficient requires:

- rational SRH coefficient propagation;
- cubic Auger products;
- radiative products;
- eight trap-state recurrences;
- temperature-dependent reciprocal mobility;
- local transport coefficients;
- Joule-like heating;
- impact-ionization generation.

At sufficiently high constitutive density, the coefficient-construction layer itself becomes the dominant cost.

## 6. Heterojunction and localized constitutive structure

The heterojunction-interface family remains near **2.64x** median speedup across the complete campaign. The strong constitutive variation is confined to a thin interface while the field operator remains finite-neighborhood.

This is structurally favorable for future spatially adaptive NRII because a thin interface can occupy a substantially smaller subset of a 3D volume than the full device. The present campaign measures full NRII only; it does not claim an active-sparse timing gain for this semiconductor case yet.

## 7. Electrothermal result

The electrothermal family is much less favorable in the current full-NRII implementation: median speedup **1.13x**. Its median q_eff is still only **0.032**, demonstrating again that low q alone is not a performance guarantee.

Temperature-dependent reciprocal mobility and local heat-generation products materially increase coefficient-algebra cost. Spatial activation or a relay that leaves the temperature field on another solver may therefore be more attractive than a monolithic all-field NRII construction.

## 8. Physical globality remains a separate issue

The global-mean negative control retains a genuine nonlocal operation. NRII can represent the coefficient sequence, but every coefficient order must evaluate the nonlocal mean. This is intentionally different from the local constitutive families.

A production drift-diffusion semiconductor model also contains the elliptic Poisson equation for electrostatic potential. The present campaign does not remove or approximate that physical globality. The intended future architecture is instead:

    global electrostatic field solve
    + local/finite-neighborhood carrier and constitutive coefficient propagation
    + solver relay when q, Gamma, or coefficient-algebra cost leaves the favorable region.

## 9. Revised semiconductor golden-region hypothesis

The data support five independent coordinates:

1. **Coefficient contraction:** q_eff controls retained order and the series-radius boundary.
2. **Transient amplification:** Gamma controls finite-precision safety.
3. **Spatial locality:** the constitutive/transport operator must remain finite-neighborhood if solver-induced globality is to be removed.
4. **Coefficient-algebra density:** the cost of one NRII order can vary by more than an order of conceptual complexity between simple mobility/radiative models and many-trap electrothermal models.
5. **Competing global-solve cost:** NRII has little to remove when Newton's correction problem is already unusually cheap.

The measured high-value semiconductor candidates are therefore not simply the most nonlinear models. They are models in which strong local constitutive structure is expensive for a global Newton correction but inexpensive to propagate coefficient-wise.

## 10. Current candidate ranking for the next stage

**Highest priority**
- heterojunction-localized carrier transport with spatial activation;
- subcritical impact-ionization front with moving activation and relay at rising q;
- field-dependent mobility in large 3D devices;
- local radiative/Auger carrier dynamics.

**Conditional**
- one/few dynamic trap states;
- moderate electrothermal coupling;
- combined electrothermal/trap/impact systems only with spatial activation or solver relay.

**Poor full-NRII target in the present implementation**
- large explicit trap-bin counts applied throughout the whole 3D volume;
- very dense local constitutive stacks in which coefficient algebra dominates;
- q_eff approaching the power-series boundary;
- physically nonlocal electrostatic formulations treated as if they were local.

## 11. Evidence boundary

This campaign is a controlled normalized constitutive benchmark. It establishes solver-mechanism trends, not semiconductor-device accuracy. Before any engineering claim, the same tests should be repeated with calibrated silicon/GaN material laws, Scharfetter-Gummel fluxes, self-consistent Poisson coupling, production preconditioned Newton baselines, and device-level verification quantities such as current-voltage curves, charge conservation, breakdown voltage, and thermal balance.
