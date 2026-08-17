# NRII Solver Relay Graph: Transcendental Fold Experiment

This branch tests a general solver-handoff architecture rather than requiring one solver to finish every stage alone.

## Test equation

We use the scalar transcendental family

F(x; mu) = x + mu + exp(-x) = 0.

The solution curve is mu(x) = -x - exp(-x). At the fold (x*,mu*)=(0,-1):

F=0, F_x=1-exp(-x)=0, F_xx=exp(-x)=1.

For mu<-1 there are two real roots; at mu=-1 they merge into a double root; for mu>-1 no real root exists. The case combines a transcendental equation, multiple roots, branch ambiguity, a singular Jacobian, and a parameter turning point.

## Relay architecture

For current anchor x_a and target mu, NRII uses the re-anchored fixed-point homotopy

X(s)=x_a+s[-mu-exp(-X(s))-x_a].

With X=sum C_k s^k and E=exp(-X)=sum E_k s^k, the forward coefficient closure is

C0=x_a, E0=exp(-C0),
C1=-mu-E0-C0,
E_n=-(1/n) sum_{k=1}^n k C_k E_{n-k},
C_{n+1}=-E_n.

No residual-to-Jacobian-to-correction loop is inserted into NRII. Fixed declared orders 8, 12, 16 are evaluated with Pade continuation. Their spread constructs candidate interval handoff data. Residual is not used to increase NRII order.

If the NRII interval passes a scalar numerical Krawczyk inclusion gate, its center is handed to damped Newton for terminal convergence. The corrected root is fed back as the next NRII anchor.

When a fixed-parameter handoff is unavailable near the fold, the scheduler switches to pseudo-arclength continuation. Once the turning region is crossed, control returns to NRII -> Krawczyk -> Newton.

Tested relay:

NRII <-> Krawczyk gate <-> Newton <-> Pseudo-arclength

Important limitation: Krawczyk uses ordinary binary64 std::exp without directed-rounding interval libraries, so it is a numerical inclusion gate, not a formal machine-verifiable certificate.

## Main result

Start: x=-1, mu=-1.718281828459045.

Final positive branch: x=1.0276300401633864, mu=-1.3854840962851955, |F|=6.66e-15.

30 accepted steps:
- 28 NRII -> Krawczyk -> Newton steps;
- 2 pseudo-arclength fold-crossing steps;
- 25 total terminal Newton updates;
- 4 augmented Newton corrector updates inside pseudo-arclength.

Closest accepted point to the singular fold: x=-0.006358907900956, mu=-1.000020260777552. The next pseudo-arclength step lands on the positive side and the normal NRII/Krawczyk/Newton route resumes.

## Controls

Direct Newton from the original negative-branch anchor at the final relay parameter converges cleanly but to the wrong negative root x=-0.766294974976151.

Naive Newton parameter continuation replaying the relay's mu schedule also remains on the negative branch and finishes at the same wrong root.

At the exact double root mu=-1, Newton requires 24 updates to reach about 3.55e-15, illustrating loss of ordinary quadratic convergence.

A direct NRII handoff state from the original negative anchor at the final parameter also remains on the negative branch. NRII is not branch-omniscient; branch identity must be preserved by the relay.

A Krawczyk box centered on the exact double root cannot satisfy a unique-root inclusion condition because the derivative vanishes there.

Pure pseudo-arclength is a strong specialist control and crosses the fold by itself. This experiment does not claim that the relay is faster or universally superior to pseudo-arclength. Its purpose is to demonstrate automatic solver-role handoff with a common structured interface.

## Interpretation

Solver A does not have to produce the final answer. It can produce the information Solver B is missing.

The handoff data can carry a state, interval or domain, branch identity, tangent, local model, approximate inverse, search subspace, or verification status. The scheduler chooses the next specialist and can feed accepted results back into NRII.

This experiment does not establish a universal solver, infinite continuation, or a theorem that NRII always finds the correct branch. True branch termination, nonexistence, an unavoidable singularity, or failure of every available handoff still terminates the relay.

## Files

- solver_relay_fold.cpp: complete C++17 reproducer
- RELAY_TRACE.csv: every accepted relay step
- RELAY_SUMMARY.csv: aggregate relay metrics
- METHOD_CONTROLS.csv: compact solver controls
- FOLD_STRESS.csv: near-double-root stress scan

Build: g++ -O3 -std=c++17 solver_relay_fold.cpp -o solver_relay_fold
