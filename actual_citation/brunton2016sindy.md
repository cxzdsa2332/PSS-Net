# brunton2016sindy

- **Title:** Discovering governing equations from data by sparse identification of nonlinear dynamical systems
- **Authors:** Steven L. Brunton, Joshua L. Proctor, J. Nathan Kutz
- **Journal/Year:** PNAS, 113(15), 3932–3937, 2016
- **DOI:** 10.1073/pnas.1517384113
- **Source consulted:** Europe PMC REST and Semantic Scholar (DOI:10.1073/pnas.1517384113)

## Verbatim opening (abstract)

"Extracting governing equations from data is a central challenge in many diverse areas
of science and engineering."

## Remainder [summary — condensed by retrieval tool; key verbatim fragments quoted]

The method combines sparsity-promoting techniques with machine learning to identify
nonlinear governing equations, assuming "only a few important terms that govern the
dynamics, so that the equations are sparse in the space of possible functions." Sparse
regression then determines the fewest terms needed. It is demonstrated on oscillators,
the Lorenz system, and fluid dynamics, and generalizes to parameterized, time-varying,
and externally forced systems.

## Supports in intro

B3 — the canonical data-driven discovery of nonlinear dynamics by sparse regression on
a library of candidate functions. Important contrast for PSS-Net: SINDy regresses the
**time derivative** dx/dt against the library, so it requires estimating derivatives
from sampled trajectories (noise-amplifying); PSS-Net instead uses the steady-state
(dx/dt = 0) algebraic equations and avoids numerical differentiation.
