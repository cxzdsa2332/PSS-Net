# Figure 2 caption

Source figure object: `sim_script/02_scaling_design/Fig2.R` (`Fig2`). The caption
is kept here rather than inside the rendered figure; panel subtitles were removed
from the individual panels for the assembled multi-panel layout.

---

**Figure 2 | Sample complexity and active perturbation design for PSS-Net**
(*p* = 8, *s* = 2 unless noted).

**(a)** Node-wise ADSIHT edge recovery (MCC) versus the rescaled budget
*N*/(*s* log *p*) for *p* = 8, 50, 100; dashed verticals mark
*N* = 1, 2, 5 × *s* log *p* and the interpolated budget reaching MCC = 0.8.

**(b)** Conceptual input-space schematic of four designs continued from a shared
random pilot (black rings): random, maximin, oracle and pilot-estimated
D-optimal. Oracle and pilot D-optimal differ only in how the response map is
obtained (true map vs noisy pilot estimate).

**(c)** Regime-by-budget map of the paired mean MCC gain over random at matched
budget; each cell also lists the percentage of seeds with a positive gain and a
dot marks the stronger structured design. D-optimal here is the oracle upper
bound (true map), discounted to a pilot estimate in (e).

**(d)** Budget *N*/(*s* log *p*) required to reach target MCC (0.5, 0.6);
leftward arrows from random indicate sample savings and ">max" marks targets
unmet within the swept budget. Exact D-optimal is defined only once *N* reaches
the estimable model rank.

**(e)** Pilot-informed exact D-optimal augmentation (`AlgDesign::optFederov`
with protected pilot runs); the cell shows oracle regret (oracle minus
pilot-estimated D-optimal MCC) and dots mark cells where pilot D-optimal still
beats random.
