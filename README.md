# Systematic Conservation Planning in Swedish Forests

R code for the MSc thesis *Systematic Conservation Planning in Swedish Forests*
(Tijmen Hennekes, MSc Environmental Resource Management, Vrije Universiteit Amsterdam).
The full methodology, data sources, and results are described in the thesis; this
repository contains the scripts that produce them.

The analysis allocates Sweden's forest to three management zones (close-to-nature,
integrated, and managed forest) using multi-zone integer linear programming with the
[prioritizr](https://prioritizr.net/) package and the Gurobi solver.

## Scripts

| Script | Purpose |
|---|---|
| `Scenario0_baseline.R` | S0 – current management, evaluated as the baseline |
| `scenario1_policy.R` | S1 – optimisation at the current budget |
| `scenario2_feasible.R` | S2 – minimum-feasible solution meeting both targets |
| `scenario3_substitution.R` | S3 – optimisation with harvested-wood substitution |
| `sensitivity1_budget.R` | Budget sensitivity (Table 8) |
| `sensitivity2_fmcost.R` | IFM-cost sensitivity (Table 9) |
| `sensitivity3_mfweight.R` | MF-weight sensitivity (Table 10) |
| `transition_intensity_geotiff.R` | Transition map between S0 and S2 (Section 4.3) |


Requires R with `prioritizr`, `terra`, `sf`, `gurobi`, `tibble`, `readxl`, and `writexl`,
plus a Gurobi licence (free for academic use). The input datasets are not included; see
the thesis for their sources. The file paths at the top of each script are absolute and
must be adjusted to your own setup. Each script is self-contained and can be run on its own.
