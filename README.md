# Labor and AI: A Simulation Exercise
 
This repository contains a **simulation exercise accompanied by a report** — not a research paper — that builds a synthetic data-generating process (DGP) to study how AI adoption reshapes regional labor market dynamics across high-, mid- and low-skill segments of the workforce.
 
## The Report
 
Simulation Report.pdf
 
## The Code
 
The exercise constructs two synthetic DGPs — a discrete (binary) AI adoption setting and a continuous (investment-dollar) treatment setting — for 1,000 regions observed over 20 years, then estimates the resulting treatment effects with four progressively more robust estimators to illustrate how staggered and continuous treatment timing biases naive difference-in-differences designs. All of the simulation, estimation and plotting is contained in a single script.
 
### Repository Structure:
 
- `Simulation Report.pdf`: The final report describing the DGP construction, the estimation methods, and the results.
- `Simulation.R`: The single script generating both DGPs and estimating all four specifications.
### Script Guide:
 
**Simulation.R**
 
- *Construction of the first (discrete) DGP: 1,000 regions over 20 years, group-specific Cobb-Douglas production functions for high-, mid- and low-skill labor, and staggered binary AI adoption tied to each region's infrastructural capacity.*
- *Baseline OLS and Two-Way Fixed Effects (TWFE) estimation, alongside a Goodman-Bacon (2021) decomposition quantifying the bias introduced by "forbidden comparisons" under staggered treatment timing.*
- *Sun & Abraham (2021) estimation producing bias-corrected, year-specific event-study coefficients and the resulting ATT.*
- *Construction of a second, continuous-treatment DGP, where AI investment — rather than a 0/1 indicator — feeds capital accumulation directly into the high-skill production function.*
- *A manually implemented Callaway-Goodman Bacon-Sant'Anna (CGBS, 2025) continuous dose-response estimator, recovering year-specific Average Causal Responses (ACRs) and an aggregate ATT evaluated at mean dosage.*
- *Side-by-side comparison of all four estimators (OLS, TWFE, SA, CGBS) across employment, wages, output and labor share, for both the discrete and continuous treatment settings.*


**Tools used:** R, Two-Way Fixed Effects (TWFE), Goodman-Bacon (2021) Decomposition, Sun & Abraham (2021) Estimator, Continuous-Treatment Difference-in-Differences (Callaway-Goodman Bacon-Sant'Anna, 2025).
 
---
 
*This project was conducted as the Final Project Report for Research Strategy (EBC4125) at Maastricht University (Academic Year 2025-2026).*
 
