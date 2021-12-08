# EQcycle_polylinear
## Background
1D spring-slider model is often used to simulate earthquake cycles. This repository contains MATLAB scripts that solve the equation of motion for 1D spring-slider system: 
![1_eq](https://user-images.githubusercontent.com/11836119/145135061-674b0371-9e54-4d69-8f0c-49f2e5f18a47.png)
with a poly-linear failure law (double slip-weakening with initial strengthening - DSWIS):
![2_eq](https://user-images.githubusercontent.com/11836119/145135080-04c2f49f-be9d-417d-aee1-5d49d15285f0.png)
The slope of the failure law segment:
![3_eq](https://user-images.githubusercontent.com/11836119/145135094-52f145d8-d61a-481a-b878-286e91e2daaa.png)

## Repository contents
- MATLAB scripts:
  - *Analytical_solutions_generic_eqn_linear_friction* - analytically solves the 1D dynamic spring slider equation with a generic linear friction segment: 
  - *DSWIS_analytical_solutions_fully_symbolic* - analytically solves the 1D dynamic spring slider equation with a double slip weakening with initial strengthening (DSWIS) failure law 
  - *DSWIS_analytical_solutions_and_plots* - mostly analytically (some variables are solved for numerically) solves the 1D dynamic spring slider equation with DSWIS and produces plots for a single set of failure law parameters
  - *DSWIS_multiple_analytical_solutions_and_plots* - mostly analytically (some variables are solved for numerically) solves the 1D dynamic spring slider equation with DSWIS and produces plots for several sets of failure law parameters (with the same scales) for comparison
- README.md
- LICENSE
