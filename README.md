# EQcycle_polylinear
## Background
1D spring-slider model is often used to simulate earthquake cycles. This repository contains MATLAB scripts that solve the equation of motion for 1D spring-slider system: 
![1_eq](https://user-images.githubusercontent.com/11836119/145133459-66b4cd4a-447a-46d4-80cc-bf8b718b6916.png)
with a poly-linear failure law (double slip-weakening with initial strengthening - DSWIS):
![2_eq](https://user-images.githubusercontent.com/11836119/145133837-886b4260-62cb-4b69-b276-87fff83b0ebb.png)
The slope of the failure law segment:
![3_eq](https://user-images.githubusercontent.com/11836119/145134347-aec04736-182a-45eb-9cf7-9424514dbd44.png)

## Repository contents
- MATLAB scripts:
  - *Analytical_solutions_generic_eqn_linear_friction* - analytically solves the 1D dynamic spring slider equation with a generic linear friction segment: 
  - *DSWIS_analytical_solutions_fully_symbolic* - analytically solves the 1D dynamic spring slider equation with a double slip weakening with initial strengthening (DSWIS) failure law 
  - *DSWIS_analytical_solutions_and_plots* - mostly analytically (some variables are solved for numerically) solves the 1D dynamic spring slider equation with DSWIS and produces plots for a single set of failure law parameters
  - *DSWIS_multiple_analytical_solutions_and_plots* - mostly analytically (some variables are solved for numerically) solves the 1D dynamic spring slider equation with DSWIS and produces plots for several sets of failure law parameters (with the same scales) for comparison
- README.md
