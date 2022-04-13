# EQcycle_polylinear
Scripts used for solving the equations and producing the figures in Bolotskaya and Hager 2022 (in prep)

## Background
1D spring-slider model is often used to simulate earthquake cycles. Failure law prescribed between the block and the rough surface plays an important role in the system behavior. This repository contains MATLAB scripts that solve the equation of motion for 1D spring-slider system: 
![1](https://user-images.githubusercontent.com/11836119/152701882-46d69950-61e6-4497-aaa0-c372362e27fa.png)
with a poly-linear failure law (double slip-weakening with initial strengthening - DSWIS):
![2](https://user-images.githubusercontent.com/11836119/152701887-23e07e9a-df1e-4040-955a-b53654c90cd4.png)
The slope of the k<sup>th</sup> failure law segment is defined as:
![3](https://user-images.githubusercontent.com/11836119/152701898-06f19781-d1c5-4f72-aec2-535ec946b251.png)
Depending on the relation between *K<sub>k</sub><sup>f</sup>* (3) and *K* (the loading slope or spring stiffness), eq. (1) has three solution regimes: harmonic oscillations, cubic growth solution, and exponential growth solution. 
By combining different slopes within one failure law and thus different stability regimes, we could observe various earthquake cycle related behaviors:
-	aseismic oscillatory creep 
-	precursor-like events prior to the main-shock 
-	different main-shock characteristics 
-	tremor-like events instead of the main shock

## Repository contents
- MATLAB scripts:
  - *Analytical_solutions_generic_eqn_linear_friction*
  - *Analytical solutions_DSWIS* 
  - *DSWIS_analytical_solutions_and_plots* 
  - *DSWIS_multiple_analytical_solutions_and_plots*
- README.md
- LICENSE

### *Analytical_solutions_generic_eqn_linear_friction*
Analytically solves non-dimensional eq. (1) with a generic linear friction segment, shows the three solution regimes: *K<sup>f</sup> < K* - harmonic oscillations, *K<sup>f</sup> = K* - cubic growth solution, and *K<sup>f</sup> > K* - exponential growth solution.
### *Analytical solutions_DSWIS*
Analytically solves the 1D dynamic spring slider equation with a double slip weakening with initial strengthening (DSWIS) failure law (2). Shows full analytic solutions for each segment with initial conditions from the previous segment.
### *DSWIS_analytical_solutions_and_plots*
Mostly analytically (the equation to find the duration of different phases does not have analytical solutions, thus we solve for them numerically) solves eq. (1) with DSWIS (2) and produces plots for a single set of failure law parameters: energy curves, phase diagrams, slip rate and slip vs. time for several earthquake cycles, slip and slip rate plot for different phases separately.
### *DSWIS_multiple_analytical_solutions_and_plots*
Mostly analytically (same as above) solves eq. (1) with DSWIS (2) and produces plots for several sets of failure law parameters (with the same axis scales) for comparison. Different sets of parameters (3 to 6 failure laws) are given as examples.

## Reference
Please refer the following article if you use EQcycle_polylinear for your research:
Bolotskaya and Hager 2022 (submitted to BSSA)  
Please contact me (bolee@mit.edu) if you decide to use it while it is still "submitted".

[![DOI](https://zenodo.org/badge/434003826.svg)](https://zenodo.org/badge/latestdoi/434003826)
