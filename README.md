# EQcycle_polylinear
Scripts used for solving the equations and produsing the figures in Bolotskaya and Hager 2022 (in prep)

## Background
1D spring-slider model is often used to simulate earthquake cycles. Failure law prescribed between the block and the rough surface plays an important role in the system behavior. This repository contains MATLAB scripts that solve the equation of motion for 1D spring-slider system: 
![1_eq_4](https://user-images.githubusercontent.com/11836119/145231080-5b7ffb60-c379-4597-8c5a-22aca580de66.png)
with a poly-linear failure law (double slip-weakening with initial strengthening - DSWIS):
![2_eq_4](https://user-images.githubusercontent.com/11836119/145231713-ec661ced-ea22-457a-8d73-8f489df6f232.png)
The slope of the failure law segment is defined as:
![3_eq_4](https://user-images.githubusercontent.com/11836119/145231745-baf84a82-6a2b-44df-bc3c-13b4c908c325.png)
Depending on the relation between *K<sub>seg</sub>* (3) and *K* (the loading slope or spring stiffness), eq. (1) has three solution regimes: harmonic oscillations, cubic growth solution, and exponential growth solution. 
By combining different slopes within one failure law and thus different stability regimes, we could observe various earthquake cycle related behaviors:
-	aseismic oscillatory creep 
-	precursor-like events prior to the main-shock 
-	different main-shock characteristics 
-	tremor-like events instead of the main shock

## Repository contents
- MATLAB scripts:
  - *Analytical_solutions_generic_eqn_linear_friction*
  - *DSWIS_analytical_solutions_fully_symbolic* 
  - *DSWIS_analytical_solutions_and_plots* 
  - *DSWIS_multiple_analytical_solutions_and_plots*
- README.md
- LICENSE

### *Analytical_solutions_generic_eqn_linear_friction*
Analytically solves eq. (1) with a generic linear friction segment, shows the three solution regimes: *K<sub>seg</sub> < K* - harmonic oscillations, *K<sub>seg</sub> = K* - cubic growth solution, and *K<sub>seg</sub> > K* - exponential growth solution.
### *DSWIS_analytical_solutions_fully_symbolic*
Analytically solves the 1D dynamic spring slider equation with a double slip weakening with initial strengthening (DSWIS) failure law (2). Shows full analytic solutions for each segment with initial conditions from the previous segment.
### *DSWIS_analytical_solutions_and_plots*
Mostly analytically (the equation to find the duration of different phases does not have analytical solutions, thus we solve for them numerically) solves eq. (1) with DSWIS (2) and produces plots for a single set of failure law parameters: energy curves, phase diagrams, slip rate and slip vs. time for several earthquake cycles, slip and slip rate plot for different phases separately.
### *DSWIS_multiple_analytical_solutions_and_plots*
Mostly analytically (same as above) solves eq. (1) with DSWIS (2) and produces plots for several sets of failure law parameters (with the same axis scales) for comparison. Different sets of parameters (3 to 6 failure laws) are given as examples.

## Reference
Please refer the following article if you use EQcycle_polylinear for your research:
Bolotskaya and Hager 2022 (in prep)  
Please contact me (bolee@mit.edu) if you decide to use it while it is still "in prep".
