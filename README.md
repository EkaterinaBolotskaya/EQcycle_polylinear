# EQcycle_polylinear
Scripts used for solving the equations and producing the figures in:

E. Bolotskaya and B.H. Hager; A 1D Spring‐Slider Model with a Simple Poly‐Linear Failure Law Produces Rich Variations in Slip Behavior. Bull. Seismol. Soc. Am. 2022; doi: https://doi.org/10.1785/0120220052

## Background
1D spring-slider model is often used to simulate earthquake cycles. Failure law prescribed between the block and the rough surface plays an important role in the system behavior. This repository contains MATLAB scripts that solve the equation of motion for 1D spring-slider system: 

$$
\ddot{D} = -\frac{1}{M}(K(D - V_{0} t) +\mu\sigma_n), \qquad \qquad \qquad \qquad (1) 
$$

where $M$ is the mass of the system, $K$ is the spring stiffness, $D$ is the slip, $V$ is the slip rate, $V_0$ is the load point velocity, $\sigma_n$ is the normal stress, and $\mu$ is the friction coefficient,

with a poly-linear failure law (double slip-weakening with initial strengthening - DSWIS):

$$\mu =
    \begin{cases}
      \mu_i - (\mu_i-\mu_s)\frac{D}{D_s}               & D \le D_s\\
      \mu_s - (\mu_s-\mu_t)\frac{D-D_s}{D_{w1}}        & D_s < D \le D_s+D_{w1}\\
      \mu_t - (\mu_t-\mu_d)\frac{D-D_s-D_{w1}}{D_{w2}} & D_s+D_{w1} < D \le D_s+D_{w1}+D_{w2}\\
      \mu_d                                            & D > D_s+D_{w1}+D_{w2}\\
    \end{cases}
    , \qquad \qquad \qquad \qquad (2)
$$

where $\mu_i$, $\mu_s$, $\mu_t$, and $\mu_d$ are the initial, static, transitional, and dynamic friction coefficients, $D_{w1}$ and $D_{w2}$ are the slip-weakening distances, and $D_s$ is the slip-strengthening distance. The last segment of the failure law is horizontal with $\mu=\mu_d$.

The slope of the $k^{th}$ failure law segment is defined as:

$$
K_k^f =  \frac{\tau_k-\tau_{k+1}}{D_k}, \qquad \qquad \qquad \qquad (3)
$$

where $\tau_k-\tau_{k+1}$ is the difference in shear stress between the beginning and the end of the segment and $D_k$ is the corresponding difference in slip. Because of the way that the stress difference is defined, “weakening” slopes are positive and “strengthening” slopes are negative.

Depending on the relation between $K_k^f$ (3) and $K$ (the loading slope or spring stiffness), equation (1) has three solution regimes: harmonic oscillations, cubic growth solution, and exponential growth solution. 
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
  - *Frequency_lower_bound_analytical*
- README.md
- LICENSE

### *Analytical_solutions_generic_eqn_linear_friction*
Analytically solves non-dimensional equation (1) with a generic linear friction segment, shows the three solution regimes: $K_k^f < K$ - harmonic oscillations, $K_k^f=K$ - cubic growth solution, and $K_k^f>K$ - exponential growth solution.
### *Analytical solutions_DSWIS*
Analytically solves the 1D dynamic spring slider equation with a double slip weakening with initial strengthening (DSWIS) failure law (2). Shows full analytic solutions for each segment with initial conditions from the previous segment.
### *DSWIS_analytical_solutions_and_plots*
Mostly analytically (the equation to find the duration of different phases does not have analytical solutions, thus we solve for them numerically) solves equation (1) with DSWIS (2) and produces plots for a single set of failure law parameters: energy curves, phase diagrams, slip rate and slip vs. time for several earthquake cycles, slip and slip rate plot for different phases separately, spectra.
### *DSWIS_multiple_analytical_solutions_and_plots*
Mostly analytically (same as above) solves equation (1) with DSWIS (2) and produces plots for several sets of failure law parameters (with the same axis scales) for comparison. Different sets of parameters (3 to 8 failure laws) are given as examples.
### *Frequency_lower_bound_analytical*
Estimates the lower bound on frequency of the oscillatoric solution of block slip with poly-linear friction for a range of fault lengths and slip-weakening distances, assuming the block goes through a single oscillation during the weakening process.

## Reference
Please refer the following article if you use EQcycle_polylinear for your research:

E. Bolotskaya and B.H. Hager; A 1D Spring‐Slider Model with a Simple Poly‐Linear Failure Law Produces Rich Variations in Slip Behavior. Bull. Seismol. Soc. Am. 2022; doi: https://doi.org/10.1785/0120220052

Release on Zenodo:

[![DOI](https://zenodo.org/badge/434003826.svg)](https://zenodo.org/badge/latestdoi/434003826)
