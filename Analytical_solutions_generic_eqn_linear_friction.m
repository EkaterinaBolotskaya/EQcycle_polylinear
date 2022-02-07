clearvars;
close all;

%% Analytical solutions for the generic equation (fully symbolic)
% by Ekaterina Bolotskaya

% 12/30/2021

% This script solves the nondim equation of motion for 1D spring-slider:
% y" + (1-K_f/K)*y = V_0*t + d_tau_i
% where y - nondim slip
% y"      - nondim acceleration
% t       - nondim time
% K_f     - failure law segment slope
% K       - elastic loading slope (spring stiffness)
% V_0     - nondim load point velocity
% d_tau_i - nondim stress difference at the beginning of the segment
% Substitute Ke = (K-K_f)/K - effective loading slope

% There are three "stability regimes" 
% and three solution types depending on the sign of Ke

%% Equation
syms y(t) 
Dy = diff(y);

syms Ke real
syms V_0 d_tau_i V_i real 
% V_i - nondim slip rate at t = 0

ode = diff(y,t,2) + Ke*y == V_0*t + d_tau_i;
cond1 = y(0) == 0;
cond2 = Dy(0) == V_i;

%% Ke > 0
assumeAlso(Ke > 0)
% Solve
ySol(t) = dsolve(ode,[cond1 cond2]);

% Display
fprintf('K_f < K \n\n y(t) = \n\n');
pretty(simplify(ySol))

fprintf('y''(t) = \n\n');
pretty(simplify(diff(ySol)))

fprintf('y"(t) = \n\n');
pretty(simplify(diff(diff(ySol))))

%% Ke < 0
syms Ke real 
assumeAlso(Ke < 0)
% Solve
ySol(t) = dsolve(ode,[cond1 cond2]);

% Display
fprintf('K_f > K \n\n y(t) = \n\n');
pretty(simplify(ySol))

fprintf('y''(t) = \n\n');
pretty(simplify(diff(ySol)))

fprintf('y"(t) = \n\n');
pretty(simplify(diff(diff(ySol))))

%% Ke = 0
% New equation
ode = diff(y,t,2) == V_0*t + d_tau_i;

% Solve
ySol(t) = dsolve(ode,[cond1 cond2]);

% Display
fprintf('K_f = K \n\n y(t) = \n\n');
pretty(simplify(ySol))

fprintf('y''(t) = \n\n');
pretty(simplify(diff(ySol)))

fprintf('y"(t) = \n\n');
pretty(simplify(diff(diff(ySol))))

