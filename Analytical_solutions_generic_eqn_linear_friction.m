clearvars;
close all;

%% Analytical solutions for the generic equation (fully symbolic)
% by Ekaterina Bolotskaya

% 11/29/2021

% This script solves the equation of motion for 1D spring-slider:
% y" + A*y = B*t + C
% where y - slip
% y" - acceleration
% t - time
% A, B, C - constants

% There are three "stability regimes" 
% and three solution types depending on the sign of A

%% Equation
syms y(t) 
Dy = diff(y);

syms A real 
syms C B v_prev real 
% v_prev - slip rate at t = 0

ode = diff(y,t,2) + A*y == B*t + C;
cond1 = y(0) == 0;
cond2 = Dy(0) == v_prev;

%% A > 0
assumeAlso(A > 0)
% Solve
ySol(t) = dsolve(ode,[cond1 cond2]);

% Display
fprintf('A > 0 \n\n y(t) = \n\n');
pretty(simplify(ySol))

fprintf('y''(t) = \n\n');
pretty(simplify(diff(ySol)))

fprintf('y"(t) = \n\n');
pretty(simplify(diff(diff(ySol))))

%% A < 0
syms A real 
assumeAlso(A < 0)
% Solve
ySol(t) = dsolve(ode,[cond1 cond2]);

% Display
fprintf('A < 0 \n\n y(t) = \n\n');
pretty(simplify(ySol))

fprintf('y''(t) = \n\n');
pretty(simplify(diff(ySol)))

fprintf('y"(t) = \n\n');
pretty(simplify(diff(diff(ySol))))

%% A = 0
% New equation
ode = diff(y,t,2) == B*t + C;

% Solve
ySol(t) = dsolve(ode,[cond1 cond2]);

% Display
fprintf('A = 0 \n\n y(t) = \n\n');
pretty(simplify(ySol))

fprintf('y''(t) = \n\n');
pretty(simplify(diff(ySol)))

fprintf('y"(t) = \n\n');
pretty(simplify(diff(diff(ySol))))

