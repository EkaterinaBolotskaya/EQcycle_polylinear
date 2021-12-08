clearvars;
close all;

%% Analytical solutions for DSWIS (fully symbolic)
% by Ekaterina Bolotskaya

% 06/15/2021

% This script solves the equation of motion for 1D spring-slider:
% a = -1/M*(K*(D-V_0*t)+mu*sigma_n)
% where a - acceleration (second derivative of slip (D))
% M - mass
% K - spring stiffness
% D - slip
% V_0 - load point velocity
% t - time
% mu - friction coefficient
% sigma_n - normal stress
% mu is determined by the failure law. In this case double slip-weakening 
% with initial strengthening (DSWIS) failure law, which consists of one
% linear strengthening and two linear weakening segments

%% Analytical solution for strengthening (1st segment)
syms us(t) 
Dus = diff(us);

syms M K V_0 d_tau_s D_s t real positive
% d_tau_s = Sn*(mu_s-mu_i) - strengthening stress increase
% D_s - slip strengthening distance

% Equation
ode = diff(us,t,2)*M == -K*(us - V_0*t) - d_tau_s*(us)/D_s;
% Initial conditions
cond1 = us(0) == 0;
cond2 = Dus(0) == 0;

% Solve
uSol(t) = dsolve(ode,[cond1 cond2]);

% Display slip
uSol = simplify(uSol, 'Steps', 20);
fprintf('Strengthening:\n u(t) = \n\n');
pretty(uSol)
% Display slip rate
vSol = diff(uSol);
fprintf('v(t) = \n\n');
pretty(vSol)
% Display acceleration
aSol = diff(vSol);
fprintf('a(t) = \n\n');
pretty(aSol)

%% Analytical solution for weakening 1 (2nd segment)
syms uw1(t) 
Duw1 = diff(uw1);

syms d_tau_w1 D_w1 real positive
syms t_s v_s real positive
% d_tau_w1 = Sn*(mu_s-mu_t) - first weakening stress drop
% D_w1 - first slip-weakening distance
% v_s - slip rate at the end of the strengthening segment
% t_s - time at the end of the strengthening segment

% Equation
ode = diff(uw1,t,2)*M == -K*(uw1 - V_0*t - V_0*t_s) + d_tau_w1*(uw1)/D_w1 - d_tau_s - K*D_s;
% Initial conditions
cond1 = uw1(0) == 0;
cond2 = Duw1(0) == v_s;

% Solve
uSol(t) = dsolve(ode,[cond1 cond2]);

% Display slip
uSol = simplify(uSol, 'Steps', 20);
fprintf('Strengthening:\n u(t) = \n\n');
pretty(uSol)
% Display slip rate
vSol = diff(uSol);
fprintf('v(t) = \n\n');
pretty(vSol)
% Display acceleration
aSol = diff(vSol);
fprintf('a(t) = \n\n');
pretty(aSol)

%% Analytical solution for weakening 2 (3rd segment)
syms uw2(t) 
Duw2 = diff(uw2);

syms d_tau_w2 D_w2 real positive
syms t_w1 v_w1 real positive
% d_tau_w2 = Sn*(mu_t-mu_d) - second weakening stress drop
% D_w2 - second slip-weakening distance
% v_w1 - slip rate at the end of the first weakening segment
% t_w1 - time at the end of the first weakening segment

% Equation
ode = diff(uw2,t,2)*M == -K*(uw2 - V_0*t - V_0*(t_s+t_w1)) + d_tau_w2*(uw2)/D_w2 - d_tau_s+d_tau_w1 - K*(D_s+D_w1);
% Initial conditions
cond1 = uw2(0) == 0;
cond2 = Duw2(0) == v_w1;

% Solve
uSol(t) = dsolve(ode,[cond1 cond2]);

% Display slip
uSol = simplify(uSol, 'Steps', 20);
fprintf('Strengthening:\n u(t) = \n\n');
pretty(uSol)
% Display slip rate
vSol = diff(uSol);
fprintf('v(t) = \n\n');
pretty(vSol)
% Display acceleration
aSol = diff(vSol);
fprintf('a(t) = \n\n');
pretty(aSol)

%% Analytical solution for flat part
syms uf(t)
Duf = diff(uf);
syms t_w2 v_w2 real positive
% v_w2 - slip rate at the end of the second weakening segment
% t_w2 - time at the end of the second weakening segment

% Equation
ode = diff(uf,t,2)*M == -K*(uf - V_0*t - V_0*(t_s+t_w1+t_w2)) - d_tau_s+d_tau_w1+d_tau_w2 - K*(D_s+D_w1+D_w2);
% Initial conditions
cond1 = uf(0) == 0;
cond2 = Duf(0) == v_w2;

% Solve
uSol(t) = dsolve(ode,[cond1 cond2]);

% Display slip
uSol = simplify(uSol, 'Steps', 20);
fprintf('Strengthening:\n u(t) = \n\n');
pretty(uSol)
% Display slip rate
vSol = diff(uSol);
fprintf('v(t) = \n\n');
pretty(vSol)
% Display acceleration
aSol = diff(vSol);
fprintf('a(t) = \n\n');
pretty(aSol)
