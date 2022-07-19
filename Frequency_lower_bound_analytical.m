clearvars;
close all;

%% Estimation of the lower bound on frequency of the oscillatoric solution 
% of block slip with poly-linear friction for a range of fault lengths and Dc
% by Ekaterina Bolotskaya

% 06/01/2022

% This script uses the oscillatory solution (K > Kf)  to the nondim equation of motion for 1D spring-slider:
% y" + (1-K_f/K)*y = V_0*t + d_tau_i
% where y - nondim slip
% y"      - nondim acceleration
% t       - nondim time
% K_f     - failure law segment slope
% K       - elastic loading slope (spring stiffness)
% V_0     - nondim load point velocity
% d_tau_i - nondim stress difference at the beginning of the segment

% to estimate the lower bound on oscillation frequency, assuming the block
% goes through a single oscillation during the weakening process.

%% Plotting and supplementary variables
lw       = 0.75;                      % line width
fs       = 11;                        % font size

%% Solve the equation, assuming K > Kf
syms uw(t)
Duw = diff(uw);

syms M K V_0 Kf t d_tau_i v_init real positive
assume(K > Kf);

ode     = diff(uw,t,2)*M == -K*(uw-V_0*t)+Kf*(uw)-d_tau_i;
cond1   = uw(0) == 0;
cond2   = Duw(0) == v_init;

uSol(t) = dsolve(ode,[cond1 cond2]);
uSol    = simplify(uSol, 'Steps', 20);
fprintf('u(t) = \n\n');
fprintf('Full solution: \n\n');
pretty(uSol)

% Displacement after one period
d_1per = subs(uSol, t, 2*pi/(sqrt((K-Kf)/M)));

fprintf('Slip of 1 period = \n\n');
pretty(simplify(d_1per))

%% Solve for Kf, if slip after one period  = Dc
syms d_tau D_c real positive
Kf1 = solve(d_1per == D_c, Kf, 'ReturnConditions', true);

%% Estimate frequency
fprintf('K_f = \n\n');
pretty(simplify(Kf1.Kf, 'Steps', 20))
fprintf('Frequency = \n\n');
syms pi
pretty(simplify(subs(1/2/pi*sqrt((K-Kf)/M), Kf, Kf1.Kf), 'Steps', 30))
clear pi

%% Fault and material parameters
G   = 30e9;                       % shear modulus
rho = 3000;                       % density
L   = logspace(0, 4.5, 50);       % fault patch length
D_c = logspace(-6, 1, 50);        % slip-weakening distance

[Lg, D_cg] = meshgrid(L, D_c);

M   = rho*Lg/pi/pi;               % effective mass (eqn. 14 from Im and Avouac, 2021) 
K   = G./Lg;                      % elastic stiffness

V_0 = 0.01/3.154e+7;              % loading rate

% Plot Kf
K_f = K - ((2*K.*M.^(1/2)*V_0*pi)./D_cg).^(2/3);

figure()
pcolor(L, D_c, K_f);
hold on;
set(gca, 'XScale', 'log', 'YScale', 'log','ColorScale','log');
set(gca,'Fontsize', fs-2)
shading interp;
c = colorbar;
xlabel('$L,\ m$', 'Interpreter', 'latex', 'FontSize', fs);
ylabel('$D_{c},\ m$', 'Interpreter', 'latex', 'FontSize', fs);
title(c, '$K_{f}$', 'Interpreter', 'latex', 'FontSize', fs);
box on;
hold off;

% Plot frequency
freq = 1/2/pi*((K-K_f)./M).^(1/2);

% Limits on stress drop
gh = log10(K_f.*D_cg);
gh(gh>8 | gh<3) = NaN;

fprintf('Minimum frequency with a reasonable stress drop = %e Hz \nMaximum frequency = %e Hz \n', min(min(freq(~isnan(gh)))), max(max(freq(~isnan(gh)))));

figure()
ax1 = axes;
pcolor(L, D_c, freq);
shading interp;
cb2 = colorbar(ax1,'Position',[.88 .10 .068 .808]);
title(cb2, '$f,\ Hz$', 'Interpreter', 'latex', 'FontSize', fs);
caxis([1e-4 10]);
set(gca, 'XScale', 'log', 'YScale', 'log','ColorScale','log');
xlabel('$L,\ m$', 'Interpreter', 'latex', 'FontSize', 11);
ylabel('$D_{c},\ m$', 'Interpreter', 'latex', 'FontSize', 11);
ax2 = axes;
loglog((K(1,:)), 1e4./(K_f(1,:)), 'k', 'Linewidth', lw);
hold on;
loglog((K(1,:)), 1e5./(K_f(1,:)), 'k', 'Linewidth', lw);
loglog((K(1,:)), 1e6./(K_f(1,:)), 'k', 'Linewidth', lw);
loglog((K(1,:)), 1e7./(K_f(1,:)), 'k', 'Linewidth', lw);
loglog((K(1,:)), 1e8./(K_f(1,:)), 'k', 'Linewidth', lw);
xlim([min(K(1,:)) max(K(1,:))]);
ylim([min(D_c) max(D_c)]);
ax2.XAxisLocation = 'top';
set(gca, 'XDir','reverse', 'color', 'none', 'Fontsize', fs-2)
xlabel('$K,\ Pa/m$', 'Interpreter', 'latex', 'FontSize', fs);
ax2.YTick = [];
set([ax1,ax2],'Position',[.15 .10 .685 .808]);
hpos = 5e9;
text(5e9,3e-2,'$\Delta\tau$ = 100\ MPa','rotation',29, 'Interpreter', 'latex', 'FontSize', fs);
text(2.2e9,7e-3,'10\ MPa','rotation',29, 'Interpreter', 'latex', 'FontSize', fs);
text(2.2e9,7e-4,'1\ MPa','rotation',29, 'Interpreter', 'latex', 'FontSize', fs);
text(2.2e9,7e-5,'0.1\ MPa','rotation',29, 'Interpreter', 'latex', 'FontSize', fs);
text(2.2e9,7e-6,'0.01\ MPa','rotation',29, 'Interpreter', 'latex', 'FontSize', fs);
box on;


