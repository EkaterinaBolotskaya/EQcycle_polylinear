clearvars;
close all;

%% Analytical solutions for DSWIS and plots
% by Ekaterina Bolotskaya

% 08/23/2021

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

%% System parameters
Sn       = 20e6;                     % normal stress 
M        = 80e6;                     % mass
K_s      = 1e6;                      % spring stiffness
V_0      = 1e-3;                     % velocity of the load point
int_vel  = 1e-7;                     % interseismic velocity

%% Failure law parameters
mu_i     = 0.69;                     % initial friction coefficient
mu_s     = 0.70;                     % static (max) friction coefficient
mu_d     = 0.60;                     % transition friction coefficient
mu_t     = 0.687;                    % dynamic friction coefficient
D_s      = 0.1052631579;             % strengthening distance
D_t      = 0.4834759358;             % weakening 1 distance
D_w      = 0.57;                     % weakening 2 distance

%% Slopes and energies
K_str  = Sn*(mu_i-mu_s)/D_s;
K_w1   = Sn*(mu_s-mu_t)/D_t;
K_w2   = Sn*(mu_t-mu_d)/D_w;
en_s   = Sn*((mu_s-mu_i)*D_s/2+(mu_i-mu_d)*D_s);
en_w   = Sn*((mu_s-mu_t)*D_t/2+(mu_t-mu_d)*D_w/2+(mu_t-mu_d)*D_t);

fprintf('K_str = %.0f K_w1 = %.0f K_w2 = %.0f \n', K_str, K_w1, K_w2);
fprintf('e_w = %.0f e_s = %.0f \n', en_w, en_s);

%% Plotting and supplementary variables
lw       = 0.75;                      % line width
fs       = 11;                       % font size
mpp      = 0.01*Sn;                  % min peak prominence for peak detection
v_pc     = 10*V_0;                    % preseismic to coseismic slip rate threshold
k        = 3;                        % number of cycles

% Colors
b_col    = [0.231372549019608 0.298039215686275 0.752941176470588];
r_col    = [0.705882352941177 0.015686274509804 0.149019607843137];
g_col    = [0 .7 0];

% Nondimentional variables for plotting
nd_t     = sqrt(M/K_s);                     % time
nd_u     = Sn*(mu_s-mu_d)/K_s;              % slip
nd_v     = nd_u/nd_t;                       % slip rate
nd_a     = nd_v/nd_t;                       % acceleration

%% Analytical solution for strengthening (segment 1)
syms u(t)
Du       = diff(u);

% Equation
ode      = diff(u,t,2)*M == -K_s*(u-V_0*t)-Sn*(mu_s-mu_i)*(u)/D_s;
% Initial conditions
cond1    = u(0) == 0;
cond2    = Du(0) == int_vel;

uSol(t)  = dsolve(ode,[cond1 cond2]);
uSol     = simplify(uSol);
vSol     = diff(uSol);
aSol     = diff(vSol);

% Cannot solve for time analytically, thus, solve numerically
t_s      = double(vpasolve(uSol == D_s, t));

st1      = pi*sqrt(D_s*M/(Sn*(mu_s-mu_i)+D_s*K_s))/40;
time_s   = double(0:st1:t_s);
disp_s   = double(subs(uSol, t, time_s));
vel_s    = double(subs(vSol, t, time_s));
ac_s     = double(subs(aSol, t, time_s));
tau_s    = double(subs(Sn*(mu_s-mu_i)*(uSol)/D_s, t, time_s));
ku_s     = double(subs(-K_s*(uSol-V_0*t), t, time_s));

%% Analytical solution for weakening 1 (segment 2)
syms u(t)
Du       = diff(u);

% Equation
ode      = diff(u,t,2)*M == -K_s*(u-V_0*t-V_0*t_s)-Sn*(mu_t-mu_s)*(u)/D_t - Sn*(mu_s-mu_i) - K_s*D_s;
% Initial conditions (from segment 1)
cond1    = u(0) == 0;
cond2    = Du(0) == vel_s(end);

uSol(t)  = dsolve(ode,[cond1 cond2]);
uSol     = simplify(uSol);
vSol     = diff(uSol);
aSol     = diff(vSol);

% Cannot solve for time analytically, thus, solve numerically
t_w1     = double(vpasolve(uSol == D_t, t, [0 Inf]));

st2      = pi*sqrt(D_t*M/(Sn*abs(mu_s-mu_t)+D_t*K_s))/40;
time_w1  = double(0:st2:t_w1);
disp_w1  = disp_s(end)+double(subs(uSol, t, time_w1));
vel_w1   = double(subs(vSol, t, time_w1));
ac_w1    = double(subs(aSol, t, time_w1));
tau_w1   = double(subs(Sn*(mu_t-mu_s)*(uSol)/D_t + Sn*(mu_s-mu_i), t, time_w1));
ku_w1    = double(subs(-K_s*(uSol-V_0*t-V_0*t_s)- K_s*D_s, t, time_w1));

time_w1  = time_s(end) + time_w1;

%% Analytical solution for weakening 2 (segment 3)
syms u(t)
Du       = diff(u);

% Equation
ode      = diff(u,t,2)*M == -K_s*(u-V_0*t-V_0*(t_s+t_w1))-Sn*(mu_d-mu_t)*(u)/D_w - Sn*(mu_t-mu_i) - K_s*(D_s+D_t);
% Initial conditions (from segment 2)
cond1    = u(0) == 0;
cond2    = Du(0) == vel_w1(end);

uSol(t)  = dsolve(ode,[cond1 cond2]);
uSol     = simplify(uSol);
vSol     = diff(uSol);
aSol     = diff(vSol);

% Cannot solve for time analytically, thus, solve numerically
t_w2     = double(vpasolve(uSol == D_w, t));

N        = 1000;
time_w2  = double(linspace(0, t_w2, N));
disp_w2  = disp_w1(end)+double(subs(uSol, t, time_w2));
vel_w2   = double(subs(vSol, t, time_w2));
ac_w2    = double(subs(aSol, t, time_w2));
tau_w2   = double(subs(Sn*(mu_d-mu_t)*(uSol)/D_w + Sn*(mu_t-mu_i), t, time_w2));
ku_w2    = double(subs(-K_s*(uSol-V_0*t-V_0*(t_s+t_w1))-K_s*(D_s + D_t), t, time_w2));

time_w2  = time_w1(end) + time_w2;

%% Analytical solution for flat part (segment 4)
syms u(t)
Du       = diff(u);

% Equation
ode      = diff(u,t,2)*M == -K_s*(u-V_0*(t)-V_0*(t_s+t_w1+t_w2)) - Sn*(mu_d-mu_i)-K_s*(D_s+D_t+D_w);
% Initial conditions (from segment 3)
cond1    = u(0) == 0;
cond2    = Du(0) == vel_w2(end);

uSol(t)  = dsolve(ode,[cond1 cond2]);
uSol     = simplify(uSol);
vSol     = diff(uSol);
aSol     = diff(vSol);

% Cannot solve for time analytically, thus, solve numerically
t_fl     = vpasolve(vSol == int_vel, t, [0 pi*sqrt(M/K_s)]);

N        = 1000;
time_fl  = double(linspace(0, t_fl, N));
disp_fl  = disp_w2(end)+double(subs(uSol, t, time_fl));
vel_fl   = double(subs(vSol, t, time_fl));
ac_fl    = double(subs(aSol, t, time_fl));
tau_fl   = double(subs(Sn*(mu_d-mu_i), t, time_fl));
ku_fl    = double(subs(-K_s*(uSol-V_0*(t)-V_0*(t_s+t_w1+t_w2))-K_s*(D_s+D_t+D_w), t, time_fl));

time_fl  = time_w2(end) + time_fl;

%% Analytical solution for early interseismic part (segment 5)
% Equation
t_in     = double(solve(K_s*(V_0*t) == 0 - ku_fl(end), t));

N        = 1000;
time_in  = double(linspace(0, t_in, N));
vel_in   = int_vel*ones(size(time_in));
disp_in  = disp_fl(end)+int_vel*time_in;
ac_in    = zeros(size(time_in));
ku_in    = double(-K_s*(disp_in - V_0*time_in - V_0*(t_s+t_w1+t_w2+t_fl)));
tau_in   = ku_in;

time_in  = time_fl(end) + time_in;

%% Several cycles patch
% One full cycle
time_1c  = [time_s time_w1 time_w2 time_fl time_in];
disp_1c  = [disp_s disp_w1 disp_w2 disp_fl disp_in]-disp_s(1);
vel_1c   = [vel_s vel_w1 vel_w2 vel_fl vel_in];
ac_1c    = [ac_s ac_w1 ac_w2 ac_fl ac_in];
tau_1c   = [tau_s tau_w1 tau_w2 tau_fl tau_in];
ku_1c    = [ku_s ku_w1 ku_w2 ku_fl ku_in];

%% Spectra
% One full cycle different order
time_1c2 = [time_in - min(time_in) [time_s time_w1 time_w2 time_fl]+max(time_in-min(time_in))];
disp_1c2 = [disp_in - min(disp_in) [disp_s disp_w1 disp_w2 disp_fl]+max(disp_in-min(disp_in))]-disp_s(1);
vel_1c2  = [vel_in vel_s vel_w1 vel_w2 vel_fl];
ac_1c2   = [ac_in ac_s ac_w1 ac_w2 ac_fl];
tau_1c2  = [tau_in tau_s tau_w1 tau_w2 tau_fl];
ku_1c2   = [ku_in ku_s ku_w1 ku_w2 ku_fl];

% Several cycles
time_f   = repmat(time_1c, 1, k) + repelem(time_1c(end)*linspace(0,k-1,k),length(time_1c));
disp_f   = repmat(disp_1c, 1, k) + repelem(disp_1c(end)*linspace(0,k-1,k),length(disp_1c));
vel_f    = repmat(vel_1c, 1, k);
ac_f     = repmat(ac_1c, 1, k);
tau_f    = repmat(tau_1c, 1, k);
ku_f     = repmat(ku_1c, 1, k);

%% Phase transition indices
% Late interseismic to preseismic
idx_lip    = find(ku_1c == max(ku_1c), 1, 'last');
% Preseimic to coseismic
[~,idx_pc] = find(vel_1c >= v_pc, 1, 'first');
if isempty(idx_pc)
    idx_pc = find(vel_1c ~= int_vel, 1, 'last') + 1;
end
% Coseismic to early interseismic
idx_cei    = find(vel_1c ~= int_vel, 1, 'last') + 1;

%% Late interseismic phase
time_li  = time_1c(1:idx_lip);
disp_li  = disp_1c(1:idx_lip);
vel_li   = vel_1c(1:idx_lip);
ac_li    = ac_1c(1:idx_lip);
tau_li   = tau_1c(1:idx_lip);
ku_li    = ku_1c(1:idx_lip);

%% Preseismic phase
time_pr  = time_1c(idx_lip:idx_pc)-time_1c(idx_lip);
disp_pr  = disp_1c(idx_lip:idx_pc)-disp_1c(idx_lip);
vel_pr   = vel_1c(idx_lip:idx_pc);
ac_pr    = ac_1c(idx_lip:idx_pc);
tau_pr   = tau_1c(idx_lip:idx_pc);
ku_pr    = ku_1c(idx_lip:idx_pc);

%% Coseimic phase
time_co  = time_1c(idx_pc:idx_cei)-time_1c(idx_pc);
disp_co  = disp_1c(idx_pc:idx_cei)-disp_1c(idx_pc);
vel_co   = vel_1c(idx_pc:idx_cei);
ac_co    = ac_1c(idx_pc:idx_cei);
tau_co   = tau_1c(idx_pc:idx_cei);
ku_co    = ku_1c(idx_pc:idx_cei);

%% Plotting
%% Stress vs. slip (several cycles)
% Choose which peak number to color
pn         = 2;
tau_d      = Sn*(mu_d-mu_i);
idx_eq     = tau_f <= tau_d;
u_abs_eq   = NaN(size(disp_f));
u_abs_eq(idx_eq) = disp_f(idx_eq);

idx_eq_ku     = ku_f <= tau_d;
u_abs_eq_ku   = NaN(size(disp_f));
u_abs_eq_ku(idx_eq_ku) = disp_f(idx_eq_ku);

% Look for peaks and crossings for energy coloring, interpolate
[~,locs_t] = findpeaks(tau_f, linspace(1,length(ku_f),length(ku_f)),'MinPeakProminence', mpp);

[~, idx_w] = min(abs(u_abs_eq - disp_f(locs_t(pn))));
u_int      = interp1([tau_f(idx_w) tau_f(idx_w+1)], [disp_f(idx_w) disp_f(idx_w+1)], tau_d, 'linear');

idx_s      = find((u_abs_eq - disp_f(locs_t(pn)))>0, 1, 'first');
idx_s_ku   = find((u_abs_eq_ku - disp_f(locs_t(pn)))>0, 1, 'first');
idx_ini    = locs_t(pn) - find(ac_f(locs_t(pn):-1:idx_w) == 0, 1, 'first');

figure();
% Fill energies
xb   = (disp_f(locs_t(pn):idx_s)/nd_u)';
y1b  = ((tau_f(locs_t(pn):idx_s)-tau_d)/Sn/(mu_s-mu_d))';
y2b  = (tau_d*ones(size(y1b))-tau_d)/Sn/(mu_s-mu_d);
fill([xb; flipud(xb)], [y2b; flipud(y1b)], 'b', 'EdgeColor', [.5 .5 .5], 'LineWidth', 0.1, 'FaceAlpha', 0.35);
hold on;
xr   = [u_int disp_f(idx_w+1:locs_t(pn))]'/nd_u;
y1r  = [tau_d-tau_d tau_f(idx_w+1:locs_t(pn))-tau_d]'/Sn/(mu_s-mu_d);
y2r  = (tau_d*ones(size(y1r))-tau_d)/Sn/(mu_s-mu_d);
fill([xr; flipud(xr)], [y2r; flipud(y1r)], 'r', 'EdgeColor', [.5 .5 .5], 'LineWidth', 0.1, 'FaceAlpha', 0.35);
% Plot Stress vs. slip curves
j(1) = plot(disp_f/nd_u, (tau_f/Sn+mu_i-mu_d)/(mu_s-mu_d), 'k', 'LineWidth', lw);
j(2) = plot(disp_f/nd_u, (ku_f/Sn+mu_i-mu_d)/(mu_s-mu_d), 'k--', 'LineWidth', lw);
legend([j(1) j(2)],'Friction','Loading', 'Interpreter', 'latex', 'FontSize', fs, 'location', 'southwest');
ylabel('$\tau_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
xlabel('$D_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
grid on;
axis tight;
hold off;

%% Stress vs. slip (zoom)
figure();
% Fill energies
fill([xb; flipud(xb)], [y2b; flipud(y1b)], 'b', 'EdgeColor', [.5 .5 .5], 'LineWidth', 0.1, 'FaceAlpha', 0.35);
hold on;
fill([xr; flipud(xr)], [y2r; flipud(y1r)], 'r', 'EdgeColor', [.5 .5 .5], 'LineWidth', 0.1, 'FaceAlpha', 0.35);
% Plot Stress vs. slip curves
j(1) = plot(disp_f/nd_u, (tau_f/Sn+mu_i-mu_d)/(mu_s-mu_d), 'k', 'LineWidth', lw);
j(2) = plot(disp_f/nd_u, (ku_f/Sn+mu_i-mu_d)/(mu_s-mu_d), 'k--', 'LineWidth', lw);
legend([j(1) j(2)],'Friction','Loading', 'Interpreter', 'latex', 'FontSize', fs, 'location', 'northeast');
ylabel('$\tau_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
xlabel('$D_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
ylim([-0.03 (max(ku_li)/Sn+mu_i-mu_d)/(mu_s-mu_d)]*1.02);
xlim([min(xr)*0.95 max(xb)*1.25]);
grid on;
hold off;

%% Phase trajectories (several cycles) fails if int_vel == 0
figure();
p(1) = plot(log(vel_1c(idx_cei:end)/nd_v), (tau_1c(idx_cei:end)/Sn+mu_i-mu_d)/(mu_s-mu_d), 'k', 'LineWidth', lw);
hold on;
p(2) = plot(log(vel_li/nd_v), (tau_li/Sn+mu_i-mu_d)/(mu_s-mu_d), 'Color', r_col, 'LineWidth', lw);
p(3) = plot(log(vel_pr/nd_v), (tau_pr/Sn+mu_i-mu_d)/(mu_s-mu_d), 'Color', g_col, 'LineWidth', lw);
p(4) = plot(log(vel_co/nd_v), (tau_co/Sn+mu_i-mu_d)/(mu_s-mu_d), 'Color', b_col, 'LineWidth', lw);

plot(log(vel_1c(idx_cei:end)/nd_v), (ku_1c(idx_cei:end)/Sn+mu_i-mu_d)/(mu_s-mu_d), 'k--', 'LineWidth', lw);
plot(log(vel_li/nd_v), (ku_li/Sn+mu_i-mu_d)/(mu_s-mu_d), '--', 'Color', r_col, 'LineWidth', lw);
plot(log(vel_pr/nd_v), (ku_pr/Sn+mu_i-mu_d)/(mu_s-mu_d), '--', 'Color', g_col, 'LineWidth', lw);
plot(log(vel_co/nd_v), (ku_co/Sn+mu_i-mu_d)/(mu_s-mu_d), '--', 'Color', b_col, 'LineWidth', lw);

plot(log(v_pc/nd_v)*ones(1,100), 1.02*(linspace(min(tau_1c), max(ku_1c), 100)/Sn+mu_i-mu_d)/(mu_s-mu_d), 'k-.', 'LineWidth', lw-0.5);
plot(log(V_0/nd_v)*ones(1,100), 1.02*(linspace(min(tau_1c), max(ku_1c), 100)/Sn+mu_i-mu_d)/(mu_s-mu_d), 'k-.', 'LineWidth', lw-0.5);
grid on;
legend(p, {'Early interseismic', 'Late interseismic', 'Preseismic', 'Coseismic'}, 'Interpreter', 'latex', 'FontSize', fs, 'location', 'southwest');
xlabel('$ln(V_{nd})$', 'Interpreter', 'latex', 'FontSize', fs);
axis tight;
ylabel('$\tau_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);

%% Slip rate and slip in one plot (several cycles)
fig1 = figure();
set(fig1, 'defaultAxesColorOrder',[r_col; b_col]);

yyaxis left;
plot(time_f/nd_t, vel_f/nd_v, 'LineWidth', lw);
ylabel('$V_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
set(gca,'FontSize',fs);
axis tight;
grid on;

yyaxis right;
p(1) = plot(time_f/nd_t, disp_f/nd_u, 'LineWidth', lw);
hold on;
p(2) = plot(time_f/nd_t, (V_0*(time_f-max(time_co) - max(time_pr) - max(time_li)) + max(disp_co) + max(disp_pr) + max(disp_li))/nd_u, 'LineWidth', lw);
legend([p(1) p(2)], 'Slider','Load point', 'Interpreter', 'latex', 'FontSize', fs, 'location', 'southeast');
ylabel('$D_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
xlabel('$t_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);

%% Late interseismic phase
fig21 = figure();
set(fig21, 'defaultAxesColorOrder',[r_col; b_col]);

% Stress and acceleration vs. time
subplot(211)
yyaxis right;
p(2) = plot(time_li/nd_t, (tau_li/Sn+mu_i-mu_d)/(mu_s-mu_d), 'LineWidth', lw);
hold on;
p(3) = plot(time_li/nd_t, (ku_li/Sn+mu_i-mu_d)/(mu_s-mu_d), 'LineWidth', lw);
title({'Late interseismic'}, 'Interpreter', 'latex', 'FontSize', fs+1);
ylabel('$\tau_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
set(gca,'xticklabel',[]);
axis tight;

yyaxis left;
p(1) = plot(time_li/nd_t, ac_li/nd_a, 'LineWidth', lw);
ylabel('$a_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
set(gca,'FontSize',fs);
axis tight;
grid on;
legend([p(2) p(3)],'Friction','Loading', 'Interpreter', 'latex', 'FontSize', fs, 'location', 'southwest');

% Slip and slip rate vs. time
subplot(212)
yyaxis right;
plot(time_li/nd_t, disp_li/nd_u, 'LineWidth', lw);
ylabel('$D_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
xlabel('$t_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);

yyaxis left;
plot(time_li/nd_t, vel_li/nd_v, 'LineWidth', lw);
ylabel('$V_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
set(gca,'FontSize',fs);
axis tight;
grid on;

%% Preseismic phase
fig22 = figure();
set(fig22, 'defaultAxesColorOrder',[r_col; b_col]);

% Stress and acceleration vs. time
subplot(211)
yyaxis right;
p(2) = plot(time_pr/nd_t, (tau_pr/Sn+mu_i-mu_d)/(mu_s-mu_d), 'LineWidth', lw);
hold on;
p(3) = plot(time_pr/nd_t, (ku_pr/Sn+mu_i-mu_d)/(mu_s-mu_d), 'LineWidth', lw);
title({'Preseismic'}, 'Interpreter', 'latex', 'FontSize', fs+1);
ylabel('$\tau_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
set(gca,'xticklabel',[]);
axis tight;

yyaxis left;
p(1) = plot(time_pr/nd_t, ac_pr/nd_a, 'LineWidth', lw);
ylabel('$a_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
set(gca,'FontSize',fs);
axis tight;
grid on;
legend([p(2) p(3)],'Friction','Loading', 'Interpreter', 'latex', 'FontSize', fs-2, 'location', 'southwest');

% Slip and slip rate vs. time
subplot(212)
yyaxis right;
plot(time_pr/nd_t, disp_pr/nd_u, 'LineWidth', lw);
ylabel('$D_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
xlabel('$t_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);

yyaxis left;
plot(time_pr/nd_t, vel_pr/nd_v, 'LineWidth', lw);
ylabel('$V_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
set(gca,'FontSize',fs);
axis tight;
grid on;

%% Coseismic phase
% Find index of the beginning of the flat segment
[~,idx_fl] = find((tau_co/Sn+mu_i-mu_d)/(mu_s-mu_d)<1e-13, 1, 'first');

fig23 = figure();
set(fig23, 'defaultAxesColorOrder',[r_col; b_col]);

% Stress and acceleration vs. time
subplot(211)
yyaxis right;
p(2) = plot(time_co/nd_t, (tau_co/Sn+mu_i-mu_d)/(mu_s-mu_d), 'LineWidth', lw);
hold on;
p(3) = plot(time_co/nd_t, (ku_co/Sn+mu_i-mu_d)/(mu_s-mu_d), 'LineWidth', lw);
title({'Coseismic'}, 'Interpreter', 'latex', 'FontSize', fs+1);
ylabel('$\tau_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
set(gca,'xticklabel',[]);
axis tight;

yyaxis left;
p(1) = plot(time_co/nd_t, ac_co/nd_a, 'LineWidth', lw);
hold on;
xline((time_co(idx_fl))/nd_t, 'k-.', 'LineWidth', lw-0.5);
ylabel('$a_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
set(gca,'FontSize',fs);
axis tight;
grid on;
legend([p(2) p(3)],'Friction','Loading', 'Interpreter', 'latex', 'FontSize', fs-2, 'location', 'southwest');

% Slip and slip rate vs. time
subplot(212)
yyaxis right;
plot(time_co/nd_t, disp_co/nd_u, 'LineWidth', lw);
ylabel('$D_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
xlabel('$t_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);

yyaxis left;
plot(time_co/nd_t, vel_co/nd_v, 'LineWidth', lw);
hold on;
xline((time_co(idx_fl))/nd_t, 'k-.', 'LineWidth', lw-0.5);
ylabel('$V_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
set(gca,'FontSize',fs);
axis tight;
grid on;

%% Spectrum
vsp = [vel_s vel_w1(2:end) vel_w2(2:end) vel_fl(2:end)];
tsp = [time_s time_w1(2:end) time_w2(2:end) time_fl(2:end)];

% Calculate frequency and padded time
t_sp_max    = 50*max(tsp);
dta         = diff(tsp);
Fs          = 1/min(dta);
time_sp     = 0:min(dta):t_sp_max;
L           = length(time_sp);
if mod(L,2) == 1
    time_sp = 0:min(dta):t_sp_max + min(dta);
    L       = L+1;
end

% Compute FFT
v_sp_pad    = interp1(tsp+t_sp_max/2, vsp, time_sp, 'linear',0);
Y           = fft(v_sp_pad,L);
P2          = abs(Y/L);
P1          = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
spec        = P1(1:L/2);
    
% Plot spectra
figure(115);
loglog((0:(Fs/L):(Fs/2-Fs/L))*2*pi/sqrt(K_s/M),spec/max(spec), 'Color', r_col, 'LineWidth', lw, 'DisplayName', 'Spectrum');
hold on;
loglog(logspace(-1,5,100),1./(logspace(-2,4,100)), 'k:', 'LineWidth', lw/2, 'DisplayName', '$f_{nd}^{\ \ -1}$');
loglog(logspace(-1,5,100),1./(logspace(-2,4,100)).^2, 'k--', 'LineWidth', lw/2, 'DisplayName', '$f_{nd}^{\ \ -2}$');
set(gca,'FontSize',fs-2);
ylabel('Amplitude', 'Interpreter', 'latex', 'FontSize', fs);
xlabel('$f_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
legend('show', 'Interpreter', 'latex', 'FontSize', fs-1, 'location', 'southwest');
xlim([1e-2 1e3]);
ylim([1e-8 2]);
grid on;

