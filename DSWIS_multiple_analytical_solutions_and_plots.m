clearvars;
close all;

%% Multiple DSWIS failure laws comparison and plots (analytical)
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

% This script compares several (3 to 8 in the examples) failure laws
% Failure law parameters are given as arrays (examples are given below).
% Please comment all but one.
% Please make sure np = length of parameter arrays, if you modify.

%% System parameters
Sn       = 20e6;                     % normal stress 
M        = 60e6;                     % mass
K_s      = 1e6;                      % spring stiffness
V_0      = 1e-3;                     % velocity of the load point
int_vel  = 1e-7;                     % interseismic velocity

%% Failure law parameters
%% Original set from Bolotskaya and Hager 2022
% np       = 4;                        % number of failure laws (plots)                       
% io       = [0 0 0 0];                % set to 1 if you want individual plots
% 
% % initial friction coefficient
% mu_i_arr = [0.66 0.67 0.68 0.69];
% % static (max) friction coefficient
% mu_s_arr = 0.70*ones(1,np);
% % transition friction coefficient
% mu_t_arr = [0.65 0.675 0.68 0.687];
% % dynamic friction coefficient
% mu_d_arr = 0.60*ones(1,np); 
% % strengthening distance
% D_s_arr  = [0.125 0.11764705882352941176470588235294 0.11111111111111111111111111111111 0.10526315789473684210526315789474];
% % weakening 1 distance
% D_t_arr  = [0.5 0.5 0.51111111111111111111111111111111 0.4834759358288770053475935828877];                        
% % weakening 2 distance
% D_w_arr  = [1.3 0.7 0.6 0.57];   

%% Original set from Bolotskaya and Hager 2022 + 3 stable ones (Supplementary)
np       = 8;                                % number of failure laws (plots)                       
io       = [0 0 0 0 1 1 0 1];                % set to 1 if you want individual plots

% initial friction coefficient
mu_i_arr = [0.66 0.67 0.68 0.69 0.69 0.69 0.69 0.69];
% static (max) friction coefficient
mu_s_arr = 0.70*ones(1,np);
% transition friction coefficient
mu_t_arr = [0.65 0.675 0.68 0.687 0.68 0.69999999999 0.676 0.6765];
% dynamic friction coefficient
mu_d_arr = 0.60*ones(1,np); 
% strengthening distance
D_s_arr  = [0.125 0.11764705882352941176470588235294 0.11111111111111111111111111111111 0.10526315789473684210526315789474 0.10526315789473684210526315789474 0.10526315789473684210526315789474 0.10526315789473684210526315789474 0.10526315789473684210526315789474];
% weakening 1 distance
D_t_arr  = [0.5 0.5 0.51111111111111111111111111111111 0.4834759358288770053475935828877 0.51111111111111111111111111111111 0.39000000005050022844596473258021 0.51111111111111111111111111111111 0.51111111111111111111111111111111];                        
% weakening 2 distance
D_w_arr  = [1.3 0.7 0.6 0.57 2.044444444444455 0.62 1.618518518518526 1.663829787234045];

%% Set of three 
% np       = 3;                        % number of failure laws (plots)                       
% io       = [0 0 0];                  % set to 1 if you want individual plots
% 
% % initial friction coefficient
% mu_i_arr = [0.66 0.675 0.69];
% % static (max) friction coefficient
% mu_s_arr = 0.70*ones(1,np);
% % transition friction coefficient
% mu_t_arr = [0.65 0.675 0.687];
% % dynamic friction coefficient
% mu_d_arr = 0.60*ones(1,np); 
% % strengthening distance
% D_s_arr  = [0.125 0.11428571428571428571428571428571 0.10526315789473684210526315789474];
% % weakening 1 distance
% D_t_arr  = [0.5 0.5 0.4834759358288770053475935828877];                        
% % weakening 2 distance
% D_w_arr  = [1.3 0.7 0.57];              

%% Gentler + horizontal set
% np       = 4;                        % number of failure laws (plots)                       
% io       = [0 0 0 0];                % set to 1 if you want individual plots
% 
% % initial friction coefficient
% mu_i_arr = [0.66 0.67 0.68 0.69];
% % static (max) friction coefficient
% mu_s_arr = 0.70*ones(1,np);
% % transition friction coefficient
% mu_t_arr = [0.65 0.6872 0.6931 0.69999999999];
% % dynamic friction coefficient
% mu_d_arr = 0.60*ones(1,np); 
% % strengthening distance
% D_s_arr  = [0.125 0.11764705882352941176470588235294 0.11111111111111111111111111111111 0.10526315789473684210526315789474];
% % weakening 1 distance
% D_t_arr  = [0.5 0.48235042735042735042735042735043 0.44055411703780424650440186431901 0.39000000005050022844596473258021];
% % weakening 2 distance
% D_w_arr  = [1.3 0.57 0.59 0.62];

%% Last one - fully stable
% np       = 4;                        % number of failure laws (plots)                       
% io       = [0 0 0 0];                % set to 1 if you want individual plots
% 
% % initial friction coefficient
% mu_i_arr = [0.66 0.67 0.68 0.69];
% % static (max) friction coefficient
% mu_s_arr = 0.70*ones(1,np);
% % transition friction coefficient
% mu_t_arr = [0.65 0.6872 0.6931 0.68];
% % dynamic friction coefficient
% mu_d_arr = 0.60*ones(1,np); 
% % strengthening distance
% D_s_arr  = [0.125 0.11764705882352941176470588235294 0.11111111111111111111111111111111 0.10526315789473684210526315789474];
% % weakening 1 distance
% D_t_arr  = [0.5 0.48235042735042735042735042735043 0.44055411703780424650440186431901 0.51111111111111111111111111111111];
% % weakening 2 distance
% D_w_arr  = [1.3 0.57 0.59 2.044444444444455];

%% 6 laws
% np       = 6;                            % number of failure laws (plots)                       
% io       = [0 0 0 0 0 0];                % set to 1 if you want individual plots
% 
% % initial friction coefficient
% mu_i_arr = [0.64 0.65 0.66 0.67 0.68 0.69];
% % static (max) friction coefficient
% mu_s_arr = 0.70*ones(1,np);
% % transition friction coefficient
% mu_t_arr = [0.65 0.675 0.68 0.687 0.6931 0.69999999999];
% % dynamic friction coefficient
% mu_d_arr = 0.60*ones(1,np); 
% % strengthening distance
% D_s_arr  = [0.14285714285714285714285714285714 0.13333333333333333333333333333333 0.125 0.11764705882352941176470588235294 0.11111111111111111111111111111111 0.10526315789473684210526315789474];
% % weakening 1 distance
% D_t_arr  = [0.5 0.5 0.51111111111111111111111111111111 0.4834759358288770053475935828877 0.44055411703780424650440186431901 0.39000000005050022844596473258021];                        
% % weakening 2 distance
% D_w_arr  = [1.3 0.7 0.6 0.57 0.59 0.62];         

%% Plotting and supplementary variables
lw       = 0.75;                      % line width
fs       = 11;                        % font size
mpp      = 0.01*Sn;                   % min peak prominence for peak detection
v_pc     = 10*V_0;                    % preseismic to coseismic slip rate threshold
k        = 3;                         % number of cycles

Nt       = 10000;
Nta      = 200;

% Colors
b_col    = [0.231372549019608 0.298039215686275 0.752941176470588];
r_col    = [0.705882352941177 0.015686274509804 0.149019607843137];
g_col    = [0 .7 0];
k_col    = [0 0 0];
lb_col   = [0, 0.4470, 0.7410];
o_col    = [0.8500, 0.3250, 0.0980];
y_col    = [0.9290, 0.6940, 0.1250];
p_col    = [0.4940, 0.1840, 0.5560];
dr_col   = [0.6350, 0.0780, 0.1840];

%% Nondimentional variables for plotting
nd_t     = sqrt(M/K_s);                        % time
nd_u     = Sn*(mu_s_arr(1)-mu_d_arr(1))/K_s;   % slip
nd_v     = nd_u/nd_t;                          % slip rate
nd_a     = nd_v/nd_t;                          % acceleration

%% Reference solution for Coulomb friction (static -> dynamic, Dc=0)
% From Madariaga (simplified)
w0       = sqrt(K_s/M);
t_sd     = linspace(0, pi/w0, Nta);
u_sd_1   = Sn*(mu_s_arr(1)-mu_d_arr(1))/K_s*(1-cos(w0*t_sd));
v_sd_1   = Sn*(mu_s_arr(1)-mu_d_arr(1))/K_s*w0*sin(w0*t_sd);
a_sd_1   = Sn*(mu_s_arr(1)-mu_d_arr(1))/K_s*w0^2*cos(w0*t_sd);

% Full
A_c      = K_s/M;
B_c      = K_s*V_0/M;
C_c      = Sn*(mu_s_arr(1)-mu_d_arr(1))/M;
t_sd     = linspace(0, pi/sqrt(A_c), Nta);
u_sd     = C_c/A_c-C_c*cos(sqrt(A_c)*t_sd)/A_c-B_c*sin(sqrt(A_c)*t_sd)/(A_c^(3/2))+B_c/A_c*t_sd;
v_sd     = C_c*sin(sqrt(A_c)*t_sd)/sqrt(A_c)-B_c*cos(sqrt(A_c)*t_sd)/A_c+B_c/A_c;

% Preallocate
uSol = sym(0);

v_sp_a   = NaN(Nt, np+1);
t_sp_a   = NaN(Nt, np+1);

% Limits
max_uli  = zeros(1, np);
max_Sli  = zeros(1, np);
min_Sli  = zeros(1, np);
max_acli = zeros(1, np);
min_acli = zeros(1, np);
        
max_upr  = zeros(1, np);
max_Spr  = zeros(1, np);
min_Spr  = zeros(1, np);
max_acpr = zeros(1, np);
min_acpr = zeros(1, np);
    
max_vco  = zeros(1, np);
max_Sco  = zeros(1, np);
min_Sco  = zeros(1, np);
max_acco = zeros(1, np);
min_acco = zeros(1, np);
        
max_uf   = zeros(1, np);
max_u1c2 = zeros(1, np);

% Figure axis
f1_ax    = gobjects(1, np);
f11_ax   = gobjects(1, np);
f12_ax   = gobjects(1, np);
f2_ax    = gobjects(1, np);
f3_ax    = gobjects(1, np);
f35_ax   = gobjects(1, np);
f36_ax   = gobjects(1, np);
f31_ax   = gobjects(1, np);
f32_ax   = gobjects(1, np);
f33_ax   = gobjects(1, np);
p        = gobjects(1, 4);
pl       = gobjects(1, np+2);

%% Plot energy curves
figure(100);
cfl      = [b_col; r_col; g_col; k_col; lb_col; o_col; y_col; p_col; dr_col];                % colors
for fl = 1:np
    sl_v = cumsum([0 D_s_arr(fl) D_t_arr(fl) D_w_arr(fl) 3*D_w_arr(fl)]);
    mu_v = [mu_i_arr(fl) mu_s_arr(fl) mu_t_arr(fl) mu_d_arr(fl) mu_d_arr(fl)];
    pl(fl) = plot(sl_v/nd_u, (mu_v-mu_d_arr(1))/(mu_s_arr(1)-mu_d_arr(1)), 'Color', cfl(fl,:), 'LineWidth', lw, 'DisplayName',lower(char(fl+64)));
    hold on;
end
pl(np+1) = plot(linspace(D_s_arr(2),1.1*nd_u,10)/nd_u, -K_s*linspace(0,1.1*nd_u-D_s_arr(2),10)/Sn/(mu_s_arr(1)-mu_d_arr(1))+1, 'k--', 'LineWidth', lw, 'DisplayName', 'Loading');
set(gca,'FontSize',fs-2);
ylabel('$\tau_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
xlabel('$D_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
legend('show', 'Interpreter', 'latex', 'FontSize', fs-1);
grid on;
xlim([0 1.1]);
ylim([-0.02 1.02]);

% Loop through np failure laws
for r = 1:np
    D_s      = D_s_arr(r);
    D_t      = D_t_arr(r);
    D_w      = D_w_arr(r);
    
    mu_i     = mu_i_arr(r);
    mu_s     = mu_s_arr(r);
    mu_t     = mu_t_arr(r);
    mu_d     = mu_d_arr(r);
    
    %% Slopes and energies
    K_str    = Sn*(mu_i-mu_s)/D_s;
    K_w1     = Sn*(mu_s-mu_t)/D_t;
    K_w2     = Sn*(mu_t-mu_d)/D_w;
    en_w     = Sn*((mu_s-mu_t)*D_t/2+(mu_t-mu_d)*D_w/2+(mu_t-mu_d)*D_t);
    en_s     = Sn*((mu_s-mu_i)*D_s/2+(mu_i-mu_d)*D_s);
    
    disp(['Case ' lower(char(r+64))]); 
    fprintf('K_str = %.0f K_w1 = %.0f K_w2 = %.0f \n', K_str, K_w1, K_w2);
    fprintf('e_w = %.0f e_s = %.0f \n', en_w, en_s);

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
    vsp      = [vel_s vel_w1(2:end) vel_w2(2:end) vel_fl(2:end)];
    tsp      = [time_s time_w1(2:end) time_w2(2:end) time_fl(2:end)];
    
    v_sp_a(1:length(vsp),r)   = vsp;
    t_sp_a(1:length(vsp),r)   = tsp;
    
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
    % Save limits for the future
    max_uli(r)  = max(disp_li);
    max_Sli(r)  = (max([ku_li tau_li])/Sn+mu_i-mu_d)/(mu_s-mu_d);
    min_Sli(r)  = (min([ku_li tau_li])/Sn+mu_i-mu_d)/(mu_s-mu_d);
    max_acli(r) = max(ac_li);
    min_acli(r) = min(ac_li);
        
    max_upr(r)  = max(disp_pr);
    max_Spr(r)  = (max([ku_pr tau_pr])/Sn+mu_i-mu_d)/(mu_s-mu_d);
    min_Spr(r)  = (min([ku_pr tau_pr])/Sn+mu_i-mu_d)/(mu_s-mu_d);
    max_acpr(r) = max(ac_pr);
    min_acpr(r) = min(ac_pr);
    
    max_vco(r)  = max(vel_co);
    max_Sco(r)  = (max([ku_co tau_co])/Sn+mu_i-mu_d)/(mu_s-mu_d);
    min_Sco(r)  = (min([ku_co tau_co])/Sn+mu_i-mu_d)/(mu_s-mu_d);
    max_acco(r) = max(ac_co);
    min_acco(r) = min(ac_co);
        
    max_uf(r)   = max(disp_f);
    max_u1c2(r) = max(disp_1c2);
      
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

    figure(1);
    f1_ax(r)   = subplot(1,np,r);
    % Fill energies
    xb         = (disp_f(locs_t(pn):idx_s)/nd_u)';
    y1b        = ((tau_f(locs_t(pn):idx_s)-tau_d)/Sn/(mu_s-mu_d))';
    y2b        = (tau_d*ones(size(y1b))-tau_d)/Sn/(mu_s-mu_d);
    fill([xb; flipud(xb)], [y2b; flipud(y1b)], 'b', 'EdgeColor', [.5 .5 .5], 'LineWidth', 0.1, 'FaceAlpha', 0.35);
    hold on;
    xr         = [u_int disp_f(idx_w+1:locs_t(pn))]'/nd_u;
    y1r        = [tau_d-tau_d tau_f(idx_w+1:locs_t(pn))-tau_d]'/Sn/(mu_s-mu_d);
    y2r        = (tau_d*ones(size(y1r))-tau_d)/Sn/(mu_s-mu_d);
    fill([xr; flipud(xr)], [y2r; flipud(y1r)], 'r', 'EdgeColor', [.5 .5 .5], 'LineWidth', 0.1, 'FaceAlpha', 0.35);
    xg         = (disp_f(locs_t(pn):idx_s_ku)/nd_u)';
    y1g        = ((ku_f(locs_t(pn):idx_s_ku)-tau_d)/Sn/(mu_s-mu_d))';
    y2g        = ((tau_f(locs_t(pn):idx_s_ku)-tau_d)/Sn/(mu_s-mu_d))';
    fill([xg; flipud(xg)], [y2g; flipud(y1g)], 'k', 'EdgeColor', [.4 .5 .5], 'LineWidth', 0.1, 'FaceAlpha', 0.2);
       % Plot Stress vs. slip curves
    j(1) = plot(disp_f/nd_u, (tau_f/Sn+mu_i-mu_d)/(mu_s-mu_d), 'k', 'LineWidth', lw);
    j(2) = plot(disp_f/nd_u, (ku_f/Sn+mu_i-mu_d)/(mu_s-mu_d), 'k--', 'LineWidth', lw);
    legend([j(1) j(2)],'Friction','Loading', 'Interpreter', 'latex', 'FontSize', fs-1, 'location', 'southwest');
    set(gca,'FontSize',fs-2);
    title(lower(char(r+64)), 'Interpreter', 'latex', 'FontSize', fs+1, 'FontWeight', 'bold');
    if r == 1
        ylabel('$\tau_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
    else
        set(gca,'yticklabel',[]);
    end
    xlabel('$D_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
    grid on;
    axis tight;
    hold off;
    % Link axis limits
    if r == np
        hold on;
        ylim([min(min_Sco) max(max_Sli)]*1.02);
        linkaxes(f1_ax, 'xy');
        hold off;
    end

    %% Stress vs. slip (zoom)
    figure(11);
    f11_ax(r)   = subplot(1,np,r);
    % Fill energies
    fill([xb; flipud(xb)], [y2b; flipud(y1b)], 'b', 'EdgeColor', [.5 .5 .5], 'LineWidth', 0.1, 'FaceAlpha', 0.35);
    hold on;
    fill([xr; flipud(xr)], [y2r; flipud(y1r)], 'r', 'EdgeColor', [.5 .5 .5], 'LineWidth', 0.1, 'FaceAlpha', 0.35);
    fill([xg; flipud(xg)], [y2g; flipud(y1g)], 'k', 'EdgeColor', [.5 .5 .5], 'LineWidth', 0.1, 'FaceAlpha', 0.2);
    
    % Plot Stress vs. slip curves
    j(1) = plot(disp_f/nd_u, (tau_f/Sn+mu_i-mu_d)/(mu_s-mu_d), 'k', 'LineWidth', lw);
    j(2) = plot(disp_f/nd_u, (ku_f/Sn+mu_i-mu_d)/(mu_s-mu_d), 'k--', 'LineWidth', lw);
    legend([j(1) j(2)],'Friction','Loading', 'Interpreter', 'latex', 'FontSize', fs-4, 'location', 'northeast');
    set(gca,'FontSize',fs-4);
    title(lower(char(r+64)), 'Interpreter', 'latex', 'FontSize', fs+1, 'FontWeight', 'bold');
    if r == 1
        ylabel('$\tau_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
    else
        set(gca,'yticklabel',[]);
    end
    xlabel('$D_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
    xlim([min(xr)*0.95 max(xb)*1.25]);
    ylim([-0.03 max_Sli(r)]*1.02);
    grid on;
    hold off;
    % Link axis limits
    if r == np
        hold on;
        linkaxes(f11_ax, 'xy');
        hold off;
    end
    
    %% Stress vs. slip (single cycle)
    figure(12);
    f12_ax(r)   = subplot(1,np,r);
    % Fill energies
    fill([xb-min(xr); flipud(xb-min(xr))], [y2b; flipud(y1b)], 'b', 'EdgeColor', [.5 .5 .5], 'LineWidth', 0.1, 'FaceAlpha', 0.35);
    hold on;
    fill([xr-min(xr); flipud(xr-min(xr))], [y2r; flipud(y1r)], 'r', 'EdgeColor', [.5 .5 .5], 'LineWidth', 0.1, 'FaceAlpha', 0.35);
    fill([xg-min(xr); flipud(xg-min(xr))], [y2g; flipud(y1g)], 'k', 'EdgeColor', [.5 .5 .5], 'LineWidth', 0.1, 'FaceAlpha', 0.2);
    
    % Plot Stress vs. slip curves
    j(1) = plot(disp_1c2/nd_u, (tau_1c2/Sn+mu_i-mu_d)/(mu_s-mu_d), 'k', 'LineWidth', lw);
    j(2) = plot(disp_1c2/nd_u, (ku_1c2/Sn+mu_i-mu_d)/(mu_s-mu_d), 'k--', 'LineWidth', lw);
    legend([j(1) j(2)],'Friction','Loading', 'Interpreter', 'latex', 'FontSize', fs-1, 'location', 'northeast');
    set(gca,'FontSize',fs-2);
    title(lower(char(r+64)), 'Interpreter', 'latex', 'FontSize', fs+1, 'FontWeight', 'bold');
    if r == 1
        ylabel('$\tau_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
    else
        set(gca,'yticklabel',[]);
    end
    xlabel('$D_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
    xlim([-0.05 max(max_u1c2/nd_u)*1.02]);
    ylim([min(min_Sco) max(max_Sli)]*1.02);
    grid on;
    hold off;
    % Link axis limits
    if r == np
        hold on;
        linkaxes(f12_ax, 'xy');
        hold off;
    end
    
    %% Phase trajectories (several cycles) fails if int_vel == 0
    figure(2);
    f2_ax(r) = subplot(1,np,r);
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
    if idx_pc ~= idx_cei
         legend(p(1:4), {'Early interseismic', 'Late interseismic', 'Preseismic', 'Coseismic'}, 'Interpreter', 'latex', 'FontSize', fs-1, 'location', 'southeast');
    else
         legend(p(1:3), {'Early interseismic', 'Late interseismic', '"Preseismic"'}, 'Interpreter', 'latex', 'FontSize', fs-1, 'location', 'southeast');
    end
    set(gca,'FontSize',fs-2);
    xlabel('$ln(V_{nd})$', 'Interpreter', 'latex', 'FontSize', fs);
    axis tight;
    title(lower(char(r+64)), 'Interpreter', 'latex', 'FontSize', fs+1, 'FontWeight', 'bold');
    if r == 1
        ylabel('$\tau_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
    else
        set(gca,'yticklabel',[]);
    end
    % Link axis limits
    if r == np
        hold on
        xlim([log(int_vel*0.6/nd_v) log(max(max_vco*1.8/nd_v))]);
        ylim([min(min_Sco) max(max_Sli)]*1.02);
        linkaxes(f2_ax, 'xy');
        hold off;
    end

    %% Vel and slip in one plot (several cycles)
    fig3 = figure(3);
    set(fig3, 'defaultAxesColorOrder',[r_col; b_col]);
    f3_ax(r) = subplot(1,np,r);
    
    yyaxis right;
    p(1) = plot(time_f/nd_t, disp_f/nd_u, 'LineWidth', lw);
    hold on;
    p(2) = plot(time_f/nd_t, (V_0*(time_f-max(time_co) - max(time_pr) - max(time_li)) + max(disp_co) + max(disp_pr) + max(disp_li))/nd_u, 'LineWidth', lw);
    set(gca,'FontSize',fs-4);
    if r == np
        ylabel('$D_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
    else
        set(gca,'yticklabel',[]);
    end
    
    yyaxis left;
    plot(time_f/nd_t, vel_f/nd_v, 'LineWidth', lw);
    legend([p(1) p(2)],'Slider','Load point', 'Interpreter', 'latex', 'FontSize', fs-4, 'location', 'southwest');
    if r == 1
        ylabel('$V_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
    else
        set(gca,'yticklabel',[]);
    end
    xlabel('$t_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
    title(lower(char(r+64)), 'Interpreter', 'latex', 'FontSize', fs+1, 'FontWeight', 'bold');
    axis tight;
    grid on;
    set(gca,'GridColor',[0 0 0]);
    % Link axis limits
    if r == np
        yyaxis right;
        Y = get(f3_ax,'YAxis');
        Y = [Y{:}];
        linkprop(Y(2,:),'Limits');
        ylim([0 max(max_uf/nd_u)*1.03]);

        yyaxis left;
        linkaxes(f3_ax, 'xy');
        ylim([0 max(max_vco/nd_v)*1.02]);
%         ylim([log(int_vel*0.6/nd_v) log(max(max_vco*1.8/nd_v))]);
        hold off;
    end
    
    %% Vel and slip in one plot (single cycle)
    fig35 = figure(35);
    set(fig35, 'defaultAxesColorOrder',[r_col; b_col]);
    f35_ax(r) = subplot(1,np,r);
    
    yyaxis right;
    p(1) = plot(time_1c2/nd_t, disp_1c2/nd_u, 'LineWidth', lw);
    hold on;
    p(2) = plot(time_1c2/nd_t, V_0*time_1c2/nd_u, 'LineWidth', lw);
    set(gca,'FontSize',fs-2);
    if r == np
        ylabel('$D_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
    else
        set(gca,'yticklabel',[]);
    end
    
    yyaxis left;
    plot(time_1c2/nd_t, vel_1c2/nd_v, 'LineWidth', lw);
    legend([p(1) p(2)],'Slider','Load point', 'Interpreter', 'latex', 'FontSize', fs-1, 'location', 'southwest');
    if r == 1
        ylabel('$V_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
    else
        set(gca,'yticklabel',[]);
    end
    xlabel('$t_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
    title(lower(char(r+64)), 'Interpreter', 'latex', 'FontSize', fs+1, 'FontWeight', 'bold');
    axis tight;
    grid on;
    set(gca,'GridColor',[0 0 0]);
    % Link axis limits
    if r == np
        yyaxis right;
        Y = get(f35_ax,'YAxis');
        Y = [Y{:}];
        linkprop(Y(2,:),'Limits');
        ylim([0 max(max_u1c2/nd_u)*1.03]);

        yyaxis left;
        linkaxes(f35_ax, 'xy');
        ylim([0 max(max_vco/nd_v)*1.02]);
        hold off;
    end

    %% Vel and slip in one plot (log zoom)
    fig36 = figure(36);
    set(fig36, 'defaultAxesColorOrder',[r_col; 0 0 0]);
    f36_ax(r) = subplot(1,np,r);
    
    yyaxis right;
    set(gca,'yticklabel',[], 'FontSize',fs-2);
    
    yyaxis left;
    semilogy(time_1c2/nd_t, vel_1c2/nd_v, 'LineWidth', lw);
    if r == 1
        ylabel('$ln(V_{nd})$', 'Interpreter', 'latex', 'FontSize', fs);
    else
        set(gca,'yticklabel',[]);
    end
    xlabel('$t_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
    title(lower(char(r+64)), 'Interpreter', 'latex', 'FontSize', fs+1, 'FontWeight', 'bold');
    axis tight;
    grid on;
    set(gca,'GridColor',[0 0 0])
    % Link axis limits
    if r == np
        yyaxis right;
        Y = get(f36_ax,'YAxis');
        Y = [Y{:}];
        linkprop(Y(2,:),'Limits');
        ylim([0 max(max_u1c2/nd_u)*1.03]);

        yyaxis left;
        linkaxes(f36_ax, 'xy');
        ylim([0 max(max_vco/nd_v)*1.02]);
        xlim([280 520]);
        hold off;
    end
    
    
    %% Late interseismic phase
    fig31 = figure(31);
    set(fig31, 'defaultAxesColorOrder',[r_col; b_col]);

    % Stress and acceleration vs. time
    f31_ax(r) = subplot(2,np,r);
    yyaxis right;
    p(2) = plot(time_li/nd_t, (tau_li/Sn+mu_i-mu_d)/(mu_s-mu_d), 'LineWidth', lw);
    hold on;
    p(3) = plot(time_li/nd_t, (ku_li/Sn+mu_i-mu_d)/(mu_s-mu_d), 'LineWidth', lw);
    set(gca,'FontSize',fs-2);
    if r == np
        ylabel('$\tau_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
    else
        set(gca,'yticklabel',[]);
    end
    axis tight;
    yyaxis left;
    p(1) = plot(time_li/nd_t, ac_li/nd_a, 'LineWidth', lw);
    if r == 1
        ylabel('$a_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
    else
        set(gca,'yticklabel',[]);
    end
    set(gca,'xticklabel',[]);
    title(lower(char(r+64)), 'Interpreter', 'latex', 'FontSize', fs+1, 'FontWeight', 'bold');
    axis tight;
    grid on;
    set(gca,'GridColor',[0 0 0])
    hold on;

    % Slip and slip rate vs. time
    f31_ax(np+r) = subplot(2,np,np+r);
    yyaxis right;
    plot(time_li/nd_t, disp_li/nd_u, 'LineWidth', lw);
    set(gca,'FontSize',fs-2);
    if r == np
        ylabel('$D_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
    else
        set(gca,'yticklabel',[]);
    end
    xlabel('$t_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
    yyaxis left;
    plot(time_li/nd_t, vel_li/nd_v, 'LineWidth', lw);
    if r == 1
        ylabel('$V_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
    else
        set(gca,'yticklabel',[]);
    end
    axis tight;
    grid on;
    set(gca,'GridColor',[0 0 0])
    hold on;
    legend([p(2) p(3)],'Friction','Loading', 'Interpreter', 'latex', 'FontSize', fs-1, 'location', 'southwest');
    % Link axis limits
    if r == np
        subplot(2,np,1)
        yyaxis right;
        ylim([min(min_Sli)*0.98 max(max_Sli)*1.02]);
        Y = get(f31_ax(1:np),'YAxis');
        Y = [Y{:}];
        linkprop(Y(2,:),'Limits');
        hold off;
        yyaxis left;
        linkaxes(f31_ax(1:np), 'xy');
        ylim([min(min_acli) max(max_acli)]*1.1/nd_a);
        hold off;

        subplot(2,np,np+1)
        yyaxis right;
        ylim([0 max(max_uli)*1.03]/nd_u);
        Y = get(f31_ax(np+1:2*np),'YAxis');
        Y = [Y{:}];
        linkprop(Y(2,:),'Limits');
        hold off;
        yyaxis left;
        linkaxes(f31_ax(np+1:2*np), 'xy');
        ylim([0 V_0*1.03]/nd_v);
        hold off;
    end
    
    %% Preseismic phase
    fig32 = figure(32);
    set(fig32, 'defaultAxesColorOrder',[r_col; b_col]);
    
    % Stress and acceleration vs. time
    f32_ax(r) = subplot(2,np,r);
    yyaxis right;
    p(2) = plot(time_pr/nd_t, (tau_pr/Sn+mu_i-mu_d)/(mu_s-mu_d), 'LineWidth', lw);
    hold on;
    p(3) = plot(time_pr/nd_t, (ku_pr/Sn+mu_i-mu_d)/(mu_s-mu_d), 'LineWidth', lw);
    set(gca,'FontSize',fs-2);
    if r == np
        ylabel('$\tau_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
    else
        set(gca,'yticklabel',[]);
    end
    axis tight;
    yyaxis left;
    p(1) = plot(time_pr/nd_t, ac_pr/nd_a, 'LineWidth', lw);
    title(lower(char(r+64)), 'Interpreter', 'latex', 'FontSize', fs+1, 'FontWeight', 'bold');
    if r == 1
        ylabel('$a_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
    else
        set(gca,'yticklabel',[]);
    end
    set(gca,'xticklabel',[]);
    axis tight;
    grid on;
    set(gca,'GridColor',[0 0 0]);

    % Slip and slip rate vs. time
    f32_ax(np+r) = subplot(2,np,np+r);
    yyaxis right;
    plot(time_pr/nd_t, disp_pr/nd_u, 'LineWidth', lw);
    set(gca,'FontSize',fs-2);
    if r == np
        ylabel('$D_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
    else
        set(gca,'yticklabel',[]);
    end
    xlabel('$t_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
    yyaxis left;
    plot(time_pr/nd_t, vel_pr/nd_v, 'LineWidth', lw);
    if r == 1
        ylabel('$V_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
    else
        set(gca,'yticklabel',[]);
    end
    axis tight;
    grid on;
    set(gca,'GridColor',[0 0 0]);
    legend([p(2) p(3)],'Friction','Loading', 'Interpreter', 'latex', 'FontSize', fs-1, 'location', 'southwest');
    % Link axis limits
    if r == np
        subplot(2,np,1)
        yyaxis right;
        ylim([min(min_Spr)*0.98 max(max_Spr)*1.02]);
        Y = get(f32_ax(1:np),'YAxis');
        Y = [Y{:}];
        linkprop(Y(2,:),'Limits');
        hold off;
        yyaxis left;
        linkaxes(f32_ax(1:np), 'xy');
        ylim([min(min_acpr) max(max_acpr)]*1.1/nd_a);
        hold off;

        subplot(2,np,np+1)
        yyaxis right;
        ylim([0 max(max_upr)*1.03]/nd_u);
        Y = get(f32_ax(np+1:2*np),'YAxis');
        Y = [Y{:}];
        linkprop(Y(2,:),'Limits');
        hold off;
        yyaxis left;
        linkaxes(f32_ax(np+1:2*np), 'xy');
        ylim([0 v_pc*1.03]/nd_v);
        hold off;
    end

    %% Coseismic phase
    % Find index of the beginning of the flat segment
    [~,idx_fl] = find((tau_co/Sn+mu_i-mu_d)/(mu_s-mu_d)<1e-13, 1, 'first');

    fig33 = figure(33);
    set(fig33, 'defaultAxesColorOrder',[r_col; b_col]);

    % Stress and acceleration vs. time
    f33_ax(r) = subplot(2,np,r);
    yyaxis right;
    p(2) = plot(time_co/nd_t, (tau_co/Sn+mu_i-mu_d)/(mu_s-mu_d), 'LineWidth', lw);
    hold on;
    p(3) = plot(time_co/nd_t, (ku_co/Sn+mu_i-mu_d)/(mu_s-mu_d), 'LineWidth', lw);
    xline((time_co(idx_fl))/nd_t, 'k-.', 'LineWidth', lw-0.5);
    set(gca,'FontSize',fs-2);
    if r == np
        ylabel('$\tau_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
    else
        set(gca,'yticklabel',[]);
    end
    axis tight;
    yyaxis left;
    p(1) = plot(time_co/nd_t, ac_co/nd_a, 'LineWidth', lw);
    hold on;
    if r == 1
        ylabel('$a_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
    else
        set(gca,'yticklabel',[]);
    end
    set(gca,'xticklabel',[]);
    title(lower(char(r+64)), 'Interpreter', 'latex', 'FontSize', fs+1, 'FontWeight', 'bold');
    axis tight;
    grid on;
    set(gca,'GridColor',[0 0 0]);

    % Slip and slip rate vs. time
    f33_ax(np+r) = subplot(2,np,np+r);
    yyaxis right;
    plot(time_co/nd_t, disp_co/nd_u, 'LineWidth', lw);
    hold on;
    [~,I_m] = max(vel_co);
    if time_co(I_m) == 0
        tc_Im = max(t_sd)/2;
    else
        tc_Im = time_co(I_m);
    end
    plot((t_sd-max(t_sd)/2+tc_Im)/nd_t, u_sd_1/nd_u, ':', 'LineWidth', lw);
    set(gca,'FontSize',fs-2);
    if r == np
        ylabel('$D_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
    else
        set(gca,'yticklabel',[]);
    end
    xlabel('$t_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
    yyaxis left;
    plot(time_co/nd_t, vel_co/nd_v, 'LineWidth', lw);
    hold on;
    plot((t_sd-max(t_sd)/2+tc_Im)/nd_t, v_sd_1/nd_v, ':', 'LineWidth', lw);
    xline((time_co(idx_fl))/nd_t, 'k-.', 'LineWidth', lw-0.5);
    if r == 1
        ylabel('$V_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
    else
        set(gca,'yticklabel',[]);
    end
    axis tight;
    grid on;
    set(gca,'GridColor',[0 0 0]);
    legend([p(2) p(3)],'Friction','Loading', 'Interpreter', 'latex', 'FontSize', fs-1, 'location', 'southwest');
    % Link axis limits
    if r == np
        subplot(2,np,1)
        yyaxis right;
        ylim([min(min_Sco) max(max_Sco)]*1.02);
        Y = get(f33_ax(1:np),'YAxis');
        Y = [Y{:}];
        linkprop(Y(2,:),'Limits');
        hold off;
        yyaxis left;
        linkaxes(f33_ax(1:np), 'xy');
        ylim([min(min_acco) max(max_acco)]*1.1/nd_a);
        hold off;

        subplot(2,np,np+1)
        yyaxis right;
        ylim([0 max(u_sd)*1.03]/nd_u);
        Y = get(f33_ax(np+1:2*np),'YAxis');
        Y = [Y{:}];
        linkprop(Y(2,:),'Limits');
        hold off;
        yyaxis left;
        linkaxes(f33_ax(np+1:2*np), 'xy');
        ylim([0 max(v_sd)*1.03]/nd_v);
        hold off;
    end

    %% Individual plots
    if io(r) == 1
    
        %% Full cycle (4 plots)
        figure(300+r);
    
        % Stress vs. slip
        f100_ax(1) = subplot(2,2,1);
        fill([xb; flipud(xb)], [y2b; flipud(y1b)], 'b', 'EdgeColor', [.5 .5 .5], 'LineWidth', 0.1, 'FaceAlpha', 0.35);
        hold on;
        fill([xr; flipud(xr)], [y2r; flipud(y1r)], 'r', 'EdgeColor', [.5 .5 .5], 'LineWidth', 0.1, 'FaceAlpha', 0.35);
        fill([xg; flipud(xg)], [y2g; flipud(y1g)], 'k', 'EdgeColor', [.4 .5 .5], 'LineWidth', 0.1, 'FaceAlpha', 0.2);
        j(1) = plot(disp_f/nd_u, (tau_f/Sn+mu_i-mu_d)/(mu_s-mu_d), 'k', 'LineWidth', lw);
        j(2) = plot(disp_f/nd_u, (ku_f/Sn+mu_i-mu_d)/(mu_s-mu_d), 'k--', 'LineWidth', lw);
        legend([j(1) j(2)],'Friction','Loading', 'Interpreter', 'latex', 'FontSize', fs-1, 'location', 'southwest');
        grid on;
        axis tight;
        hold off;
        set(gca,'FontSize',fs-2);
        ylabel('$\tau_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
        xlabel('$D_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
        ylim([min_Sco(r) max_Sli(r)]*1.02);

        % Stress vs. slip (zoom)
        f100_ax(2) = subplot(2,2,2);
        fill([xb; flipud(xb)], [y2b; flipud(y1b)], 'b', 'EdgeColor', [.5 .5 .5], 'LineWidth', 0.1, 'FaceAlpha', 0.35);
        hold on;
        fill([xr; flipud(xr)], [y2r; flipud(y1r)], 'r', 'EdgeColor', [.5 .5 .5], 'LineWidth', 0.1, 'FaceAlpha', 0.35);
        fill([xg; flipud(xg)], [y2g; flipud(y1g)], 'k', 'EdgeColor', [.4 .5 .5], 'LineWidth', 0.1, 'FaceAlpha', 0.2);
        j(1) = plot(disp_f/nd_u, (tau_f/Sn+mu_i-mu_d)/(mu_s-mu_d), 'k', 'LineWidth', lw);
        j(2) = plot(disp_f/nd_u, (ku_f/Sn+mu_i-mu_d)/(mu_s-mu_d), 'k--', 'LineWidth', lw);
        grid on;
        set(gca,'FontSize',fs-2);
        ylabel('$\tau_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
        xlabel('$D_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
        ylim([-0.03 max_Sli(r)]*1.02);
        xlim([min(xr)*0.95 max(xb)*1.05]);
        hold off;

        % Phase diagrams
        f100_ax(3) = subplot(2,2,3);
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
        if idx_pc ~= idx_cei
            legend(p(1:4), {'Early interseismic', 'Late interseismic', 'Preseismic', 'Coseismic'}, 'Interpreter', 'latex', 'FontSize', fs-1, 'location', 'southeast');
        else
            legend(p(1:3), {'Early interseismic', 'Late interseismic', '"Preseismic"'}, 'Interpreter', 'latex', 'FontSize', fs-1, 'location', 'southeast');
        end
        set(gca,'FontSize',fs-2);
        xlabel('$ln(V_{nd})$', 'Interpreter', 'latex', 'FontSize', fs);
        ylabel('$\tau_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
        xlim([log(int_vel*0.6/nd_v) log(max([max_vco(r) v_pc])*1.8/nd_v)]);
        ylim([min_Sco(r) max_Sli(r)]*1.02);

        % Slip and slip rate (several cycles)
        f100_ax(4) = subplot(2,2,4);
        yyaxis right;
        p(1) = plot(time_f/nd_t, disp_f/nd_u, 'Color', b_col, 'LineWidth', lw);
        hold on;
        p(2) = plot(time_f/nd_t, (V_0*(time_f-max(time_co) - max(time_pr) - max(time_li)) + max(disp_co) + max(disp_pr) + max(disp_li))/nd_u, 'Color', b_col, 'LineWidth', lw);
        set(gca,'FontSize',fs-2);
        ylabel('$D_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
        axis tight;
        ylim([0 max_uf(r)/nd_u*1.03]);
        yyaxis left;
        plot(time_f/nd_t, vel_f/nd_v, 'Color', r_col, 'LineWidth', lw);
        legend([p(1) p(2)],'Slider','Load point', 'Interpreter', 'latex', 'FontSize', fs-1, 'location', 'southeast');
        ylim([0 max([max_vco(r) v_pc])/nd_v*1.02]);
        ylabel('$V_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
        xlabel('$t_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
        ax = gca;
        ax.YAxis(1).Color = r_col;
        ax.YAxis(2).Color = b_col;
        grid on;
        set(gca,'GridColor',[0 0 0]);
        hold off;
        
        %% Full cycle (3 plots)
        figure(400+r);
    
        % Stress vs. slip
        f400_ax(1) = subplot(2,2,1);
        fill([xb-min(xr); flipud(xb-min(xr))], [y2b; flipud(y1b)], 'b', 'EdgeColor', [.5 .5 .5], 'LineWidth', 0.1, 'FaceAlpha', 0.35);
        hold on;
        fill([xr-min(xr); flipud(xr-min(xr))], [y2r; flipud(y1r)], 'r', 'EdgeColor', [.5 .5 .5], 'LineWidth', 0.1, 'FaceAlpha', 0.35);
        fill([xg-min(xr); flipud(xg-min(xr))], [y2g; flipud(y1g)], 'k', 'EdgeColor', [.4 .5 .5], 'LineWidth', 0.1, 'FaceAlpha', 0.2);
        j(1) = plot(disp_1c2/nd_u, (tau_1c2/Sn+mu_i-mu_d)/(mu_s-mu_d), 'k', 'LineWidth', lw);
        j(2) = plot(disp_1c2/nd_u, (ku_1c2/Sn+mu_i-mu_d)/(mu_s-mu_d), 'k--', 'LineWidth', lw);
        legend([j(1) j(2)],'Friction','Loading', 'Interpreter', 'latex', 'FontSize', fs-1, 'location', 'southwest');
        grid on;
        axis tight;
        hold off;
        set(gca,'FontSize',fs-2);
        ylabel('$\tau_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
        xlabel('$D_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
        ylim([min_Sco(r) max_Sli(r)]*1.02);
        xlim([-0.05 max(max_u1c2/nd_u)*1.02]);

        % Phase diagrams
        f400_ax(3) = subplot(2,2,2);
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
        if idx_pc ~= idx_cei
                legend(p(1:4), {'Early interseismic', 'Late interseismic', 'Preseismic', 'Coseismic'}, 'Interpreter', 'latex', 'FontSize', fs-1, 'location', 'southeast');
        else
                legend(p(1:3), {'Early interseismic', 'Late interseismic', '"Preseismic"'}, 'Interpreter', 'latex', 'FontSize', fs-1, 'location', 'southeast');
        end
        set(gca,'FontSize',fs-2);
        xlabel('$ln(V_{nd})$', 'Interpreter', 'latex', 'FontSize', fs);
        ylabel('$\tau_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
        xlim([log(int_vel*0.6/nd_v) log(max([max_vco(r) v_pc])*1.8/nd_v)]);
        ylim([min_Sco(r) max_Sli(r)]*1.02);

        % Slip and slip rate (single cycle)
        f400_ax(4) = subplot(2,2,[3 4]);
        yyaxis right;
        p(1) = plot(time_1c2/nd_t, disp_1c2/nd_u, 'Color', b_col, 'LineWidth', lw);
        hold on;
        p(2) = plot(time_1c2/nd_t, V_0*time_1c2/nd_u, 'Color', b_col, 'LineWidth', lw);
        set(gca,'FontSize',fs-2);
        ylabel('$D_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
        axis tight;
        ylim([0 max_u1c2(r)/nd_u*1.03]);
        yyaxis left;
        plot(time_1c2/nd_t, vel_1c2/nd_v, 'Color', r_col, 'LineWidth', lw);
        legend([p(1) p(2)],'Slider','Load point', 'Interpreter', 'latex', 'FontSize', fs-1, 'location', 'northwest');
        ylim([0 max([max_vco(r) v_pc])/nd_v*1.02]);
        ylabel('$V_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
        xlabel('$t_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
        ax = gca;
        ax.YAxis(1).Color = r_col;
        ax.YAxis(2).Color = b_col;
        grid on;
        set(gca,'GridColor',[0 0 0]);
        hold off;
        
        %% Phases separately
        fig200 = figure(200+r);
        set(fig200, 'defaultAxesColorOrder',[r_col; b_col]);

        % Check if there's coseismic phase
        if idx_pc ~= idx_cei
                spn = 3;
        else
                spn = 2;
        end

        % Late interseismic
        % Stress and acceleration vs. time
        f200_ax(1) = subplot(2,spn,1);
        yyaxis right;
        p(2) = plot(time_li/nd_t, (tau_li/Sn+mu_i-mu_d)/(mu_s-mu_d), 'LineWidth', lw);
        hold on;
        p(3) = plot(time_li/nd_t, (ku_li/Sn+mu_i-mu_d)/(mu_s-mu_d), 'LineWidth', lw);
        axis tight;
        ylim([min_Sli(r)*0.98 max_Sli(r)*1.02]);
        yyaxis left;
        p(1) = plot(time_li/nd_t, ac_li/nd_a, 'LineWidth', lw);
        set(gca,'FontSize',fs-2);
        set(gca,'xticklabel',[]);
        title('Late\ interseismic', 'Interpreter', 'latex', 'FontSize', fs+1, 'FontWeight', 'bold');
        ylabel('$a_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
        axis tight;
        ylim([min_acli(r) max_acli(r)]*1.1/nd_a);
        grid on;
        set(gca,'GridColor',[0 0 0]);
        
        % Slip and slip rate vs. time
        f200_ax(4) = subplot(2,spn,spn+1);
        yyaxis right;
        plot(time_li/nd_t, disp_li/nd_u, 'LineWidth', lw);
        xlabel('$t_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
        ylim([0 max_uli(r)*1.03]/nd_u);
        yyaxis left;
        plot(time_li/nd_t, vel_li/nd_v, 'LineWidth', lw);
        set(gca,'FontSize',fs-2);
        ylabel('$V_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
        axis tight;
        grid on;
        set(gca,'GridColor',[0 0 0]);
        ylim([0 V_0*1.03]/nd_v);
        legend([p(2) p(3)],'Friction','Loading', 'Interpreter', 'latex', 'FontSize', fs-1, 'location', 'southwest');

        % Preseismic
        % Stress and acceleration vs. time
        f200_ax(2) = subplot(2,spn,2);
        yyaxis right;
        p(2) = plot(time_pr/nd_t, (tau_pr/Sn+mu_i-mu_d)/(mu_s-mu_d), 'LineWidth', lw);
        hold on;
        p(3) = plot(time_pr/nd_t, (ku_pr/Sn+mu_i-mu_d)/(mu_s-mu_d), 'LineWidth', lw);
        axis tight;
        ylim([min_Spr(r)*0.98 max_Spr(r)*1.02]);
        set(gca,'FontSize',fs-2);
        if idx_pc == idx_cei
            ylabel('$\tau_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
        end    
        yyaxis left;
        p(1) = plot(time_pr/nd_t, ac_pr/nd_a, 'LineWidth', lw);
        set(gca,'xticklabel',[]);
        if idx_pc ~= idx_cei
                title('Preseismic', 'Interpreter', 'latex', 'FontSize', fs+1, 'FontWeight', 'bold');
        else
                title('"Preseismic"', 'Interpreter', 'latex', 'FontSize', fs+1, 'FontWeight', 'bold');
        end
        axis tight;
        ylim([min_acpr(r) max_acpr(r)]*1.1/nd_a);
        grid on;
        set(gca,'GridColor',[0 0 0]);

        % Slip and slip rate vs. time
        f200_ax(5) = subplot(2,spn,spn+2);
        yyaxis right;
        plot(time_pr/nd_t, disp_pr/nd_u, 'LineWidth', lw);
        xlabel('$t_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
        ylim([0 max_upr(r)*1.03]/nd_u);
        set(gca,'FontSize',fs-2);
        if idx_pc == idx_cei
            ylabel('$D_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
        end 
        yyaxis left;
        plot(time_pr/nd_t, vel_pr/nd_v, 'LineWidth', lw);
        axis tight;
        grid on;
        set(gca,'GridColor',[0 0 0]);
        ylim([0 v_pc*1.03]/nd_v);
        legend([p(2) p(3)],'Friction','Loading', 'Interpreter', 'latex', 'FontSize', fs-1, 'location', 'southwest');

        % Coseismic
        if idx_pc ~= idx_cei

            % Stress and acceleration vs. time
            f200_ax(3) = subplot(2,spn,3);
            yyaxis right;
            p(2) = plot(time_co/nd_t, (tau_co/Sn+mu_i-mu_d)/(mu_s-mu_d), 'LineWidth', lw);
            hold on;
            p(3) = plot(time_co/nd_t, (ku_co/Sn+mu_i-mu_d)/(mu_s-mu_d), 'LineWidth', lw);
            xline((time_co(idx_fl))/nd_t, 'k-.', 'LineWidth', lw-0.5);
            set(gca,'FontSize',fs-2);
            ylabel('$\tau_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
            ylim([min_Sco(r) max_Sco(r)]*1.02);
            axis tight;
            yyaxis left;
            p(1) = plot(time_co/nd_t, ac_co/nd_a, 'LineWidth', lw);
            ylim([min_acco(r) max_acco(r)]*1.1/nd_a);
            set(gca,'xticklabel',[]);
            hold on;
            title('Coseismic', 'Interpreter', 'latex', 'FontSize', fs+1, 'FontWeight', 'bold');
            axis tight;
            grid on;
            set(gca,'GridColor',[0 0 0]);

            % Slip and slip rate vs. time
            f200_ax(6) = subplot(2,spn,2*spn);
            yyaxis right;
            plot(time_co/nd_t, disp_co/nd_u, 'LineWidth', lw);
            hold on;
            plot((t_sd-max(t_sd)/2+time_co(I_m))/nd_t, u_sd_1/nd_u, ':', 'LineWidth', lw);
            ylim([0 max(u_sd)*1.03]/nd_u);
            ylabel('$D_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
            xlabel('$t_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
            yyaxis left;
            plot(time_co/nd_t, vel_co/nd_v, 'LineWidth', lw);
            hold on;
            plot((t_sd-max(t_sd)/2+time_co(I_m))/nd_t, v_sd_1/nd_v, ':', 'LineWidth', lw);
            xline((time_co(idx_fl))/nd_t, 'k-.', 'LineWidth', lw-0.5);
            set(gca,'FontSize',fs-2);
            ylim([0 max(v_sd)*1.03]/nd_v);
            axis tight;
            grid on;
            set(gca,'GridColor',[0 0 0]);
            legend([p(2) p(3)],'Friction','Loading', 'Interpreter', 'latex', 'FontSize', fs-1, 'location', 'southwest');
        end

    end

end

%% Spectra
% Add harmonic solution to arrays for spectra
v_sp_a(1:length(v_sd_1), np+1) = v_sd_1;
t_sp_a(1:length(v_sd_1), np+1) = t_sd;

% Calculate frequency and padded time
t_sp_max    = 50*max(max(t_sp_a));
dta         = diff(t_sp_a,1);
dta(dta<=0) = NaN;
Fs          = 1/min(min(dta));
time_sp     = 0:min(min(dta)):t_sp_max;
L           = length(time_sp);
if mod(L,2) == 1
    time_sp = 0:min(dta):t_sp_max + min(dta);
    L       = L+1;
end

v_sp_pad    = zeros(length(time_sp),np+1);
spec        = zeros(L/2,np+1);

% Compute FFT
for i = 1:np+1
    v_sp_pad(:,i) = interp1(t_sp_a(~isnan(t_sp_a(:,i)),i)+t_sp_max/2, v_sp_a(~isnan(v_sp_a(:,i)),i), time_sp, 'linear',0);
    Y             = fft(v_sp_pad(:,i),L);
    P2            = abs(Y/L);
    P1            = P2(1:L/2+1);
    P1(2:end-1)   = 2*P1(2:end-1);
    spec(:,i)     = P1(1:L/2);
end

% Plot spectra
figure(115);
pl(np+1) = loglog((0:(Fs/L):(Fs/2-Fs/L))*2*pi/sqrt(K_s/M),spec(:,np+1)/max(spec(:,np+1)), 'Color', [0.75, 0.75, 0], 'LineWidth', lw, 'DisplayName', '"Reference" $D_{c}=0$');
hold on;
for fl = 1:np
    pl(fl) = loglog((0:(Fs/L):(Fs/2-Fs/L))*2*pi/sqrt(K_s/M),spec(:,fl)/max(spec(:,np+1)), 'Color', cfl(fl,:), 'LineWidth', lw, 'DisplayName',lower(char(fl+64)));
end
loglog(logspace(-1,5,100),1./(logspace(-2,4,100)), 'k:', 'LineWidth', lw/2, 'DisplayName', '$f_{nd}^{\ \ -1}$');
loglog(logspace(-1,5,100),1./(logspace(-2,4,100)).^2, 'k--', 'LineWidth', lw/2, 'DisplayName', '$f_{nd}^{\ \ -2}$');
set(gca,'FontSize',fs-2);
ylabel('Amplitude', 'Interpreter', 'latex', 'FontSize', fs);
xlabel('$f_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
legend('show', 'Interpreter', 'latex', 'FontSize', fs-1, 'location', 'southwest');
xlim([1e-2 1e3]);
ylim([1e-8 2]);
grid on;


%% Plot energy curves
figure(112);
cfl      = [b_col; r_col; g_col; k_col; lb_col; o_col; y_col; p_col; dr_col];                % colors
pl(np+1)  = plot([D_s_arr(2),D_s_arr(2),1.3*nd_u]/nd_u, [1 0 0], 'Color', [0.75, 0.75, 0], 'LineWidth', lw, 'DisplayName', '"Reference" $D_{c}=0$');
hold on;
for fl = 1:np
    sl_v  = cumsum([0 D_s_arr(fl) D_t_arr(fl) D_w_arr(fl) 3*D_w_arr(fl)]);
    mu_v  = [mu_i_arr(fl) mu_s_arr(fl) mu_t_arr(fl) mu_d_arr(fl) mu_d_arr(fl)];
    pl(fl) = plot(sl_v/nd_u, (mu_v-mu_d_arr(1))/(mu_s_arr(1)-mu_d_arr(1)), 'Color', cfl(fl,:), 'LineWidth', lw, 'DisplayName',lower(char(fl+64)));
    
end
pl(np+2) = plot(linspace(D_s_arr(2),1.3*nd_u,10)/nd_u, -K_s*linspace(0,1.3*nd_u-D_s_arr(2),10)/Sn/(mu_s_arr(1)-mu_d_arr(1))+1, 'k--', 'LineWidth', lw, 'DisplayName', 'Loading');
set(gca,'FontSize',fs-2);
ylabel('$\tau_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
xlabel('$D_{nd}$', 'Interpreter', 'latex', 'FontSize', fs);
hl = legend('show', 'Interpreter', 'latex', 'FontSize', fs-1);
grid on;
xlim([0 1.3]);
ylim([-0.02 1.02]);
