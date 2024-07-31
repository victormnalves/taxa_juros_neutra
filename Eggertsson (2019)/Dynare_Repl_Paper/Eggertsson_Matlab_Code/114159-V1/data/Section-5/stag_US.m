% Secular stagnation model with productivity growth (US calibration)
%
% Authors: Gauti Eggertsson, Neil Mehrotra, and Jacob Robbins
% Date: June 6, 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
clc;

% Aggregate demand at ZLB
syms beta popgrw y_fe mu_0 tau_ss bg_ss sy sm
syms d df dl dlf out outf v vf

f1 = (sy*d*vf)+bg_ss - (beta/(1+beta))*sm*(outf-tau_ss-(dl/mu_0));

% Aggregate supply with downward wage rigidity
syms phi gamma pi_tar

f2 = (outf^phi) - ((gamma*(pi_tar/mu_0))*(out^phi)/v) - ((1-gamma)*(y_fe^phi));

% Law of motion for debt constraint
syms d_ss rho_d

f3 = log(df/d_ss) - rho_d*log(d/d_ss);
f4 = dlf - d;

f = [f1;f2;f3;f4];

% Define state and control variable
y = v;
yf = vf;

x = [out,d,dl];
xf = [outf,df,dlf];

f = subs(f,[x,y,xf,yf],exp([x,y,xf,yf]));

fx  = jacobian(f,x);
fxf = jacobian(f,xf);
fy  = jacobian(f,y);
fyf = jacobian(f,yf);

% Calibration of initial steady state

beta   = 0.96^20;
mu_0   = 1.021^20;
alpha  = 0.7;
popgrw = 1.007^20;
rho_d  = 1;

phi    = (alpha-1)/alpha;

% Targets for calibration
y_ss   = 0.87;              % US output gap
pi_ss  = 1.014^20;          % inflation rate (core PCE growth) in 2014
pi_tar = (1.02*1.021)^20;   % wages indexed by inflation target and productivity

mu_ss = 1.021^20;
kappa = log(mu_ss/mu_0)/log(y_ss/y_fe);

% Calibration of income distribution
sy = (popgrw^2)/(1+popgrw+popgrw^2);
sm = popgrw/(1+popgrw+popgrw^2);

g_ss  = 0.2;                % government spending (% of GDP)
bg_ss = 1/20;               % public debt (1x annual GDP)
bp_ss = 0.05;               % private debt (1x annual GDP)

tau_ss  = (g_ss + bg_ss*((1/(pi_ss*popgrw*mu_0))-1))/sm;
d_ss    = bp_ss/(pi_ss*sy);
ym_ss   = (bp_ss + bg_ss + sm*(beta/(1+beta))*(tau_ss+(d_ss/mu_0)))/(sm*beta/(1+beta));
y_fe    = ym_ss/y_ss;

gamma = (1-(y_ss^phi))/(1-((y_ss^phi)*pi_tar/(pi_ss*mu_ss)));

% Steady state variable
out  = log(ym_ss);
v    = log(pi_ss);
d    = log(d_ss);
dl   = log(d_ss);
%mu   = log(mu_ss);

outf = out;
vf   = v;
df   = d;
dlf  = dl;
%muf  = mu;

fx = eval(fx);
fy = eval(fy);
fxf = eval(fxf);
fyf = eval(fyf);

f = eval(f);

[gx,hx,exitflag] = gx_hx(fy,fx,fyf,fxf,1.01)    % solution for linearized system

% Initial state

out_ini  = log(y_fe/ym_ss);
nat_rate = 1.009^20;

tau_ini  = (g_ss + bg_ss*((nat_rate/(popgrw*mu_0))-1))/sm;
con2     = (sy/nat_rate)+(sm*beta/((1+beta)*mu_0));
d_ini    = (((beta/(1+beta))*sm*(y_fe-tau_ini))-bg_ss)/con2;

dl_ini = log(d_ini/d_ss);

x0 = zeros(3,1);
x0(1) = out_ini;
x0(3) = dl_ini;

IR = ir(gx,hx,x0,21);               % impulse response

pi_path = pi_ss*exp(IR(:,1));       % inflation transition path (US spreadsheet, column N, Figure 5.xlsx)
out_path = y_ss*exp(IR(:,2));       % output transition path (US spreadsheet, column N, Figure 5.xlsx

IR = 100*IR;

% Income distribution and consumption

yy_ss = (y_ss-(sm*ym_ss))/sy;
cy_ss = yy_ss + (d_ss*pi_ss);
cm_ss = (1/(1+beta))*(ym_ss-tau_ss-(d_ss/mu_0));
bm_ss = beta*cm_ss;
co_ss = bm_ss/(pi_ss*mu_0);