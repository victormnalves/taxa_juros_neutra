% Secular stagnation model with productivity growth (Eurozone calibration)
%
% Authors: Gauti Eggertsson, Neil Mehrotra, and Jacob Robbins
% Date: June 6, 2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
clc;

% Aggregate demand at ZLB
syms beta popgrw 
syms d df dl dlf mu muf out outf v vf

f1 = (1/vf) - (((1+beta)/beta)*popgrw*d)/(outf-(dl/mu));

% Hysteresis equation
syms kappa y_fe mu_0;

f2 = mu - mu_0*(outf/y_fe)^kappa;

% Aggregate supply with downward wage rigidity
syms phi gamma pi_tar

f3 = (outf^phi) - ((gamma*pi_tar/mu)*(out^phi)/v) - ((1-gamma)*(y_fe^phi));

% Law of motion for debt constraint
syms d_ss rho_d

f4 = log(df/d_ss) - rho_d*log(d/d_ss);
f5 = dlf - d;

f = [f1;f2;f3;f4;f5];

% Define state and control variable
y = [v,mu];
yf = [vf,muf];

x = [out,d,dl];
xf = [outf,df,dlf];

f = subs(f,[x,y,xf,yf],exp([x,y,xf,yf]));

fx  = jacobian(f,x);
fxf = jacobian(f,xf);
fy  = jacobian(f,y);
fyf = jacobian(f,yf);

% Calibration of initial steady state

beta  = 0.96^20;     % rate of time preference
mu_0  = 1.018^20;    % gross productivity growth (pre-shock)
alpha = 0.7;         % labor share
popgrw = 1;          % gross population growth
y_fe   = 1;          % full-employment output
rho_d  = 1;          % persistence of debt constraint shock

phi = (alpha-1)/alpha;

% Calibration targets
y_ss = 0.9;          % output gap in secular stagnation 
pi_ss = 1^20;        % gross inflation rate in secular stagnation
pi_tar = 1.02^20;    % gross inflation target

mu_ss = 1.0022^20;   % targeted post-hysteresis productivity growth rate
kappa = log(mu_ss/mu_0)/log(y_ss/y_fe);
con1 = ((1+beta)/beta)*popgrw + (1/(pi_ss*mu_ss));

d_ss = (y_ss/pi_ss)/con1;
gamma = (1-(y_ss^phi))/(1-((y_ss^phi)*pi_tar/(pi_ss*mu_ss)));

% Steady state variable
out  = log(y_ss);
v    = log(pi_ss);
d    = log(d_ss);
dl   = log(d_ss);
mu   = log(mu_ss);

outf = out;
vf   = v;
df   = d;
dlf  = dl;
muf  = mu;

fx = eval(fx);
fy = eval(fy);
fxf = eval(fxf);
fyf = eval(fyf);

f = eval(f);

[gx,hx,exitflag] = gx_hx(fy,fx,fyf,fxf,1.01) % solution for linearized system

% Initial state
out_ini = log(1/y_ss);

nat_rate = 1.005^20;

con2 = ((1+beta)/beta)*popgrw + (nat_rate/mu_0);
d_ini = nat_rate/con2;

dl_ini = log(d_ini/d_ss);

x0 = zeros(3,1);
x0(1) = out_ini;
x0(3) = dl_ini;

IR = ir(gx,hx,x0,21);               % impulse response

pi_path = pi_ss*exp(IR(:,1));       % inflation transition path (EU spreadsheet, column N, Figure 5.xlsx)
prod_path = mu_ss*exp(IR(:,2));     % productivity transition path (EU spreadsheet, column K, Figure 5.xlsx)
out_path = y_ss*exp(IR(:,3));       % output transition path (EU spreadsheet, column N, Figure 5.xlsx)

IR = 100*IR;