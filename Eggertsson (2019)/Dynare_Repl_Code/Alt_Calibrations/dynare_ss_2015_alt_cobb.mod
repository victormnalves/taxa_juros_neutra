%--------------------------------------------------------------------------
% Authors: Alex Crescentini & Federico Giri.
% Università Politecnica delle Marche, Ancona, Italy.
% January 2023.
%--------------------------------------------------------------------------
% Dynare Replication Code for:
% "A Model of Secular Stagnation: Theory and Quantitative Evaluation"
% Eggertsson, Mehrotra and Robbins (2019), American Economic Journal.
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
//DESCRIPTION
%--------------------------------------------------------------------------

%Code for Steady State at 2015

global run_alt_number

%--------------------------------------------------------------------------
//LOAD DATA TO BE USED IN THE MODEL
%--------------------------------------------------------------------------

datasurvival1_1970 = readmatrix('data/data_survival_1_1970.xlsx');
datasurvival2_1970 = readmatrix('data/data_survival_2_1970.xlsx');
data_uncond_survival_1970 = readmatrix('data/data_uncondit_survival_1970.xlsx');
datasurvival1_2015 = readmatrix('data/data_survival_1_2015.xlsx');
datasurvival2_2015 = readmatrix('data/data_survival_2_2015.xlsx');
data_uncond_survival_2015 = readmatrix('data/data_uncondit_survival_2015.xlsx');
datahc = readmatrix('data/data_hc.xlsx');
data_debt_1970 = readmatrix('data/data_debt_1970.xlsx');
data_debt_2015 = readmatrix('data/data_debt_2015.xlsx');
ss_1970 = xlsread('data/ss_1970.xlsx');
ss_2015 = xlsread('data/ss_2015.xlsx');
popgr_1970 = xlsread('data/popgr_1970.xlsx');
%popgr_2015 = xlsread('data/popgr_2015.xlsx');

%--------------------------------------------------------------------------
//GLOBAL VARIABLES
%--------------------------------------------------------------------------

@#define J = 56
@#define T = 152

%--------------------------------------------------------------------------
//ENDOGENOUS VARIABLES
%--------------------------------------------------------------------------

var 

@#for j in 1:J
	n@{j}
@#endfor

@#for j in 1:J
	c@{j}
@#endfor

@#for j in 1:J+1
	a@{j}
@#endfor

@#for j in 1:J
	lambda@{j}
@#endfor

@#for j in 1:J-16
	pi@{j}
@#endfor

x56, q32, w, rk, r, price, tau, PI, Y, N, L, C, K, gov_rev, gov_deficit
gov_debt, I, IY, LS

;

%--------------------------------------------------------------------------
//PARAMETERS
%--------------------------------------------------------------------------

parameters

alpha, beta, delta, gamma, sigma, mu, g,

@#for j in 1:J-16
	hc@{j}
@#endfor

@#for j in 1:J
    s@{j}
@#endfor

@#for j in 1:J
    sv@{j}
@#endfor

@#for j in 1:J
    su@{j}
@#endfor

n, fert_26, e, b, D, theta, AL_growth, AL, AK, A_adj, popgr1970
    
;

%--------------------------------------------------------------------------
//PARAMETERS' CALIBRATION
%-------------------------------------------------------------------------- 

alpha=param_nc(1);
beta=param_nc(2);
delta=param_nc(3);
gamma=param_nc(4);
sigma=param_nc(5); %IF WE SET IT = 1, IT BECOMES A COBB-DOUGLAS 
mu=param_nc(6);
g=param_nc(7);     %government spend as % of gdp

@#for j in 1:40
	hc@{j}=datahc(1,@{j});
@#endfor

@#for j in 1:J
    s@{j}=datasurvival1_2015(1,@{j});
@#endfor
    
@#for j in 1:J
    sv@{j}=datasurvival2_2015(1,@{j});
@#endfor
        
@#for j in 1:J
    su@{j}=data_uncond_survival_2015(1,@{j});
@#endfor
            
n=param_2015(1);                %2015
fert_26=param_2015(2);          %2015 ENDOG.FERTILITY.m
e=param_2015(3);                %2015 e=1/z with z=1
b=param_2015(4);                %2015 government debt as % of gdp
D=param_2015(5);                %2015
theta=param_2015(6);            %2015
AL_growth= param_2015(7);       %2015
AL=param_2015(8);               %2015
AK=param_2015(9);               %2015
A_adj=ss_1970(283);             %1970
popgr1970=popgr_1970;

%--------------------------------------------------------------------------
//MODEL BLOCK: EQUATIONS
%--------------------------------------------------------------------------

model;

//Demographics
%--------------------------------------------------------------------------

[name = 'Demographics Equation 1']
%n1=(1+n)^151;
n1=0.965569227819326; %directly from eggertsson for simplicity

@#for j in 1:J-1
    [name = 'Demographics Equation @{j+1}']
    n@{j+1} = (s@{j}*n@{j})/(1+n);
@#endfor

//Households: First-Order conditions
%--------------------------------------------------------------------------

@#for j in 1:J-1
    [name = 'Households First-order Condition @{j}']
    (1/beta) = ((c@{j+1}/c@{j})^(-1/gamma))*(1+r) + (lambda@{j+1})*((c@{j})^(1/gamma))*sv@{j+1}/(su@{j}*beta^@{j});
@#endfor

    [name = 'Households First-order Condition 56']
    x56 = (fert_26/mu)^(-gamma)*c56;
    
//Households: Budget Constraints
%--------------------------------------------------------------------------

[name = 'Households Budget Constraint: initial condition']
a1=0;

@#for j in 1:J-26
    [name = 'Households Budget Constraint @{j}']
    a@{j+1} = (1+r)/(sv@{j})*a@{j} + (pi@{j}+(1-tau)*w*hc@{j})*(1+AL_growth)^@{j} - c@{j};
@#endfor

[name = 'Households Budget Constraint 31']
a32 = (1+r)/sv31*a31 + (pi31+(1-tau)*w*hc31)*(1+AL_growth)^@{31} + q32*(1+AL_growth)^@{32} - c31;

[name = 'Bequest equation']
q32 = ((x56*fert_26*n56)/(n32)*(1/(1+AL_growth)^57))*1/(1+n);
%q32 = ((x56*fert_26*n56)/(n32)*(1/(1+AL_growth)^57))*1/(1+popgr1970);

@#for j in 32:J-16
    [name = 'Households Budget Constraint @{j}']
    a@{j+1} = (1+r)/sv@{j}*a@{j} + (pi@{j}+(1-tau)*w*hc@{j})*(1+AL_growth)^@{j} - c@{j};
@#endfor

@#for j in 41:J-1
    [name = 'Households Budget Constraint @{j}']
    a@{j+1} = (1+r)/sv@{j}*a@{j} - c@{j};
@#endfor
    
    [name = 'Households Budget Constraint 56']
    a57 = (1+r)/sv56*a56 - (fert_26*x56) - c56;

[name = 'Households Budget Constraint: terminal condition']
a57 = 0;

//Households: Financial Constraints (OBCs)
%--------------------------------------------------------------------------

@#for j in 1:J-16
    [name = 'Households Financial Constraint @{j}']
    min(lambda@{j},(a@{j}+((D*w*hc@{j})*(1+AL_growth)^@{j}))) = 0;  
@#endfor
@#for j in 41:J
    [name = 'Households Financial Constraint @{j}']
    min(lambda@{j},(a@{j})) = 0;  
@#endfor

//Households: Profits from Firms
%--------------------------------------------------------------------------

@#for j in 1:J-16
    [name = 'Households Profits from Firms @{j}']
    pi@{j}=hc@{j}*PI/L;
@#endfor

//Firms
%--------------------------------------------------------------------------

price = (theta-1)/theta;
PI = Y/theta;

%%COBB-DOUGLAS    
w = ((price)*(1-alpha)*AL*((AK*K)/(AL*L))^(alpha))/A_adj;
rk = ((price)*alpha*AK*((AK*K)/(AL*L))^(alpha-1))/A_adj;
Y = ((AK*K)^(alpha)*(AL*L)^(1-alpha))/A_adj;
r = rk/e + (1-delta)*e/e - 1;

 //Government
%--------------------------------------------------------------------------

[name = 'Government']
gov_deficit*gov_rev = ((1+AL_growth)*(1+n)-1)*(gov_debt*K);
gov_debt=b*Y/K;
gov_rev = (g*Y+r*gov_debt*K);
tau*w = ((gov_rev)*(1-gov_deficit))/(L);

//Aggregates
%--------------------------------------------------------------------------

[name = '(N) Population']
N = (
@#for j in 1:J
    + n@{j}
@#endfor
);

[name = '(L) Labor']
L = (
@#for j in 1:J-16
    + n@{j}*hc@{j}
@#endfor
);

[name = '(C) Consumption']
C = (
@#for j in 1:J
    + n@{j}*c@{j}/(1+AL_growth)^(@{j})
@#endfor
);

[name = '(K) Kapital']
K = ((
@#for j in 1:J
    + n@{j}*a@{j}/sv@{j}*1/(1+AL_growth)^(@{j})
@#endfor
)/(e+gov_debt));

[name = '(I) Investment']
I = (1+AL_growth)*(1+n)*e*K - (1-delta)*e*K;

[name = 'Investment Output Ratio (I/Y)']
IY = I/Y;

[name = 'Labor Share']
LS = L/(L+r*K+PI+delta*K);

end;

%--------------------------------------------------------------------------
//STEADY-STATE BLOCK: Initial Values
%--------------------------------------------------------------------------

initval;

@#for j in 1:J
    n@{j}=ss_2015(@{j});
@#endfor
    
@#for j in 1:J
	c@{j}=ss_2015(@{j+56});
@#endfor
        
@#for j in 1:J+1
	a@{j}=ss_2015(@{j+112});
@#endfor
            
@#for j in 1:J
	lambda@{j}=ss_2015(@{j+169});
@#endfor
                
@#for j in 1:J-16
	pi@{j}=ss_2015(@{j+225});
@#endfor

x56=ss_2015(266);
q32=ss_2015(267);
w=ss_2015(268);
rk=ss_2015(269);
r=ss_2015(270);
price=ss_2015(271);
tau=ss_2015(272);
PI=ss_2015(273);
Y=ss_2015(274);
N=ss_2015(275);
L=ss_2015(276);
C=ss_2015(277);
K=ss_2015(278);
gov_rev=ss_2015(279);
gov_deficit=ss_2015(280);
gov_debt=ss_2015(281);
I=ss_2015(282);

end;

%--------------------------------------------------------------------------
//STEADY-STATE BLOCK: Computation
%--------------------------------------------------------------------------

%options_.debug=1;
steady(maxit=1000, solve_algo=0);
%resid;
%check;

%--------------------------------------------------------------------------
//FILES PREPERATION FOR TRANSITION DYNAMICS: "dynare_transition.mod"
%--------------------------------------------------------------------------

%%2015: SCALING FOR PRODUCTIVITY GROWTH%%

J=56;
AL_growth=0.006460448230438;
%AL_growth=0; 

cons=oo_.steady_state(57:112,1);
capi=oo_.steady_state(113:169,1);
for j = 1:J
	c_scaled_2015(j,1)=cons(j,1)/(1+AL_growth)^(j);
end
filename='data/c_scaled_2015.xlsx';
writematrix(c_scaled_2015,filename,'WriteMode','replacefile');
J=57;
for j = 1:J
	a_scaled_2015(j,1)=capi(j,1)/(1+AL_growth)^(j-1);
end
filename='data/a_scaled_2015.xlsx';
writematrix(a_scaled_2015,filename,'WriteMode','replacefile');
x56_scaled_2015=oo_.steady_state(266,1)/(1+AL_growth)^56;
filename='data/x56_scaled_2015.xlsx';
writematrix(x56_scaled_2015,filename,'WriteMode','replacefile');
q32_scaled_2015=oo_.steady_state(267,1);
%q32_scaled_2015=oo_.steady_state(267,1)/(1+AL_growth);
filename='data/q32_scaled_2015.xlsx';
writematrix(q32_scaled_2015,filename,'WriteMode','replacefile');
filename='data/ss_2015.xlsx';
writematrix(oo_.steady_state,filename,'WriteMode','replacefile');
%filename='data/popgr_2015.xlsx';
%writematrix(M_.params(216),filename,'WriteMode','replacefile');

%Consumer Debt to GDP
datasurvival2_2015 = readmatrix('data/data_survival_2_2015.xlsx');
popgen=oo_.steady_state(1:56,1);
output=oo_.steady_state(274,1);
personal_debt=0;
J=56;
for j = 1:J
	a_scaled2_2015(j,1)=capi(j,1)/(1+AL_growth)^(j);
end
for j = 1:J
	personal_debt=personal_debt+(a_scaled2_2015(j,1)<0)*a_scaled2_2015(j,1)*popgen(j,1)/datasurvival2_2015(1,j);
end
CDY=-personal_debt/output;
filename='data/CDY_2015.xlsx';
writematrix(CDY,filename,'WriteMode','replacefile');

%--------------------------------------------------------------------------
//LATEX FILES
%--------------------------------------------------------------------------

%write_latex_dynamic_model;
%write_latex_static_model;
%write_latex_definitions;
%write_latex_parameter_table;
%collect_latex_files;