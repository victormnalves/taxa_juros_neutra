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

%Code for Steady State at 1970

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

x56, q32, w, rk, r, price, tau, PI, Y, N, L, C, K, gov_rev, gov_deficit,
gov_debt, I, A_adj

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

n, fert_26, e, b, D, theta, AL_growth, AL, AK

;

%--------------------------------------------------------------------------
//PARAMETERS' CALIBRATION
%-------------------------------------------------------------------------- 

alpha=0.238808876581195;
beta=0.979996368002075;
delta=0.1244;
gamma=0.75;
sigma=0.6; %IF WE SET IT = 1, IT BECOMES A COBB-DOUGLAS AND WE NEED TO CHANGE SOME EQUATIONS
mu=21.629099653217530;
g=0.2128;     %government spend as % of gdp

@#for j in 1:40
	hc@{j}=datahc(1,@{j});
@#endfor

@#for j in 1:J
    s@{j}=datasurvival1_1970(1,@{j});
@#endfor
    
@#for j in 1:J
    sv@{j}=datasurvival2_1970(1,@{j});
@#endfor
        
@#for j in 1:J
    su@{j}=data_uncond_survival_1970(1,@{j});
@#endfor
            
n=0.013549868016229;               %1970
fert_26=1.4;                       %1970 ENDOG.FERTILITY.m
e=1.3;                             %1970 e=1/z with z=0.7692
b=0.422844989041041;               %1970 government debt as % of gdp
D=0.144020808265552;               %1970
theta=7.928343604268710;           %1970
AL_growth = 0.020209202089143; %1970 
AL=1;                              %1970
AK=1;                              %1970
    
%--------------------------------------------------------------------------
//MODEL BLOCK: EQUATIONS
%--------------------------------------------------------------------------

model;

//Demographics
%--------------------------------------------------------------------------

[name = 'Demographics Equation 1']
n1=1;

@#for j in 1:J-1
    [name = 'Demographics Equation @{j+1}']
    n@{j+1} = (s@{j}*n@{j})/(1+n);
@#endfor

//Households: First-Order conditions
%--------------------------------------------------------------------------

@#for j in 1:J-1
    [name = 'Households First-order Condition @{j}']
    (1/beta) = ((c@{j+1}/c@{j})^(-1/gamma))*(1+r) + (lambda@{j+1})*((c@{j})^(1/gamma))*sv@{j+1}/(su@{j+1}*beta^@{j}*e);
@#endfor

    [name = 'Households First-order Condition 56']
    x56 = (fert_26/mu)^(-gamma)*c56;

//Households: Budget Constraints
%--------------------------------------------------------------------------

[name = 'Households Budget Constraint: Initial Condition']
a1=0;

@#for j in 1:J-26
    [name = 'Households Budget Constraint @{j}']
    e*a@{j+1} = (rk+e*(1-delta))*a@{j}/sv@{j} + (pi@{j}+(1-tau)*w*hc@{j})*(1+AL_growth)^@{j} - c@{j};
@#endfor

[name = 'Households Budget Constraint 31']
e*a32 = (rk+e*(1-delta))*(a31/sv31+q32*(1+AL_growth)^@{32}) + (pi31+(1-tau)*w*hc31)*(1+AL_growth)^@{31} - c31;

[name = 'Bequest equation']
q32 = ((x56*fert_26*n56)/(n32)*(1/(1+AL_growth)^57))*1/(1+n);

@#for j in 32:J-16
    [name = 'Households Budget Constraint @{j}']
    e*a@{j+1} = (rk+e*(1-delta))*a@{j}/sv@{j} + (pi@{j}+(1-tau)*w*hc@{j})*(1+AL_growth)^@{j} - c@{j};
@#endfor

@#for j in 41:J-1
    [name = 'Households Budget Constraint @{j}']
    e*a@{j+1} = (rk+e*(1-delta))*a@{j}/sv@{j} - c@{j};
@#endfor
    
[name = 'Households Budget Constraint 56']
e*a57 = (rk+e*(1-delta))*a56/sv56 - (fert_26*x56) - c56;

[name = 'Households Budget Constraint: Terminal Condition']
a57 = 0;

//Households: Financial Constraints (OBCs)
%--------------------------------------------------------------------------

@#for j in 1:J-16
    [name = 'Households Financial Constraint @{j}']
    min(lambda@{j},a@{j}+(D/(1+r))*(1+AL_growth)^@{j}) = 0; 
@#endfor

@#for j in 41:J
    [name = 'Households Financial Constraint @{j}']
    min(lambda@{j},a@{j}) = 0;
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
A_adj = (price*(((alpha)*(AK*K)^((sigma-1)/sigma)+(1-alpha)*(AL*L)^((sigma-1)/sigma))^(1/(sigma-1)))*(1-alpha)*AL^((sigma-1)/sigma)*L^(-1/sigma));
w  = 1;
rk = (price*(((alpha)*(AK*K)^((sigma-1)/sigma)+(1-alpha)*(AL*L)^((sigma-1)/sigma))^(1/(sigma-1)))*(alpha)*AK^((sigma-1)/sigma)*K^(-1/sigma))/A_adj;
Y  = (((alpha)*(AK*K)^((sigma-1)/sigma)+(1-alpha)*(AL*L)^((sigma-1)/sigma))^(sigma/(sigma-1)))/A_adj;
r = rk/e + (1-delta)*e/e - 1;
PI = Y/theta;

 //Government
%--------------------------------------------------------------------------

[name = 'Government']
b*Y*((1+AL_growth)*(1+n)) = g*Y + (1+r)*b*Y - tau*w*L;
gov_rev=g*Y+r*b*Y;
gov_deficit=((b*Y)*((1+AL_growth)*(1+n)-1))/(gov_rev);
gov_debt=b*Y/K;

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
e*K = (
@#for j in 1:J
    + n@{j}*e*a@{j}*1/(1+AL_growth)^(@{j})
@#endfor
)-(b*Y);

[name = '(I) Investment']
I = (1+AL_growth)*(1+n)*e*K - (1-delta)*e*K;

end;

%--------------------------------------------------------------------------
//STEADY-STATE BLOCK: Initial Values
%--------------------------------------------------------------------------

initval;

@#for j in 1:J
    n@{j}=ss_1970(@{j});
@#endfor

@#for j in 1:J
	c@{j}=ss_1970(@{j+56});
@#endfor

@#for j in 1:J+1
	a@{j}=ss_1970(@{j+112});
@#endfor

@#for j in 1:J
	lambda@{j}=ss_1970(@{j+169});
@#endfor

@#for j in 1:J-16
	pi@{j}=ss_1970(@{j+225});
@#endfor

x56=ss_1970(266);
q32=ss_1970(267);
w=ss_1970(268);
rk=ss_1970(269);
r=ss_1970(270);
price=ss_1970(271);
tau=ss_1970(272);
PI=ss_1970(273);
Y=ss_1970(274);
N=ss_1970(275);
L=ss_1970(276);
C=ss_1970(277);
K=ss_1970(278);
gov_rev=ss_1970(279);
gov_deficit=ss_1970(280);
gov_debt=ss_1970(281);
I=ss_1970(282);
A_adj=ss_1970(283);

end;

%--------------------------------------------------------------------------
//STEADY-STATE BLOCK: Computation
%--------------------------------------------------------------------------

%options_.debug=1;
steady(maxit=1000, solve_algo=0);
%resid;
%check;

%--------------------------------------------------------------------------
//OUTPUT AS IN AUERBACH & KOTLIKOFF (1987) AND EGGERTSSON ET AL. (2019)
%--------------------------------------------------------------------------

disp(oo_.steady_state(266:278,1));

disp(['   '])
disp(['DYNARE-Statistics'])
disp(['   '])
disp(['Year  ' 'Capital (K)  ' 'Labor (L)  ' 'Population (N)'])
disp([num2str(0) '     ' num2str(oo_.steady_state(278,1)) '       ' num2str(oo_.steady_state(276,1)) '     ' num2str(oo_.steady_state(275,1))])

disp(['   '])
disp (['Income (Y)  ' 'Consumption (C)  ' 'Agg. Profit (PI)' 'Investment (I)'])
disp([num2str(oo_.steady_state(274,1)) '     ' num2str(oo_.steady_state(277,1)) '          ' num2str(oo_.steady_state(273,1)) '    ' num2str(oo_.steady_state(282,1))])

disp(['   '])
disp(['Wage (w)  ' 'Rental K (rk)  ' 'Interest (r)  '   'Wage Tax (tau)'])
disp([num2str(oo_.steady_state(268,1)) '         ' num2str(oo_.steady_state(269,1)) '        ' num2str(oo_.steady_state(270,1)) '      ' num2str(oo_.steady_state(272,1))])

disp(['   '])
disp(['Bequest (q32*fert_26)  ' 'q32     ' 'x56/(1+AL_growth)^56'])
disp([num2str(oo_.steady_state(267,1)*fert_26) '                ' num2str(oo_.steady_state(267,1)) '  ' num2str(oo_.steady_state(266,1)/(1+AL_growth)^56)])

disp(['   '])
disp(['Pop Growth (n)  '])
disp([num2str(n)])

disp(['   '])
disp(['Debt (b*Y/K)  '])
disp([num2str(b*(oo_.steady_state(274,1)/oo_.steady_state(278,1)))])

disp(['   '])
disp(['G (g*Y)  '])
disp([num2str(g*oo_.steady_state(274,1))])

disp(['   '])
disp(['Y = C + I + g*Y '])
disp([num2str(oo_.steady_state(277,1) + oo_.steady_state(282,1) + g*oo_.steady_state(274,1))])

%--------------------------------------------------------------------------
//FILES PREPERATION FOR TRANSITION DYNAMICS: "dynare_transition.mod"
%--------------------------------------------------------------------------

%%1970: SCALING FOR PRODUCTIVITY GROWTH%%

J=56;
AL_growth=0.020209202089143; 

cons=oo_.steady_state(57:112,1);
capi=oo_.steady_state(113:169,1);
for j = 1:J
	c_scaled_1970(j,1)=cons(j,1)/(1+AL_growth)^(j);
end
filename='data/c_scaled_1970.xlsx';
writematrix(c_scaled_1970,filename,'WriteMode','replacefile');
J=57;
for j = 1:J
	a_scaled_1970(j,1)=capi(j,1)/(1+AL_growth)^(j-1);
end
filename='data/a_scaled_1970.xlsx';
writematrix(a_scaled_1970,filename,'WriteMode','replacefile');
x56_scaled_1970=oo_.steady_state(266,1)/(1+AL_growth)^56;
filename='data/x56_scaled_1970.xlsx';
writematrix(x56_scaled_1970,filename,'WriteMode','replacefile');
q32_scaled_1970=oo_.steady_state(267,1);
%q32_scaled_1970=oo_.steady_state(267,1)/(1+AL_growth);
filename='data/q32_scaled_1970.xlsx';
writematrix(q32_scaled_1970,filename,'WriteMode','replacefile');
filename='data/ss_1970.xlsx';
writematrix(oo_.steady_state,filename,'WriteMode','replacefile');
filename='data/popgr_1970.xlsx';
writematrix(M_.params(216),filename,'WriteMode','replacefile');

%--------------------------------------------------------------------------
//LATEX FILES
%--------------------------------------------------------------------------

%write_latex_dynamic_model;
%write_latex_static_model;
%write_latex_definitions;
%write_latex_parameter_table;
%collect_latex_files;