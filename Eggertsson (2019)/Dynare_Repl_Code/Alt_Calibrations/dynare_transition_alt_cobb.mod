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

%Code for Transitional Dynamics

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
ss_1970 = readmatrix('data/ss_1970.xlsx');
ss_2015 = readmatrix('data/ss_2015.xlsx');
%exo_matrix = readmatrix('data/exo_matrix_start_alt1.xlsx');
%exo_matrix = readmatrix('data/exo_matrix_start_alt2.xlsx');
%exo_matrix = readmatrix('data/exo_matrix_start_alt3.xlsx');
%exo_matrix = readmatrix('data/exo_matrix_start_alt4.xlsx');
c_scaled_1970 = readmatrix('data/c_scaled_1970.xlsx');
a_scaled_1970 = readmatrix('data/a_scaled_1970.xlsx');
x56_scaled_1970 = readmatrix('data/x56_scaled_1970.xlsx');
q32_scaled_1970 = readmatrix('data/q32_scaled_1970.xlsx');
c_scaled_2015 = readmatrix('data/c_scaled_2015.xlsx');
a_scaled_2015 = readmatrix('data/a_scaled_2015.xlsx');
x56_scaled_2015 = readmatrix('data/x56_scaled_2015.xlsx');
q32_scaled_2015 = readmatrix('data/q32_scaled_2015.xlsx');

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
gov_debt,

@#for j in 1:J-16
	opt_inc@{j}
@#endfor

;

%--------------------------------------------------------------------------
//EXOGENOUS VARIABLES
%--------------------------------------------------------------------------

varexo 

fert_26, e, b, theta, AL,

@#for j in 1:J
    s@{j}
@#endfor

@#for j in 1:J
    sv@{j}
@#endfor

@#for j in 1:J
    su@{j}
@#endfor

@#for j in 1:40
    D@{j}
@#endfor

;

%--------------------------------------------------------------------------
//PARAMETERS
%--------------------------------------------------------------------------

parameters

alpha, beta, delta, gamma, sigma, mu, g

@#for j in 1:J-16
	hc@{j}
@#endfor

AK, A_adj

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

AK=1;
A_adj = ss_1970(283); %IT COMES FROM THE 1970 SS


%--------------------------------------------------------------------------
//MODEL BLOCK: EQUATIONS
%--------------------------------------------------------------------------

model;

//Demographics
%--------------------------------------------------------------------------

[name = 'Demographics Equation 1']
n1 = 1/su25(-1) * n25(-1) * fert_26;

@#for j in 1:J-1
    [name = 'Demographics Equation @{j+1}']
    n@{j+1} = (s@{j}(-1)*n@{j}(-1));
@#endfor

//Households: First-Order conditions
%--------------------------------------------------------------------------

@#for j in 1:J-1
    [name = 'Households First-order Condition @{j}']
    (1/beta) = ((c@{j+1}(+1)/c@{j})^(-1/gamma))*(1+r(+1)) + (lambda@{j+1})*((c@{j})^(1/gamma))*sv@{j+1}/(su@{j}*beta^@{j}); 
@#endfor

    [name = 'Households First-order Condition 56']
    x56 = ((fert_26(-55)/mu)^(-gamma))*c56;

//Households: Budget Constraints
%--------------------------------------------------------------------------

[name = 'Households Budget Constraint: initial condition']
a1=0;

@#for j in 1:J-26
    [name = 'Households Budget Constraint @{j}']
    a@{j+1} = (1+r)/(sv@{j})*a@{j}(-1) + (pi@{j}+(1-tau)*w*hc@{j}) - c@{j};
@#endfor

    [name = 'Households Budget Constraint 31']
    a32 = (1+r)/sv31*a31(-1) + (pi31+(1-tau)*w*hc31) + q32(+1) - c31;
    
[name = 'Bequest equation']
    q32 = (x56(-1)*fert_26(-56)*n56(-1))/(n32);
    %q32 = (x56*fert_26(-55)*n56)/(sv32*n32);

@#for j in 32:J-16
    [name = 'Households Budget Constraint @{j}']
    a@{j+1} = (1+r)/sv@{j}*a@{j}(-1) + (pi@{j}+(1-tau)*w*hc@{j}) - c@{j};
@#endfor

@#for j in 41:J-1
    [name = 'Households Budget Constraint @{j}']
    a@{j+1} = (1+r)/sv@{j}*a@{j}(-1) - c@{j};
@#endfor
    
    [name = 'Households Budget Constraint 56']
    a57 = (1+r)/sv56*a56(-1) - (fert_26(-55)*x56) - c56;

[name = 'Households Budget Constraint: terminal condition']
a57 = 0;

//Households: Financial Constraints (OBCs)
%--------------------------------------------------------------------------

@#for j in 1:J-16
    [name = 'Households Financial Constraint @{j}']
    min(lambda@{j},(a@{j}+((D@{j}(+1)*opt_inc@{j}(+1))))) = 0;
@#endfor
@#for j in 41:J
    [name = 'Households Financial Constraint @{j}']
    min(lambda@{j},(a@{j})) = 0;
@#endfor

@#for j in 1:J-16
    opt_inc@{j} = w*hc@{j};
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
r = rk/e(-1) + (1-delta)*e/e(-1) - 1;

//Government
%--------------------------------------------------------------------------

[name = 'Government']
gov_debt*K=(gov_debt(-1)*K(-1)*(1+r(-1))+g(-1)*Y(-1)-(gov_rev(-1)*(1-gov_deficit(-1))));
gov_deficit = (b(+1)*Y(+1)-gov_debt*K)/gov_rev;
gov_rev = (g*Y+r*gov_debt*K);
tau*w*L = ((gov_rev)*(1-gov_deficit));

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
    + n@{j}*c@{j}
@#endfor
);

[name = '(K) Kapital'] %TRY WITH e AND NOT e(-1)
K = (
@#for j in 1:J
   + n@{j}*a@{j}(-1)/sv@{j}
@#endfor
)/(e(-1)+gov_debt);

end;

%--------------------------------------------------------------------------
//TRANSITIONAL DYNAMICS BLOCK: Initial and Ending Values
%--------------------------------------------------------------------------

initval;

%----------------------------Endog.Variables-------------------------------

@#for j in 1:56
    n@{j}=ss_1970(@{j},1);
@#endfor
 
@#for j in 1:56
    c@{j}=c_scaled_1970(@{j},1);
@#endfor
 
@#for j in 1:57
    a@{j}=a_scaled_1970(@{j},1);
@#endfor
 
@#for j in 170:225
    lambda@{j-169}=ss_1970(@{j},1);
@#endfor
 
@#for j in 226:265
    pi@{j-225}=ss_1970(@{j},1);
@#endfor
 
x56=x56_scaled_1970(1,1);
%q32=ss_1970(267,1);
q32=q32_scaled_1970(1,1);
w=ss_1970(268,1);
rk=ss_1970(269,1);
r=ss_1970(270,1);
price=ss_1970(271,1);
tau=ss_1970(272,1);
PI=ss_1970(273,1);
Y=ss_1970(274,1);
N=ss_1970(275,1);
L=ss_1970(276,1);
C=ss_1970(277,1);
K=ss_1970(278,1);
gov_rev=ss_1970(279,1);
gov_deficit=ss_1970(280,1);
gov_debt=ss_1970(281,1);
@#for j in 1:J-16
    opt_inc@{j} = ss_1970(268,1)*hc@{j};
@#endfor

%----------------------------Exog.Variables--------------------------------

fert_26=exo_matrix(1,1);
e=exo_matrix(1,2);
b=exo_matrix(1,3);
theta=exo_matrix(1,4);
AL=exo_matrix(1,5);

@#for j in 6:61
	s@{j-5}=exo_matrix(1,@{j});
@#endfor

@#for j in 62:117
	sv@{j-61}=exo_matrix(1,@{j});
@#endfor

@#for j in 118:173
	su@{j-117}=exo_matrix(1,@{j});
@#endfor

@#for j in 174:213
	D@{j-173}=exo_matrix(1,@{j});
@#endfor

end;

%options_.debug=1;
%resid;
%steady(maxit=1000, solve_algo=0);
%check;

%--------------------------------------------------------------------------

endval;

%----------------------------Endog.Variables-------------------------------

@#for j in 1:56
    n@{j}=ss_2015(@{j},1);
@#endfor

@#for j in 1:56
    c@{j}=c_scaled_2015(@{j},1);
@#endfor

@#for j in 1:57
    a@{j}=a_scaled_2015(@{j},1);
@#endfor

@#for j in 170:225
    lambda@{j-169}=ss_2015(@{j},1);
@#endfor
 
@#for j in 226:265
    pi@{j-225}=ss_2015(@{j},1);
@#endfor

x56=x56_scaled_2015(1,1);
%q32=ss_2015(267,1);
q32=q32_scaled_2015(1,1);
w=ss_2015(268,1);
rk=ss_2015(269,1);
r=ss_2015(270,1);
price=ss_2015(271,1);
tau=ss_2015(272,1);
PI=ss_2015(273,1);
Y=ss_2015(274,1);
N=ss_2015(275,1);
L=ss_2015(276,1);
C=ss_2015(277,1);
K=ss_2015(278,1);
gov_rev=ss_2015(279,1);
gov_deficit=ss_2015(280,1);
gov_debt=ss_2015(281,1);
@#for j in 1:J-16
    opt_inc@{j} = ss_2015(268,1)*hc@{j};
@#endfor

%----------------------------Exog.Variables--------------------------------
 
fert_26=exo_matrix(152,1);
e=exo_matrix(152,2);
b=exo_matrix(152,3);
theta=exo_matrix(152,4);
AL=exo_matrix(152,5);

@#for j in 6:61
	s@{j-5}=exo_matrix(152,@{j});
@#endfor

@#for j in 62:117
	sv@{j-61}=exo_matrix(152,@{j});
@#endfor

@#for j in 118:173
	su@{j-117}=exo_matrix(152,@{j});
@#endfor

@#for j in 174:213
	D@{j-173}=exo_matrix(152,@{j});
@#endfor
end;
 
%--------------------------------------------------------------------------
//TRANSITIONAL DYNAMICS BLOCK: COMPUTATION
%--------------------------------------------------------------------------

perfect_foresight_setup(periods=150);
oo_.exo_simul=exo_matrix;
perfect_foresight_solver(stack_solve_algo=7, solve_algo=9); 

%--------------------------------------------------------------------------
//FILES PREPERATION FOR COMPARISON WITH EGGERTSSON ET AL. (2019)
%--------------------------------------------------------------------------

save('dynare_transition_alt1.mat');

%--------------------------------------------------------------------------
//LATEX FILES
%--------------------------------------------------------------------------

%write_latex_dynamic_model;
%write_latex_static_model;
%write_latex_definitions;
%write_latex_parameter_table;
%collect_latex_files;
