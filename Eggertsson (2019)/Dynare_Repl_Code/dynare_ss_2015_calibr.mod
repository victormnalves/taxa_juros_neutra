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

%Code used for calibration targets at 2015

%--------------------------------------------------------------------------
//LOAD DATA TO BE USED IN THE MODEL
%--------------------------------------------------------------------------

ss_2015 = xlsread('data/ss_2015_calibr.xlsx');

%Parameters to be Calibrated at 2015 Targets
alpha_par_par = xlsread('data/alpha.xlsx');
beta_par_par = xlsread('data/beta.xlsx');
D_par_par = xlsread('data/D_2015.xlsx');
theta_par_par = xlsread('data/theta_2015.xlsx');
mu_par_par = xlsread('data/mu.xlsx');

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
gov_debt, I, A_adj, IY, LS, BY 

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

%%Parameters that do not change with SS
	@#for j in 1:40
        hc@{j}=datahc(1,@{j});
    @#endfor
    delta=delta_p;
    gamma=gamma_p;
    sigma=sigma_p; 
    g=g_p;     

%%Parameters that change with SS
    @#for j in 1:J
        s@{j}=datasurvival1_2015(1,@{j});
    @#endfor
    @#for j in 1:J
        sv@{j}=datasurvival2_2015(1,@{j});
    @#endfor
    @#for j in 1:J
        su@{j}=data_uncond_survival_2015(1,@{j});
    @#endfor     
    n=n_2015;               
    fert_26=fert_26_2015;                                    
    b=b_2015;               
    AL_growth=AL_growth_2015;
    e=1;       
    AL=1;                              
    AK=1;                               

%To be Calibrated at 2015 Targets
    alpha=alpha_par_par;
    beta=beta_par_par;
    mu=mu_par_par;
    D=D_par_par;
    theta=theta_par_par;

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
    min(lambda@{j},a@{j}+(D*w*hc@{j})*(1+AL_growth)^@{j}) = 0;  
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

[name = 'Bequest to output']
BY=q32*n32/Y;

end;

%--------------------------------------------------------------------------
//STEADY-STATE BLOCK: Initial Values
%--------------------------------------------------------------------------

initval;

@#for j in 1:J
    n@{j}=0.1;
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
A_adj=ss_2015(283);
IY=ss_2015(284);
LS=ss_2015(285);
BY=ss_2015(286);

end;

%--------------------------------------------------------------------------
//STEADY-STATE BLOCK: Computation
%--------------------------------------------------------------------------

%options_.debug=1;
%options_.dynatol.f=1e-3;
steady(maxit=1000, solve_algo=0);
%resid;
%check;

%--------------------------------------------------------------------------
//NEEDED FOR CALIBRATION
%--------------------------------------------------------------------------

AL_growth=0.006460448230438;
capi=oo_.steady_state(113:169,1);
filename='data/ss_2015_calibr.xlsx';
writematrix(oo_.steady_state,filename,'WriteMode','replacefile');

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
global CDY
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