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
%exo_matrix = readmatrix('data/exo_matrix_start10.xlsx');
exo_matrix = readmatrix('data/exo_matrix_start100.xlsx');
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
gov_debt

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

AK, A, A_adj

@#for j in 41:56
	hc@{j}
@#endfor

;

%--------------------------------------------------------------------------
//PARAMETERS' CALIBRATION
%-------------------------------------------------------------------------- 

alpha=0.238808876581195;
beta=0.979996368002075; %beta=1/(1+0.020411945034762)
delta=0.1244;
gamma=0.75;
sigma=0.6; %IF WE SET IT = 1, IT BECOMES A COBB-DOUGLAS AND WE NEED TO CHANGE SOME EQUATIONS
mu=21.629099653217530;
g=0.2128; %government spend as % of gdp

@#for j in 1:56
	hc@{j}=datahc(1,@{j});
@#endfor

AK=1;
A=1;
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
    (1/beta) = ((c@{j+1}(+1)/c@{j})^(-1/gamma))*(1+r(+1)) + (lambda@{j+1})*((c@{j})^(1/gamma))*sv@{j+1}/(su@{j+1}*beta^@{j}*e); 
@#endfor

    [name = 'Households First-order Condition 56']
    x56 = ((fert_26(-55)/mu)^(-gamma))*c56;

//Households: Budget Constraints
%--------------------------------------------------------------------------

[name = 'Households Budget Constraint: initial condition']
a1=0;

@#for j in 1:J-26
    [name = 'Households Budget Constraint @{j}']
    e*a@{j+1} = (rk+(1-delta)*e)*a@{j}(-1)/(sv@{j}) + (pi@{j}+(1-tau)*w*hc@{j}) - c@{j};
@#endfor

    [name = 'Households Budget Constraint 31']
    e*a32 = (rk+(1-delta)*e)*(a31(-1)/sv31+q32(+1))  + (pi31+(1-tau)*w*hc31) - c31;
    
[name = 'Bequest equation']
    q32 = (x56(-1)*fert_26(-56)*n56(-1))/(n32);

@#for j in 32:J-16
    [name = 'Households Budget Constraint @{j}']
    e*a@{j+1} = (rk+(1-delta)*e)*a@{j}(-1)/sv@{j} + (pi@{j}+(1-tau)*w*hc@{j}) - c@{j};
@#endfor

@#for j in 41:J-1
    [name = 'Households Budget Constraint @{j}']
    e*a@{j+1} = (rk+(1-delta)*e)*a@{j}(-1)/sv@{j} - c@{j};
@#endfor
    
    [name = 'Households Budget Constraint 56']
    a57 = (rk+(1-delta)*e)*a56(-1)/sv56 - (fert_26(-55)*x56) - c56;

[name = 'Households Budget Constraint: terminal condition']
a57 = 0;

//Households: Financial Constraints (OBCs)
%--------------------------------------------------------------------------

@#for j in 1:J-16
    [name = 'Households Financial Constraint @{j}']
    min(lambda@{j},a@{j}+D@{j}/(1+r)) = 0;
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
w  = (price*(((alpha)*(AK*K)^((sigma-1)/sigma)+(1-alpha)*(AL*L)^((sigma-1)/sigma))^(1/(sigma-1)))*(1-alpha)*AL^((sigma-1)/sigma)*L^(-1/sigma))/A_adj;
rk = (price*(((alpha)*(AK*K)^((sigma-1)/sigma)+(1-alpha)*(AL*L)^((sigma-1)/sigma))^(1/(sigma-1)))*(alpha)*AK^((sigma-1)/sigma)*K^(-1/sigma))/A_adj;             
Y  = (((alpha)*(AK*K)^((sigma-1)/sigma)+(1-alpha)*(AL*L)^((sigma-1)/sigma))^(sigma/(sigma-1)))/A_adj;
r = rk/e(-1) + (1-delta)*e/e(-1) - 1;
PI = Y/theta;

//Government
%--------------------------------------------------------------------------

[name = 'Government']
b(+1)*Y(+1) = g*Y + (1+r)*b*Y - tau*w*L;
gov_rev=g*Y+r*b*Y;
gov_deficit=(b(+1)*Y(+1)-b*Y)/(gov_rev);
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
    + n@{j}*c@{j}
@#endfor
);

[name = '(K) Kapital']
e*K = (
@#for j in 1:J
    + n@{j}*e*a@{j}(-1)
@#endfor
)-(b*Y);

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

save('dynare_transition.mat');

%--------------------------------------------------------------------------
//PLOTS
%--------------------------------------------------------------------------

% %population of each generation
% pop_indiv_dyn=Simulated_time_series.data(1:152,1:56);
% %consumption of each generation
% c_indiv_dyn=Simulated_time_series.data(1:152,57:112);
% %asset of each generation
% a_indiv_dyn=Simulated_time_series.data(1:152,113:169);
% % %income of each generation
% % %inc_indiv_dyn=;
% %profit of each generation
% profit_indiv_dyn=Simulated_time_series.data(1:152,226:265);
% %br of each generation
% br_indiv_dyn=Simulated_time_series.data(1:152,267);
% % %bg of each generation
% % bg_indiv_dyn=ps_full.opt.bg;
% %bgo of each generation
% bgo_indiv_dyn=Simulated_time_series.data(1:152,266);
% % %bro of each generation
% % bro_indiv_dyn=ps_full.opt.bro;
% % %bg_adj of each generation
% % bg_adj_indiv_dyn=ps_full.opt.bg_adj;
% % %br_adj of each generation
% % br_adj_indiv_dyn=ps_full.opt.br_Adj;
% % %brt of each generation
% % brt_indiv_dyn=ps_full.opt.brt;
% % %bgt of each generation
% % bgt_indiv_dyn=ps_full.opt.bgt;
% % %debtlimit of each generation
% % dlim_indiv_dyn=ps_full.opt.dlim;
% 
% %%Economy Variables economy_full
% 
% % %Labor productivity
% % AL_dyn=economy_full.ag.AL;
% %Capital
% K_dyn=Simulated_time_series.data(1:152,278);
% %Total Labor
% L_dyn=Simulated_time_series.data(1:152,276);
% %Total Income
% Y_dyn=Simulated_time_series.data(1:152,274);
% %Total Consumption
% C_dyn=Simulated_time_series.data(1:152,277);
% %Total Profit
% PI_dyn=Simulated_time_series.data(1:152,273);
% % %Wbase
% % wbase_dyn=economy_full.ag.wbase;
% % %Total brt
% % brt_dyn=economy_full.ag.brt;
% % %Total bgt
% % bgt_dyn=economy_full.ag.bgt;
% % %Total Population
% % pop_dyn=Simulated_time_series.data(1:152,275);
% % %A_adj
% % A_adj_dyn=economy_full.ag.A_adj;
% 
% %%Government Variables gov_dull
% 
% %Gov Debt
% gov_debt_dyn=Simulated_time_series.data(1:152,281);
% %Gov Deficit
% gov_deficit_dyn=Simulated_time_series.data(1:152,280);
% % %Gov Debt/GDP
% % gov_debt_gdp_dyn=gov_full.debt_gdp;
% % %Gov Spending
% % gov_spend_amt_gdp_dyn=gov_full.spend.amt_gdp;
% %Gov Revenues
% gov_tax_revt_dyn=Simulated_time_series.data(1:152,279);
% %Gov Taxes
% gov_tax_dyn=Simulated_time_series.data(1:152,272);
% 
% %%Prices Variables prices_full
% 
% % %Mpk
% % mpk_dyn=prices_full.mpk;
% %Rental k
% rentk_dyn=Simulated_time_series.data(1:152,269);
% %Interest Rate
% r_dyn=Simulated_time_series.data(1:152,270);
% %Wages
% wages_dyn=Simulated_time_series.data(1:152,268);
% 
% %%%%Individual Variables%%%%
% 
% %Population of each generation (n(j))
% figure;
% for j=1:56
%     plot(pop_indiv_dyn(:,j))
%     hold on
% end
% title('n(j)')
% 
% %Consumption of each generation (c(j))
% figure;
% for j=1:56
%     %figure;
%     plot(c_indiv_dyn(:,j))
%     hold on
% end
% title('c(j)')
% 
% %Asset of each generation (a(j))
% figure;
% for j=1:56
%     %figure;
%     plot(a_indiv_dyn(:,j))
%     hold on
% end
% title('a(j)')
% 
% %Profit of each generation (pi(j))
% figure;
% for j=1:40
%     plot(profit_indiv_dyn(:,j))
%     hold on
% end
% title('pi(j)')
% 
% %q32 and ps_full.opt.br(1,32)
% figure;
% plot(br_indiv_dyn);
% legend('Eggertsson')
% title('q32')
% 
% %x56 and ps_full.opt.bgo(1,56)
% figure;
% plot(bgo_indiv_dyn);
% legend('Eggertsson')
% title('x56')
% 
% %%%%Economy Variables%%%%
% 
% %Aggregate Capital (K)
% figure;
% plot(K_dyn);
% title('K')
% %Aggregate Labor (L)
% figure;
% plot(L_dyn);
% title('L')
% %Aggregate Income (Y)
% figure;
% plot(Y_dyn);
% title('Y')
% %Aggregate Consumption (C)
% figure;
% plot(C_dyn);
% title('C')
% %Aggregate Profit (PI)
% figure;
% plot(PI_dyn);
% title('PI')
% 
% %%Government Variables%%%%
% 
% %Government Debt
% figure;
% plot(gov_debt_dyn);
% title('Gov Debt')
% %Government Deficit
% figure;
% plot(gov_deficit_dyn);
% title('Gov Deficit')
% %Government Tax Revenues
% figure;
% plot(gov_tax_revt_dyn);
% title('Gov Tax Revenues')
% 
% %%%%Prices Variables%%%%
% 
% %Rental k (rentk)
% figure;
% plot(rentk_dyn);
% title('Rental k')
% %Interest rate (r)
% figure;
% plot(r_dyn);
% title('Interest Rate (r)')
% %Wage (w)
% figure;
% plot(wages_dyn);
% title('Wage (w)')

%--------------------------------------------------------------------------
//LATEX FILES
%--------------------------------------------------------------------------

%write_latex_dynamic_model;
%write_latex_static_model;
%write_latex_definitions;
%write_latex_parameter_table;
%collect_latex_files;
