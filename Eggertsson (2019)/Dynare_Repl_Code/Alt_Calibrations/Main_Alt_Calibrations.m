%--------------------------------------------------------------------------
% Authors: Alex Crescentini & Federico Giri.
% Università Politecnica delle Marche, Ancona, Italy.
% January 2023.
%--------------------------------------------------------------------------
% Dynare Replication Code for:
% "A Model of Secular Stagnation: Theory and Quantitative Evaluation"
% Eggertsson, Mehrotra and Robbins (2019), American Economic Journal.
%--------------------------------------------------------------------------

%% Description

%The following file compares the output of the model coming from the
%original Matlab code by Eggertsson, Mehrotra and Robbins (2019) and
%by the replication with Dynare from us for the four alternative 
%calibrations.

%The various output from Eggertsson, Mehrotra and Robbins (2019) have been
%already runned and saved in:
%matlab_transition_alt1.mat,
%matlab_transition_alt2.mat, 
%matlab_transition_alt3.mat,
%matlab_transition_alt4.mat.

%% Choose which alternative calibration you want to compare (1,2,3,4)

clear all
close all

global run_alt_number

run_alt_number = 1; %choose 1,2,3 or, 4

%% Alternative Calibrations

if run_alt_number==1
    
    
    %Matrix of exogenous variables path alt1
    exo_matrix = readmatrix('data/exo_matrix_start_alt1.xlsx');
    
    %Parameters alt1
    [param_nc] = parameters_notchange;
    [param_1970] = parameters_1970;
    [param_2015] = parameters_2015;
    
    %Steady State 1970
    dynare dynare_ss_1970_alt;
    %Steady State 2015
    dynare dynare_ss_2015_alt;
    %Transitional Dynamics
    dynare dynare_transition_alt nostrict;

    clear all
    close all

    %dynare output
    load dynare_transition_alt1.mat 
    %matlab output
    load matlab_transition_alt1.mat %TOLERANCE 10, MAXITER = 100

elseif run_alt_number==2

    %Matrix of exogenous variables path alt2
    exo_matrix = readmatrix('data/exo_matrix_start_alt2.xlsx');
    
    %Parameters alt2
    [param_nc] = parameters_notchange;
    [param_1970] = parameters_1970;
    [param_2015] = parameters_2015;
    
    %Steady State 1970
    dynare dynare_ss_1970_alt;
    %Steady State 2015
    dynare dynare_ss_2015_alt;
    %Transitional Dynamics
    dynare dynare_transition_alt nostrict;

    clear all
    close all

    %dynare output
    load dynare_transition_alt2.mat 
    %matlab output
    load matlab_transition_alt2.mat %TOLERANCE 10, MAXITER = 100
    
elseif run_alt_number==3

    %Matrix of exogenous variables path alt3
    exo_matrix = readmatrix('data/exo_matrix_start_alt3.xlsx');
    
    %Parameters alt3
    [param_nc] = parameters_notchange;
    [param_1970] = parameters_1970;
    [param_2015] = parameters_2015;
    
    %Steady State 1970
    dynare dynare_ss_1970_alt_cobb;
    %Steady State 2015
    dynare dynare_ss_2015_alt_cobb;
    %Transitional Dynamics
    dynare dynare_transition_alt_cobb nostrict;

    clear all
    close all

    %dynare output
    load dynare_transition_alt3.mat 
    %matlab output
    load matlab_transition_alt3.mat %TOLERANCE 10, MAXITER = 100
    
elseif run_alt_number==4

    %Matrix of exogenous variables path alt4
    exo_matrix = readmatrix('data/exo_matrix_start_alt4.xlsx');
    
    %Parameters alt4
    [param_nc] = parameters_notchange;
    [param_1970] = parameters_1970;
    [param_2015] = parameters_2015;
    
    %Steady State 1970
    dynare dynare_ss_1970_alt;
    %Steady State 2015
    dynare dynare_ss_2015_alt;
    %Transitional Dynamics
    dynare dynare_transition_alt nostrict;

    clear all
    close all

    %dynare output
    load dynare_transition_alt4.mat 
    %matlab output
    load matlab_transition_alt4.mat %TOLERANCE 10, MAXITER = 100
    
end

%% ASSIGN VARIABLE: ALTERNATIVE CALIBRATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%%Individual Variables [ps_full in Eggertsson et al. (2019)]%%%%%%%%%%%%%%%

%Population of each generation (n{j})
pop_indiv_dyn=Simulated_time_series.data(1:152,1:56);
pop_indiv_egg=ps_full.demog.pop;
%Consumption of each generation (c{j})
c_indiv_dyn=Simulated_time_series.data(1:152,57:112);
c_indiv_egg=ps_full.opt.C;
%Asset of each generation (a{j})
a_indiv_dyn=Simulated_time_series.data(1:152,113:169);
a_indiv_egg=ps_full.opt.a;
%Profit of each generation (p{j})
profit_indiv_dyn=Simulated_time_series.data(1:152,226:265);
profit_indiv_egg=ps_full.opt.profit;
%Bequest received by generation 32 (q32)
br_indiv_dyn=Simulated_time_series.data(1:152,267);
br_indiv_egg=ps_full.opt.br;
%Bequest given by generation 56 (x56)
bgo_indiv_dyn=Simulated_time_series.data(1:152,266);
bgo_indiv_egg=ps_full.opt.bgo;
%Income of each generation (opt_inc{j})
inc_indiv_dyn=Simulated_time_series.data(1:152,282:321);
inc_indiv_egg=ps_full.opt.inc;

%%Economy Variables [economy_full in Eggertsson et al. (2019)]%%%%%%%%%%%%%

%Total Kapital (K)
K_dyn=Simulated_time_series.data(1:152,278);
K_egg=economy_full.ag.K;
%Total Labor (L)
L_dyn=Simulated_time_series.data(1:152,276);
L_egg=economy_full.ag.L;
%Total Income (Y)
Y_dyn=Simulated_time_series.data(1:152,274);
Y_egg=economy_full.ag.Y;
%Total Consumption (C)
C_dyn=Simulated_time_series.data(1:152,277);
C_egg=economy_full.ag.C;
%Total Profit (PI)
PI_dyn=Simulated_time_series.data(1:152,273);
PI_egg=economy_full.ag.profit;
%Total Population (N)
N_dyn=Simulated_time_series.data(1:152,275);
N_egg=economy_full.ag.pop;

%%Government Variables [gov_full in Eggertsson et al. (2019)]%%%%%%%%%%%%%%

%Government Debt
gov_debt_dyn=Simulated_time_series.data(1:152,281);
gov_debt_egg=gov_full.debt;
%Government Deficit
gov_deficit_dyn=Simulated_time_series.data(1:152,280);
gov_deficit_egg=gov_full.deficit;
%Government Revenues 
gov_tax_revt_dyn=Simulated_time_series.data(1:152,279);
gov_tax_revt_egg=gov_full.tax.revt;
%Government Taxation (tau)
gov_tax_dyn=Simulated_time_series.data(1:152,272);
gov_tax_egg=ps_full.tax.wa(:,1);

%%Prices Variables [prices_full in Eggertsson et al. (2019)]%%%%%%%%%%%%%%%

%Rental k (rk)
rentk_dyn=Simulated_time_series.data(1:152,269);
rentk_egg=prices_full.rentk;
%Interest Rate (r)
r_dyn=Simulated_time_series.data(1:152,270);
r_egg=prices_full.r;
%Wages (w)
wages_dyn=Simulated_time_series.data(1:152,268);
wages_egg=prices_full.wages;

%% PLOT GRAPHS: ALTERNATIVE CALIBRATIONS 

set(0,'defaultfigurecolor',[1 1 1 ])
set( gca                       , ...
    'FontName'   , 'Arial' );

%%%%Individual Variables%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Population of each generation (n(j))
figure;
for j=1:56
    %figure;
    plot(pop_indiv_egg(:,j))
    hold on
    plot(pop_indiv_dyn(:,j))
    hold on
end
title('n(j)')

%Consumption of each generation (c(j))
figure;
for j=1:56
    %figure;
    plot(c_indiv_egg(:,j))
    hold on
    plot(c_indiv_dyn(:,j))
    hold on
end
title('c(j)')

%Asset of each generation (a(j))
figure;
for j=1:56
    %figure;
    plot(a_indiv_egg(:,j))
    hold on
    plot(a_indiv_dyn(:,j))
    hold on
end
title('a(j)')

%Profit of each generation (pi(j))
figure;
for j=1:40
    plot(profit_indiv_egg(:,j))
    hold on
    plot(profit_indiv_dyn(:,j))
    hold on
end
title('pi(j)')

%q32 and ps_full.opt.br(1,32)
figure;
plot(br_indiv_egg(1:152,32),'-b','LineWidth', 2.5);
hold on
plot(br_indiv_dyn(1:152),'--r','LineWidth', 2.5);
%legend('Eggertsson et al. (2019)', 'Dynare')
title('q32 and ps_full.opt.br(1,32)')
%xlim([t(1) t(end)])
xtickangle(45);

set(gca, ...
  'Box'         , 'off'     , ...
  'fontsize'    , 14        , ...
  'FontWeight'  , 'bold'    , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'off'      , ...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'LineWidth'   , 1.0         );

if run_alt_number==1
    saveas(gcf,'Figures/Bequest_alt1','epsc')
elseif run_alt_number==2
    saveas(gcf,'Figures/Bequest_alt2','epsc')
elseif run_alt_number==3
    saveas(gcf,'Figures/Bequest_alt3','epsc')
elseif run_alt_number==4
    saveas(gcf,'Figures/Bequest_alt4','epsc')
end

%x56 and ps_full.opt.bgo(1,56)
figure;
plot(bgo_indiv_egg(:,56),'-b','LineWidth', 2.5);
hold on
plot(bgo_indiv_dyn,'--r','LineWidth', 2.5);
%legend('Eggertsson et al. (2019)', 'Dynare')
title('x56 and ps_full.opt.bgo(1,56)')
%xlim([t(1) t(end)])
xtickangle(45);

set(gca, ...
  'Box'         , 'off'     , ...
  'fontsize'    , 14        , ...
  'FontWeight'  , 'bold'    , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'off'      , ...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'LineWidth'   , 1.0         );

if run_alt_number==1
    saveas(gcf,'Figures/x56_alt1','epsc')
elseif run_alt_number==2
    saveas(gcf,'Figures/x56_alt2','epsc')
elseif run_alt_number==3
    saveas(gcf,'Figures/x56_alt3','epsc')
elseif run_alt_number==4
    saveas(gcf,'Figures/x56_alt4','epsc')
end


%%%%Economy Variables%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Aggregate Capital (K)
figure;
plot(K_egg,'-b','LineWidth', 2.5);
hold on
plot(K_dyn,'--r','LineWidth', 2.5);
%legend('Eggertsson et al. (2019)', 'Dynare')
title('Aggregate Kapital (K)')
%xlim([t(1) t(end)])
xtickangle(45);

set(gca, ...
  'Box'         , 'off'     , ...
  'fontsize'    , 14        , ...
  'FontWeight'  , 'bold'    , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'off'      , ...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'LineWidth'   , 1.0         );

if run_alt_number==1
    saveas(gcf,'Figures/Aggregate_Kapital_alt1','epsc')
elseif run_alt_number==2
    saveas(gcf,'Figures/Aggregate_Kapital_alt2','epsc')
elseif run_alt_number==3
    saveas(gcf,'Figures/Aggregate_Kapital_alt3','epsc')
elseif run_alt_number==4
    saveas(gcf,'Figures/Aggregate_Kapital_alt4','epsc')
end


%Aggregate Labor (L)
figure;
plot(L_egg,'-b','LineWidth', 2.5);
hold on
plot(L_dyn,'--r','LineWidth', 2.5);
%legend('Eggertsson et al. (2019)', 'Dynare')
title('Aggregate Labor (L)')
%xlim([t(1) t(end)])
xtickangle(45);

set(gca, ...
  'Box'         , 'off'     , ...
  'fontsize'    , 14        , ...
  'FontWeight'  , 'bold'    , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'off'      , ...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'LineWidth'   , 1.0         );

if run_alt_number==1
    saveas(gcf,'Figures/Aggregate_Labor_alt1','epsc')
elseif run_alt_number==2
    saveas(gcf,'Figures/Aggregate_Labor_alt2','epsc')
elseif run_alt_number==3
    saveas(gcf,'Figures/Aggregate_Labor_alt3','epsc')
elseif run_alt_number==4
    saveas(gcf,'Figures/Aggregate_Labor_alt4','epsc')
end


%Aggregate Income (Y)
figure;
plot(Y_egg,'-b','LineWidth', 2.5);
hold on
plot(Y_dyn,'--r','LineWidth', 2.5);
%legend('Eggertsson et al. (2019)', 'Dynare')
title('Aggregate Income (Y)')
%xlim([t(1) t(end)])
xtickangle(45);

set(gca, ...
  'Box'         , 'off'     , ...
  'fontsize'    , 14        , ...
  'FontWeight'  , 'bold'    , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'off'      , ...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'LineWidth'   , 1.0         );


if run_alt_number==1
    saveas(gcf,'Figures/Aggregate_Income_alt1','epsc')
elseif run_alt_number==2
    saveas(gcf,'Figures/Aggregate_Income_alt2','epsc')
elseif run_alt_number==3
    saveas(gcf,'Figures/Aggregate_Income_alt3','epsc')
elseif run_alt_number==4
    saveas(gcf,'Figures/Aggregate_Income_alt4','epsc')
end


%Aggregate Consumption (C)
figure;
plot(C_egg,'-b','LineWidth', 2.5);
hold on
plot(C_dyn,'--r','LineWidth', 2.5);
%legend('Eggertsson et al. (2019)', 'Dynare')
title('Aggregate Consumption (C)')
%xlim([t(1) t(end)])
xtickangle(45);

set(gca, ...
  'Box'         , 'off'     , ...
  'fontsize'    , 14        , ...
  'FontWeight'  , 'bold'    , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'off'      , ...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'LineWidth'   , 1.0         );

if run_alt_number==1
    saveas(gcf,'Figures/Aggregate_Consumption_alt1','epsc')
elseif run_alt_number==2
    saveas(gcf,'Figures/Aggregate_Consumption_alt2','epsc')
elseif run_alt_number==3
    saveas(gcf,'Figures/Aggregate_Consumption_alt3','epsc')
elseif run_alt_number==4
    saveas(gcf,'Figures/Aggregate_Consumption_alt4','epsc')
end


%Aggregate Profit (PI)
figure;
plot(PI_egg,'-b','LineWidth', 2.5);
hold on
plot(PI_dyn,'--r','LineWidth', 2.5);
%legend('Eggertsson et al. (2019)', 'Dynare')
title('Aggregate Profit (PI)')
%xlim([t(1) t(end)])
xtickangle(45);

set(gca, ...
  'Box'         , 'off'     , ...
  'fontsize'    , 14        , ...
  'FontWeight'  , 'bold'    , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'off'      , ...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'LineWidth'   , 1.0         );

if run_alt_number==1
    saveas(gcf,'Figures/Aggregate_Profit_alt1','epsc')
elseif run_alt_number==2
    saveas(gcf,'Figures/Aggregate_Profit_alt2','epsc')
elseif run_alt_number==3
    saveas(gcf,'Figures/Aggregate_Profit_alt3','epsc')
elseif run_alt_number==4
    saveas(gcf,'Figures/Aggregate_Profit_alt4','epsc')
end


%Aggregate Population (N)
figure;
plot(N_egg,'-b','LineWidth', 2.5);
hold on
plot(N_dyn,'--r','LineWidth', 2.5);
%legend('Eggertsson et al. (2019)', 'Dynare')
title('Aggregate Population (N)')
%xlim([t(1) t(end)])
xtickangle(45);

set(gca, ...
  'Box'         , 'off'     , ...
  'fontsize'    , 14        , ...
  'FontWeight'  , 'bold'    , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'off'      , ...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'LineWidth'   , 1.0         );

if run_alt_number==1
    saveas(gcf,'Figures/Aggregate_Population_alt1','epsc')
elseif run_alt_number==2
    saveas(gcf,'Figures/Aggregate_Population_alt2','epsc')
elseif run_alt_number==3
    saveas(gcf,'Figures/Aggregate_Population_alt3','epsc')
elseif run_alt_number==4
    saveas(gcf,'Figures/Aggregate_Population_alt4','epsc')
end


%%Government Variables%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Government Debt
figure;
plot(gov_debt_egg,'-b','LineWidth', 2.5);
hold on
plot(gov_debt_dyn,'--r','LineWidth', 2.5);
%legend('Eggertsson et al. (2019)', 'Dynare')
title('Government Debt')
%xlim([t(1) t(end)])
xtickangle(45);

set(gca, ...
  'Box'         , 'off'     , ...
  'fontsize'    , 14        , ...
  'FontWeight'  , 'bold'    , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'off'      , ...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'LineWidth'   , 1.0         );

if run_alt_number==1
    saveas(gcf,'Figures/Gov_debt_alt1','epsc')
elseif run_alt_number==2
    saveas(gcf,'Figures/Gov_debt_alt2','epsc')
elseif run_alt_number==3
    saveas(gcf,'Figures/Gov_debt_alt3','epsc')
elseif run_alt_number==4
    saveas(gcf,'Figures/Gov_debt_alt4','epsc')
end


%Government Deficit
figure;
plot(gov_deficit_egg,'-b','LineWidth', 2.5);
hold on
plot(gov_deficit_dyn,'--r','LineWidth', 2.5);
%legend('Eggertsson et al. (2019)', 'Dynare')
title('Government Deficit')
%xlim([t(1) t(end)])
xtickangle(45);

set(gca, ...
  'Box'         , 'off'     , ...
  'fontsize'    , 14        , ...
  'FontWeight'  , 'bold'    , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'off'      , ...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'LineWidth'   , 1.0         );

if run_alt_number==1
    saveas(gcf,'Figures/Gov_deficit_alt1','epsc')
elseif run_alt_number==2
    saveas(gcf,'Figures/Gov_deficit_alt2','epsc')
elseif run_alt_number==3
    saveas(gcf,'Figures/Gov_deficit_alt3','epsc')
elseif run_alt_number==4
    saveas(gcf,'Figures/Gov_deficit_alt4','epsc')
end


%Government Revenues
figure;
plot(gov_tax_revt_egg,'-b','LineWidth', 2.5);
hold on
plot(gov_tax_revt_dyn,'--r','LineWidth', 2.5);
%legend('Eggertsson et al. (2019)', 'Dynare')
title('Government Tax Revenues')
%xlim([t(1) t(end)])
xtickangle(45);

set(gca, ...
  'Box'         , 'off'     , ...
  'fontsize'    , 14        , ...
  'FontWeight'  , 'bold'    , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'off'      , ...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'LineWidth'   , 1.0         );

if run_alt_number==1
    saveas(gcf,'Figures/Gov_revenues_alt1','epsc')
elseif run_alt_number==2
    saveas(gcf,'Figures/Gov_revenues_alt2','epsc')
elseif run_alt_number==3
    saveas(gcf,'Figures/Gov_revenues_alt3','epsc')
elseif run_alt_number==4
    saveas(gcf,'Figures/Gov_revenues_alt4','epsc')
end


%Government Taxation (tau)
figure;
plot(gov_tax_egg,'-b','LineWidth', 2.5);
hold on
plot(gov_tax_dyn,'--r','LineWidth', 2.5);
%legend('Eggertsson et al. (2019)', 'Dynare')
title('Government Taxation (tau)')
%xlim([t(1) t(end)])
xtickangle(45);

set(gca, ...
  'Box'         , 'off'     , ...
  'fontsize'    , 14        , ...
  'FontWeight'  , 'bold'    , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'off'      , ...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'LineWidth'   , 1.0         );

if run_alt_number==1
    saveas(gcf,'Figures/Gov_Tax_alt1','epsc')
elseif run_alt_number==2
    saveas(gcf,'Figures/Gov_Tax_alt2','epsc')
elseif run_alt_number==3
    saveas(gcf,'Figures/Gov_Tax_alt3','epsc')
elseif run_alt_number==4
    saveas(gcf,'Figures/Gov_Tax_alt4','epsc')
end


%%%%Prices Variables%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Rental k (rentk)
figure;
plot(rentk_egg,'-b','LineWidth', 2.5);
hold on
plot(rentk_dyn,'--r','LineWidth', 2.5);
legend('Eggertsson et al. (2019)', 'Dynare')
title('Rental k (rk)')
%xlim([t(1) t(end)])
xtickangle(45);

set(gca, ...
  'Box'         , 'off'     , ...
  'fontsize'    , 14        , ...
  'FontWeight'  , 'bold'    , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'off'      , ...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'LineWidth'   , 1.0         );

if run_alt_number==1
    saveas(gcf,'Figures/rk_alt1','epsc')
elseif run_alt_number==2
    saveas(gcf,'Figures/rk_alt2','epsc')
elseif run_alt_number==3
    saveas(gcf,'Figures/rk_alt3','epsc')
elseif run_alt_number==4
    saveas(gcf,'Figures/rk_alt4','epsc')
end


%Interest rate (r)
figure;
plot(r_egg,'-b','LineWidth', 2.5);
hold on
plot(r_dyn,'--r','LineWidth', 2.5);
%legend('Eggertsson et al. (2019)', 'Dynare')
title('Interest Rate (r)')
%xlim([t(1) t(end)])
xtickangle(45);

set(gca, ...
  'Box'         , 'off'     , ...
  'fontsize'    , 14        , ...
  'FontWeight'  , 'bold'    , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'off'      , ...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'LineWidth'   , 1.0         );

if run_alt_number==1
    saveas(gcf,'Figures/r_alt1','epsc')
elseif run_alt_number==2
    saveas(gcf,'Figures/r_alt2','epsc')
elseif run_alt_number==3
    saveas(gcf,'Figures/r_alt3','epsc')
elseif run_alt_number==4
    saveas(gcf,'Figures/r_alt4','epsc')
end


%Wage (w)
figure;
plot(wages_egg,'-b','LineWidth', 2.5);
hold on
plot(wages_dyn,'--r','LineWidth', 2.5);
%legend('Eggertsson et al. (2019)', 'Dynare')
title('Wage (w)')
%xlim([t(1) t(end)])
xtickangle(45);

set(gca, ...
  'Box'         , 'off'     , ...
  'fontsize'    , 14        , ...
  'FontWeight'  , 'bold'    , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'off'      , ...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'LineWidth'   , 1.0         );

if run_alt_number==1
    saveas(gcf,'Figures/wage_alt1','epsc')
elseif run_alt_number==2
    saveas(gcf,'Figures/wage_alt2','epsc')
elseif run_alt_number==3
    saveas(gcf,'Figures/wage_alt3','epsc')
elseif run_alt_number==4
    saveas(gcf,'Figures/wage_alt4','epsc')
end

%% HP Filtered Transition Path of the Natural Rate of Interest 

r_egg_cut=r_egg(1:61);
r_dyn_cut=r_dyn(1:61);

[Trend_r_egg,Cyclical_r_egg] = hpfilter(r_egg_cut, 6.25);
%[Trend_r_egg,Cyclical_r_egg] = one_sided_hp_filter_kalman(r_egg_cut, 6.25);
[Trend_r_dyn,Cyclical_r_dyn] = hpfilter(r_dyn_cut, 6.25);
%[Trend_r_dyn,Cyclical_r_dyn] = one_sided_hp_filter_kalman(r_dyn_cut, 6.25);

%HP Trend Interest rate (r)
year_vec = [1:61]' + 1969;

figure;
plot(year_vec, Trend_r_egg,'-b','LineWidth',2);
%plot(year_vec, r_egg_cut,'LineWidth',2);

hold on
plot(year_vec, Trend_r_dyn,'--r','LineWidth',2);
%plot(year_vec, r_dyn_cut,'LineWidth',2);

hold on;
yline(0); hold on;
%legend('Eggertsson et al. (2019)', 'Dynare');
xlabel('Period');
ylabel('Rates');
title('HP Trend Interest Rate % (r)');
legend('Eggertsson et al. (2019)', 'Dynare')

xtickangle(45);

set(gca, ...
  'Box'         , 'off'     , ...
  'fontsize'    , 14        , ...
  'FontWeight'  , 'bold'    , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'YGrid'       , 'off'      , ...
  'XColor'      , [.3 .3 .3], ...
  'YColor'      , [.3 .3 .3], ...
  'LineWidth'   , 1.0         );

if run_alt_number==1
    saveas(gcf,'Figures/int_rate_hp_alt1','epsc')
elseif run_alt_number==2
    saveas(gcf,'Figures/int_rate_hp_alt2','epsc')
elseif run_alt_number==3
    saveas(gcf,'Figures/int_rate_hp_alt3','epsc')
elseif run_alt_number==4
    saveas(gcf,'Figures/int_rate_hp_alt4','epsc')
end
