%--------------------------------------------------------------------------
% Authors: Alex Crescentini & Federico Giri.
% Università Politecnica delle Marche, Ancona, Italy.
% January 2023.
%--------------------------------------------------------------------------
% Dynare Replication Code for:
% "A Model of Secular Stagnation: Theory and Quantitative Evaluation"
% Eggertsson, Mehrotra and Robbins (2019), American Economic Journal.
%--------------------------------------------------------------------------

%% DESCRIPTION

%The following file compares the output of the model coming from the
%original Matlab code by Eggertsson, Mehrotra and Robbins (2019) and
%by the replication with Dynare from us.

%In detail: First it runs the main calibration and produces the 
%steady state of the model at 1970, at 2015 and the transitional dynamics,
%with both out Dynare code (point A) and the original Matlab code (point B) 
%of the authors. Second (point C), it compares the output with the one 
%obtained by Eggertsson, Mehrotra and Robbins (2019).

%In doing so, you need to change your path and choose if you want to
%recalibrate the model or not.

%% A) DYNARE (our replication code)

clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------------BASELINE CALIBRATION--------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

run_calibration=0; %Choose run_calibration=1 if you want to recalibrate the model, otherwise load the already calibrated parameters

%%Parameters that do not change with SS

    datahc = readmatrix('data/data_hc.xlsx');     %human capital profile
    delta_p=0.1244;                               %depreciation rate
    gamma_p=0.75;                                 %elasticity intert. consumpt.
    sigma_p=0.6;                                  %elasticity K/L: IF WE SET IT = 1, it becomes a Cobb-Douglas and we need to change some equations
    g_p=0.2128;                                   %government spend as % of gdp (see jr_gdebt.xlsx)

%%Parameters that change with SS

    %2015
    datasurvival1_2015 = readmatrix('data/data_survival_1_2015.xlsx');                  %2015 survival probabilities (see gens.m and genmhs.m and jr_survival.csv)
    datasurvival2_2015 = readmatrix('data/data_survival_2_2015.xlsx');                  %2015 survival probabilities (see gens.m and genmhs.m and jr_survival.csv)
    data_uncond_survival_2015 = readmatrix('data/data_uncondit_survival_2015.xlsx');    %2015 survival probabilities (see gens.m and genmhs.m and jr_survival.csv)
    n_2015=-0.002578211534468;                                                          %2015 Population Growth Rate 2015 (see endog_fertility.m)
    fert_26_2015=0.9375;                                                                %2015 Fertility Rate (see endog_fertility.m)
    b_2015=1.172880806665463;                                                           %2015 government debt as % of gdp (see jr_gdebt.xlsx)
    e_2015=1.009368277002316;                                                           %2015 relative price of capital goods
    AL_growth_2015= 0.006460448230438;                                                  %2015 Productivity Growth (see jr_tfp.xlsx)
    AL_2015=2.99049071248668;                                                           
    AK_2015=1;
    %1970
    datasurvival1_1970 = readmatrix('data/data_survival_1_1970.xlsx');                  %1970 survival probabilities (see gens.m and genmhs.m and jr_survival.csv)
    datasurvival2_1970 = readmatrix('data/data_survival_2_1970.xlsx');                  %1970 survival probabilities (see gens.m and genmhs.m and jr_survival.csv)
    data_uncond_survival_1970 = readmatrix('data/data_uncondit_survival_1970.xlsx');    %1970 survival probabilities (see gens.m and genmhs.m and jr_survival.csv)
    n_1970=0.013549868016229;                                                           %1970 Population Growth Rate 2015 (see endog_fertility.m)
    fert_26_1970=1.4;                                                                   %1970 Fertility Rate (see endog_fertility.m)
    b_1970=0.422844989041041;                                                           %1970 government debt as % of gdp (see jr_gdebt.xlsx)
    e_1970=1.3;                                                                         %1970 relative price of capital goods
    AL_growth_1970=0.020209202089143;                                                   %1970 Productivity Growth (see jr_tfp.xlsx)
    AL_1970=1;
    AK_1970=1;
    
if run_calibration==1
    
    %%CALIBRATION for 2015 TARGETS
    
        %Targets: these are the moments we match for 2015
        moments_2015.IY = .159;         %driven by alpha
        moments_2015.r = -.0147;        %driven by beta
        moments_2015.debt_inc = .0633;  %driven by debt limit D
        moments_2015.ls = 0.6599;       %driven by theta
        moments_2015.beq_inc = .03;     %driven by mu
        %Parameters bounds and guesses
        lb = [.20;1/(1+0.051967178624027);.05;3;10];
        ub = [.26;1;.4;12;50];
        param_guess = [.2380;0.971250971250971;.2014;5.1009;30];
        %Run Calibration
        [calibration_2015_baseline_dynare] = fmincon(@(param) calibration_ss_2015(param,moments_2015),param_guess,[],[],[],[],lb,ub);
        %[calibration_2015_baseline_dynare] = fminsearch(@(param) calibration_ss_2015(param,moments_2015),param_guess);
        %Save Results
        writematrix(calibration_2015_baseline_dynare,'calibration_2015_baseline_dynare.xlsx');
        calibrated_parameters_2015=xlsread('calibration_2015_baseline_dynare.xlsx');
        
    %%CALIBRATION for 1970 TARGETS
    
        %Targets: these are the moments we match for 1970
        moments_1970.debt_inc = .0421;  %driven by debt limit D
        moments_1970.ls = 0.7240;       %driven by theta
        %Parameters bounds and guesses
        lb = [.05;2];
        ub = [.5;20];
        param_guess = [.2014;5.1009];
        %Run Calibration
        options = optimoptions('fmincon','Display','iter','DiffMinChange',0.00005,'OptimalityTolerance',1e-10);
        [calibration_1970_baseline_dynare] = fmincon(@(param) calibration_ss_1970(param,moments_1970),param_guess,[],[],[],[],lb,ub,[],options);
        %[calibration_1970_baseline_dynare] = fminsearch(@(param) calibration_ss_1970(param,moments_1970),param_guess);
        %Save Results
        writematrix(calibration_1970_baseline_dynare,'calibration_1970_baseline_dynare.xlsx');
        calibrated_parameters_1970=xlsread('calibration_1970_baseline_dynare.xlsx');
        
else 
    
    calibrated_parameters_1970=xlsread('calibration_1970_baseline_dynare.xlsx');
    calibrated_parameters_2015=xlsread('calibration_2015_baseline_dynare.xlsx');

end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------------------SIMULATION---------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Steady State 1970
dynare dynare_ss_1970;
%Steady State 2015
dynare dynare_ss_2015;
%Transitional Dynamics
dynare dynare_transition nostrict;


%% B) MATLAB (original code from Eggertsson, Mehrotra and Robbins (2019))

clear all
close all

%%%%%%%%%%%%PUT YOUR PATH HERE AND BELOW%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%We have set TOLERANCE (run_schedule.tol) = 1e-5 and MAXITER (run_schedule.maxiter) = 200
%Windows (Federico)
run('.....add your path of this folder here......\EI_ReplicEggertsson2019_CrescentiniGiri_2023\Dynare_Repl_Code\Eggertsson_Matlab_Code\114159-V1\data\Section-8\sec_stag_calibration_control_panel.m');

%macOS (Alex)
%run('/....add your path of this folder here...../EI_ReplicEggertsson2019_CrescentiniGiri_2023/Dynare_Repl_Code/Eggertsson_Matlab_Code/114159-V1/data/Section-8/sec_stag_calibration_control_panel.m');


%% C) LOAD OUTPUT FROM DYNARE AND MATLAB CODES and PRODUCES PLOTS FOR COMPARISON
 
clear all
close all

%dynare output
load dynare_transition.mat 

%%%%%%%%%%%%PUT YOUR PATH HERE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%matlab output
%Windows (Federico)
load '.....add your path of this folder here......\EI_ReplicEggertsson2019_CrescentiniGiri_2023\Dynare_Repl_Code\Eggertsson_Matlab_Code\114159-V1\data\Section-8\matlab_transition.mat'

%macOS (Alex)
%load '/....add your path of this folder here...../EI_ReplicEggertsson2019_CrescentiniGiri_2023/Dynare_Repl_Code/Eggertsson_Matlab_Code/114159-V1/data/Section-8/matlab_transition.mat'


%%Output Comparison

run('plot_comparison.m');