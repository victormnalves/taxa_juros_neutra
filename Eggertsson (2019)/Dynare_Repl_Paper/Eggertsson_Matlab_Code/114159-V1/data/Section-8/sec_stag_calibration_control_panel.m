%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Secular Stagnation Control Panel             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% for Eggertsson, Mehrotra, & Robbins, A Model of Secular Stagnation:
% Theory and Quantitative Evaluation

% This program generates all the results we need for the secular stagnation
% paper. It combines calibration with running the results and storing them.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1: Calibration for 2015 Non Sec Stag Economy    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   % This section calibrates the economy to the US economy in 2015. We do
   % this by choosing parameters to minimize an objective function. 

    clear all
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SET UP OPTIONS                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    calibrate_model = 0; % Set to 1 if you want to recalibrate model
    run_transition = 1; % Set to 1 to run with transition results
    run_alt = 0; % Set to 1 to run alternative calibrations
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set up moments to match for 2015             %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % These are the moments we match in the objective function. 
    
        moments_2015.IY = .159;
        moments_2015.r = -.0147;
        moments_2015.debt_inc = .0633;
        moments_2015.ls = 65.99;
        moments_2015.beq_inc = .03;
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set up other params for 2015                 %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        % These are parameters we do not calibrate ourselves through the
        % objective function, but take from other literature. 
        
        other_params_main.gamma = .75;
        other_params_main.deprec = .1244;
        other_params_main.sigma = .6;
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Reminder of how param vector is set up       %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        % The order of the vector we are minimizing matters -- here are
        % where all of the parameters are stored in the vector
    
        %     delta = param_vector(1);
        %     epsilon = param_vector(2);
        %     dlim_pct = param_vector(3);
        %     theta = param_vector(4);
        %     mu = param_vector(5);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set up Minimization                          %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        % Generate lower and upper bounds for minimization procedure

        lb = [1^(-1)-1;.20;.05;3;10];
        ub = [.95060^(-1)-1;.26;.4;12;50];
        
        % Our initial guess
        
        param_guess = [.0296;.2380;.2014;5.1009;30];

        % Minimize the objective function and store it

        if calibrate_model==1 % only run of the code if we are calibrating the model
            
            calibration_main_2015 = fmincon(@(param) sec_stag_calibration(param,other_params_main,moments_2015),param_guess,[],[],[],[],lb,ub);

            save('matlab_results/calibration_main_2015','calibration_main_2015')

        else % else simply load in the pre-stored sults
            
            load('matlab_results/calibration_main_2015.mat') 
        end

        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set up other params for 1970                 %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % We also calibrate a few parameters for the US economy in 1970.
        % Here we give the program some of the parameters we just
        % calibrated. 
        
        other_params_main.delta = calibration_main_2015(1);
        other_params_main.epsilon = calibration_main_2015(2);
        other_params_main.mu = calibration_main_2015(5);
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set up moments to match for 1970             %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % These are the moments we match for 1970
        
        moments_1970.debt_inc = .0421;
        moments_1970.ls = 72.40;

        % Generate lower and upper bounds for minimization

        lb = [.05;2];
        ub = [.5;20];

        param_guess = [.2014;5.1009];

        % Minimize objective function and store results
        
        % Only calibrate if calibrate_model is set to 1. Else, read in
        % pre-calibrated parameters. 
        if calibrate_model==1

            calibration_main_1970 = fmincon(@(param) sec_stag_calibration_1970(param,other_params_main,moments_1970),param_guess,[],[],[],[],lb,ub);

            save('matlab_results/calibration_main_1970','calibration_main_1970')

        else
            
            load('matlab_results/calibration_main_1970.mat')
            
            
        end
        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 2: Run Full Set of Results for Main Calibration     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Now that we have our calibration we can run the model and store the
    % results. 
        
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Use the calibrated parameters                %%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        params_main.gamma = other_params_main.gamma;
        params_main.delta = calibration_main_2015(1); 
        params_main.deprec = other_params_main.deprec;
        params_main.epsilon = calibration_main_2015(2);
        params_main.sigma = other_params_main.sigma;
        params_main.mu = calibration_main_2015(5);
        
        params_main.dlim_pct_1970 = calibration_main_1970(1);
        params_main.dlim_pct_2015 = calibration_main_2015(3);
        
        params_main.theta_1970 = calibration_main_1970(2);
        params_main.theta_2015 = calibration_main_2015(4);
        
        params_main.delta_ss = .9849^(-1) - 1;
        params_main.gamma_as = 0.9045;
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Create Table of Calibration Results                      %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        
        resultsfile = fopen('text_results/main_calib.txt','w');
        
        fprintf(resultsfile,['\n']);
        fprintf(resultsfile,['2015 Main Calibration Results']);
        fprintf(resultsfile,['\n']);
        fprintf(resultsfile,['\n']);
        fprintf(resultsfile,['1.Beta..............,,,,,,,,,,,,,,,,,,..........' num2str(1/(1+params_main.delta),'%04.4f')  ]);
        fprintf(resultsfile,['\n']);
        fprintf(resultsfile,['2.Alpha.........................................' num2str(params_main.epsilon,'%04.4f')  ]);
        fprintf(resultsfile,['\n']);
        fprintf(resultsfile,['3.Mu............................................' num2str(params_main.mu,'%04.4f')  ]);
        fprintf(resultsfile,['\n']);
        fprintf(resultsfile,['4.Theta.........................................' num2str(params_main.theta_2015,'%04.4f')  ]);
        fprintf(resultsfile,['\n']);
        fprintf(resultsfile,['5.Debt limit....................................' num2str(params_main.dlim_pct_2015,'%04.4f')  ]);

        fprintf(resultsfile,['\n']);
        fprintf(resultsfile,['\n']);
        fprintf(resultsfile,['\n']);
        fprintf(resultsfile,['1970 Main Calibration Results']);
        fprintf(resultsfile,['\n']);
        fprintf(resultsfile,['\n']);
        fprintf(resultsfile,['1.Theta.........................................' num2str(params_main.theta_1970,'%04.4f')  ]);
        fprintf(resultsfile,['\n']);
        fprintf(resultsfile,['2.Debt limit....................................' num2str(params_main.dlim_pct_1970,'%04.4f')  ]);

        fclose(resultsfile);
        
        resultsfile = fopen('text_results/moments.txt','w');
        
        fprintf(resultsfile,['\n']);
        fprintf(resultsfile,['\n']);
        fprintf(resultsfile,['2015 Moments to Match']);
        fprintf(resultsfile,['\n']);
        fprintf(resultsfile,['\n']);
        fprintf(resultsfile,['1.Natural Interest Rate........................' num2str(moments_2015.r,'%04.4f')  ]);
        fprintf(resultsfile,['\n']);
        fprintf(resultsfile,['2.Nom Investment Output........................' num2str(moments_2015.IY,'%04.4f')  ]);
        fprintf(resultsfile,['\n']);
        fprintf(resultsfile,['3.Labor Share...................................' num2str(moments_2015.ls,'%04.4f')  ]);
        fprintf(resultsfile,['\n']);
        fprintf(resultsfile,['4.Debt to income................................' num2str(moments_2015.debt_inc,'%04.4f')  ]);
        fprintf(resultsfile,['\n']);
        fprintf(resultsfile,['5.Bequests to income............................' num2str(moments_2015.beq_inc,'%04.4f')  ]);
        fprintf(resultsfile,['\n']);
        
        fprintf(resultsfile,['\n']);
        fprintf(resultsfile,['\n']);
        fprintf(resultsfile,['1970 Moments to Match']);
        fprintf(resultsfile,['\n']);
        fprintf(resultsfile,['\n']);
        fprintf(resultsfile,['1.Natural Interest Rate........................' num2str(2.62,'%04.4f')  ]);
        fprintf(resultsfile,['\n']);
        fprintf(resultsfile,['2.Nom Investment Output........................' num2str(16.8,'%04.4f')  ]);
        fprintf(resultsfile,['\n']);
        fprintf(resultsfile,['3.Labor Share...................................' num2str(moments_1970.ls,'%04.4f')  ]);
        fprintf(resultsfile,['\n']);
        fprintf(resultsfile,['4.Debt to income................................' num2str(moments_1970.debt_inc,'%04.4f')  ]);
        fprintf(resultsfile,['\n']);
        fprintf(resultsfile,['5.Bequests to income............................' num2str(3.0,'%04.4f')  ]);
        fprintf(resultsfile,['\n']);
        
        fclose(resultsfile);
       
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Set up controller, run it                    %%
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
	% The controller controls which results are actually run. Set
	% transition = 0 for a short run. 
    
        controller_main.transition = run_transition; % Set equal to 0 for a shorter run
        controller_main.ff = 0; % Financial frictions
        controller_main.elasticities = 0; 
        controller_main.decomposition = 0;
        controller_main.sec_stag = 0; % Secular stag equilibrium
        
        controller_main.savename = 'main_calibration';
        
        % This runs all of the results, creates charts, graphs, etc. 
        %sec_stag_runmachine(params_main,controller_main)
        [ps_full,gov_full,economy_full,prices_full,nipa_bank_full,results]=sec_stag_runmachine(params_main,controller_main);
        
        save('matlab_transition.mat');

% 
%         
%         controller_test = controller_main;
%         controller_test.elasticities = 0;
%         controller_test.decomposition = 0;
%         controller_test.sec_stag = 0;
% 
%         % sec_stag_runmachine(params_main,controller_test)
%         
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Step 3: Calibration for Alt 1                        %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%     % Now we do the same as the above, except we do different calibrations.
% 
%     % Change depreciation
%     
%         other_params_alt1 = other_params_main;
%         
%         other_params_alt1.deprec = .08;
% 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % Set up Minimization                          %
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%         % Generate lower and upper bounds for minimization
% 
%         lb = [1^(-1)-1;.20;.05;3;10];
%         ub = [.95060^(-1)-1;.30;.4;12;50];
% 
%         param_guess = [.0296;.2380;.2014;5.1009;30];
% 
%         % Minimize 
% 
%         if calibrate_model==1
%             
%             calibration_alt1_2015 = fmincon(@(param) sec_stag_calibration(param,other_params_alt1,moments_2015),param_guess,[],[],[],[],lb,ub);
% 
%             save('matlab_results/calibration_alt1_2015','calibration_alt1_2015')
% 
%         else
%             
%             load('matlab_results/calibration_alt1_2015')
%         
%         end
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % Set up other params for 1970                 %
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         
%         other_params_alt1.delta = calibration_alt1_2015(1);
%         other_params_alt1.epsilon = calibration_alt1_2015(2);
%         other_params_alt1.mu = calibration_alt1_2015(5);
%         
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % Set up moments to match for 1970             %
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%         % Generate lower and upper bounds for minimization
% 
%         lb = [.05;2];
%         ub = [.5;20];
% 
%         param_guess = [.2014;5.1009];
% 
%         % Minimize 
% 
%         if calibrate_model==1
%             
%             calibration_alt1_1970 = fmincon(@(param) sec_stag_calibration_1970(param,other_params_alt1,moments_1970),param_guess,[],[],[],[],lb,ub);
% 
%             save('matlab_results/calibration_alt1_1970','calibration_alt1_1970')
% 
%         else
%             
%             load('matlab_results/calibration_alt1_1970')
%         
%         end
%         
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Step 4: Run Full Set of Results for Alt1                     %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   % Set up the calibrationparameters             %%
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%         params_alt1.gamma = other_params_alt1.gamma;
%         params_alt1.delta = calibration_alt1_2015(1); 
%         params_alt1.deprec = other_params_alt1.deprec;
%         params_alt1.epsilon = calibration_alt1_2015(2);
%         params_alt1.sigma = other_params_alt1.sigma;
%         params_alt1.mu = calibration_alt1_2015(5);
%         
%         params_alt1.dlim_pct_1970 = calibration_alt1_1970(1);
%         params_alt1.dlim_pct_2015 = calibration_alt1_2015(3);
%         
%         params_alt1.theta_1970 = calibration_alt1_1970(2);
%         params_alt1.theta_2015 = calibration_alt1_2015(4);
%         
%         params_alt1.delta_ss = .9895^(-1)-1;
%         params_alt1.gamma_as = .9195;
%         
%         
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % Create Table of Calibration Results                      %
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%         
%         resultsfile = fopen('text_results/alt1_calib.txt','w');
%         
%         fprintf(resultsfile,['\n']);
%         fprintf(resultsfile,['2015 Alt1 Calibration Results']);
%         fprintf(resultsfile,['\n']);
%         fprintf(resultsfile,['\n']);
%         fprintf(resultsfile,['1.Beta..............,,,,,,,,,,,,,,,,,,..........' num2str(1/(1+params_alt1.delta),'%04.4f')  ]);
%         fprintf(resultsfile,['\n']);
%         fprintf(resultsfile,['2.Alpha.........................................' num2str(params_alt1.epsilon,'%04.4f')  ]);
%         fprintf(resultsfile,['\n']);
%         fprintf(resultsfile,['3.Mu............................................' num2str(params_alt1.mu,'%04.4f')  ]);
%         fprintf(resultsfile,['\n']);
%         fprintf(resultsfile,['4.Theta.........................................' num2str(params_alt1.theta_2015,'%04.4f')  ]);
%         fprintf(resultsfile,['\n']);
%         fprintf(resultsfile,['5.Debt limit....................................' num2str(params_alt1.dlim_pct_2015,'%04.4f')  ]);
% 
%         fprintf(resultsfile,['\n']);
%         fprintf(resultsfile,['\n']);
%         fprintf(resultsfile,['\n']);
%         fprintf(resultsfile,['1970 Alt1 Calibration Results']);
%         fprintf(resultsfile,['\n']);
%         fprintf(resultsfile,['\n']);
%         fprintf(resultsfile,['1.Theta.........................................' num2str(params_alt1.theta_1970,'%04.4f')  ]);
%         fprintf(resultsfile,['\n']);
%         fprintf(resultsfile,['2.Debt limit....................................' num2str(params_alt1.dlim_pct_1970,'%04.4f')  ]);
% 
%         fclose(resultsfile);
%         
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   % Set up controller, run it                    %%
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%         controller_alt1 = controller_main;
%         controller_alt1.sec_stag = 0;
%         
%         controller_alt1.savename = 'alt1_calibration';
%         
%         sec_stag_runmachine(params_alt1,controller_alt1)
%         
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Step 5: Calibration for Alt 2                        %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%     % Change gamma to equal 1
%     
%         other_params_alt2 = other_params_main;
%         
%         other_params_alt2.gamma = .99;
% 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % Set up Minimization                          %
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%         % Generate lower and upper bounds for minimization
% 
%         lb = [1^(-1)-1;.20;.05;3;10];
%         ub = [.95060^(-1)-1;.30;.4;12;50];
% 
%         param_guess = [.0296;.2380;.2014;5.1009;30];
% 
%         % Minimize 
% 
%         
%         if calibrate_model==1
%             
%             calibration_alt2_2015 = fmincon(@(param) sec_stag_calibration(param,other_params_alt2,moments_2015),param_guess,[],[],[],[],lb,ub);
% 
%             save('matlab_results/calibration_alt2_2015','calibration_alt2_2015')
% 
%         else
%             
%             load('matlab_results/calibration_alt2_2015')
%             
%         end
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % Set up other params for 1970                 %
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         
%         other_params_alt2.delta = calibration_alt2_2015(1);
%         other_params_alt2.epsilon = calibration_alt2_2015(2);
%         other_params_alt2.mu = calibration_alt2_2015(5);
%         
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % Set up moments to match for 1970             %
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%         % Generate lower and upper bounds for minimization
% 
%         lb = [.05;2];
%         ub = [.5;20];
% 
%         param_guess = [.2014;5.1009];
% 
%         % Minimize 
% 
%         if calibrate_model==1
%             
%             calibration_alt2_1970 = fmincon(@(param) sec_stag_calibration_1970(param,other_params_alt2,moments_1970),param_guess,[],[],[],[],lb,ub);
% 
%             save('matlab_results/calibration_alt2_1970','calibration_alt2_1970')
% 
%         else
%             
%             load('matlab_results/calibration_alt2_1970')
%             
%         end
%         
%         
% % Only run if this parameter is set to 1.     
% if run_alt==1 
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % Step 6: Run Full Set of Results for alt2                     %
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       % Set up the calibrationparameters             %%
%       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%             params_alt2.gamma = other_params_alt2.gamma;
%             params_alt2.delta = calibration_alt2_2015(1); 
%             params_alt2.deprec = other_params_alt2.deprec;
%             params_alt2.epsilon = calibration_alt2_2015(2);
%             params_alt2.sigma = other_params_alt2.sigma;
%             params_alt2.mu = calibration_alt2_2015(5);
% 
%             params_alt2.dlim_pct_1970 = calibration_alt2_1970(1);
%             params_alt2.dlim_pct_2015 = calibration_alt2_2015(3);
% 
%             params_alt2.theta_1970 = calibration_alt2_1970(2);
%             params_alt2.theta_2015 = calibration_alt2_2015(4);
% 
%             params_alt2.delta_ss = .9895^(-1)-1;
%             params_alt2.gamma_as = .9195;
% 
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         % Create Table of Calibration Results                      %
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 
%             resultsfile = fopen('text_results/alt2_calib.txt','w');
% 
%             fprintf(resultsfile,['\n']);
%             fprintf(resultsfile,['2015 Alt2 Calibration Results']);
%             fprintf(resultsfile,['\n']);
%             fprintf(resultsfile,['\n']);
%             fprintf(resultsfile,['1.Beta..............,,,,,,,,,,,,,,,,,,..........' num2str(1/(1+params_alt2.delta),'%04.4f')  ]);
%             fprintf(resultsfile,['\n']);
%             fprintf(resultsfile,['2.Alpha.........................................' num2str(params_alt2.epsilon,'%04.4f')  ]);
%             fprintf(resultsfile,['\n']);
%             fprintf(resultsfile,['3.Mu............................................' num2str(params_alt2.mu,'%04.4f')  ]);
%             fprintf(resultsfile,['\n']);
%             fprintf(resultsfile,['4.Theta.........................................' num2str(params_alt2.theta_2015,'%04.4f')  ]);
%             fprintf(resultsfile,['\n']);
%             fprintf(resultsfile,['5.Debt limit....................................' num2str(params_alt2.dlim_pct_2015,'%04.4f')  ]);
% 
%             fprintf(resultsfile,['\n']);
%             fprintf(resultsfile,['\n']);
%             fprintf(resultsfile,['\n']);
%             fprintf(resultsfile,['1970 Alt2 Calibration Results']);
%             fprintf(resultsfile,['\n']);
%             fprintf(resultsfile,['\n']);
%             fprintf(resultsfile,['1.Theta.........................................' num2str(params_alt2.theta_1970,'%04.4f')  ]);
%             fprintf(resultsfile,['\n']);
%             fprintf(resultsfile,['2.Debt limit....................................' num2str(params_alt2.dlim_pct_1970,'%04.4f')  ]);
% 
%             fclose(resultsfile);
% 
% 
%       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       % Set up controller, run it                    %%
%       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%             controller_alt2 = controller_main;
%             controller_alt2.sec_stag = 0;
% 
%             controller_alt2.savename = 'alt2_calibration';
% 
%             sec_stag_runmachine(params_alt2,controller_alt2)
% 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % Step 7: Calibration for Alt 3                        %
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%         % Change sigma to equal to 1
% 
%             other_params_alt3 = other_params_main;
% 
%             other_params_alt3.sigma = 1;
% 
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         % Set up Minimization                          %
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%             % Generate lower and upper bounds for minimization
% 
%             lb = [1^(-1)-1;.20;.05;3;10];
%             ub = [.95060^(-1)-1;.30;.4;12;50];
% 
%             param_guess = [.0296;.2380;.2014;5.1009;30];
% 
%             % Minimize 
% 
%             if calibrate_model==1
% 
%                 calibration_alt3_2015 = fmincon(@(param) sec_stag_calibration(param,other_params_alt3,moments_2015),param_guess,[],[],[],[],lb,ub);
% 
%                 save('matlab_results/calibration_alt3_2015','calibration_alt3_2015')
% 
%             else
% 
%                 load('matlab_results/calibration_alt3_2015')
% 
%             end
% 
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         % Set up other params for 1970                 %
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%             other_params_alt3.delta = calibration_alt3_2015(1);
%             other_params_alt3.epsilon = calibration_alt3_2015(2);
%             other_params_alt3.mu = calibration_alt3_2015(5);
% 
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         % Set up moments to match for 1970             %
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%             % Generate lower and upper bounds for minimization
% 
%             lb = [.05;2];
%             ub = [.5;20];
% 
%             param_guess = [.2014;5.1009];
% 
%             % Minimize 
% 
%             if calibrate_model==1
% 
%                 calibration_alt3_1970 = fmincon(@(param) sec_stag_calibration_1970(param,other_params_alt3,moments_1970),param_guess,[],[],[],[],lb,ub);
% 
%                 save('matlab_results/calibration_alt3_1970','calibration_alt3_1970')
% 
%             else
% 
%                 load('matlab_results/calibration_alt3_1970')
% 
%             end
% 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % Step 8: Run Full Set of Results for alt3                     %
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       % Set up the calibrationparameters             %%
%       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%             params_alt3.gamma = other_params_alt3.gamma;
%             params_alt3.delta = calibration_alt3_2015(1); 
%             params_alt3.deprec = other_params_alt3.deprec;
%             params_alt3.epsilon = calibration_alt3_2015(2);
%             params_alt3.sigma = other_params_alt3.sigma;
%             params_alt3.mu = calibration_alt3_2015(5);
% 
%             params_alt3.dlim_pct_1970 = calibration_alt3_1970(1);
%             params_alt3.dlim_pct_2015 = calibration_alt3_2015(3);
% 
%             params_alt3.theta_1970 = calibration_alt3_1970(2);
%             params_alt3.theta_2015 = calibration_alt3_2015(4);
% 
%             params_alt3.delta_ss = .9895^(-1)-1;
%             params_alt3.gamma_as = .9195;
% 
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         % Create Table of Calibration Results                      %
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 
%             resultsfile = fopen('text_results/alt3_calib.txt','w');
% 
%             fprintf(resultsfile,['\n']);
%             fprintf(resultsfile,['2015 Alt3 Calibration Results']);
%             fprintf(resultsfile,['\n']);
%             fprintf(resultsfile,['\n']);
%             fprintf(resultsfile,['1.Beta..............,,,,,,,,,,,,,,,,,,..........' num2str(1/(1+params_alt3.delta),'%04.4f')  ]);
%             fprintf(resultsfile,['\n']);
%             fprintf(resultsfile,['2.Alpha.........................................' num2str(params_alt3.epsilon,'%04.4f')  ]);
%             fprintf(resultsfile,['\n']);
%             fprintf(resultsfile,['3.Mu............................................' num2str(params_alt3.mu,'%04.4f')  ]);
%             fprintf(resultsfile,['\n']);
%             fprintf(resultsfile,['4.Theta.........................................' num2str(params_alt3.theta_2015,'%04.4f')  ]);
%             fprintf(resultsfile,['\n']);
%             fprintf(resultsfile,['5.Debt limit....................................' num2str(params_alt3.dlim_pct_2015,'%04.4f')  ]);
% 
%             fprintf(resultsfile,['\n']);
%             fprintf(resultsfile,['\n']);
%             fprintf(resultsfile,['\n']);
%             fprintf(resultsfile,['1970 Alt3 Calibration Results']);
%             fprintf(resultsfile,['\n']);
%             fprintf(resultsfile,['\n']);
%             fprintf(resultsfile,['1.Theta.........................................' num2str(params_alt3.theta_1970,'%04.4f')  ]);
%             fprintf(resultsfile,['\n']);
%             fprintf(resultsfile,['2.Debt limit....................................' num2str(params_alt3.dlim_pct_1970,'%04.4f')  ]);
% 
%             fclose(resultsfile);
% 
%       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       % Set up controller, run it                    %%
%       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%             controller_alt3 = controller_main;
%             controller_alt3.sec_stag = 0;
% 
%             controller_alt3.savename = 'alt3_calibration';
% 
%             sec_stag_runmachine(params_alt3,controller_alt3)
% 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % Step 7: Calibration for Alt 4                        %
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%         % Change gamma to .5
% 
%             other_params_alt4 = other_params_main;
% 
%             other_params_alt4.gamma = .5;
% 
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         % Set up Minimization                          %
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%             % Generate lower and upper bounds for minimization
% 
%             lb = [1^(-1)-1;.20;.05;3;10];
%             ub = [.95060^(-1)-1;.30;.4;12;50];
% 
%             param_guess = [.0296;.2380;.2014;5.1009;30];
% 
%             % Minimize 
% 
%             if calibrate_model==1
% 
%                 calibration_alt4_2015 = fmincon(@(param) sec_stag_calibration(param,other_params_alt4,moments_2015),param_guess,[],[],[],[],lb,ub);
% 
%                 save('matlab_results/calibration_alt4_2015','calibration_alt4_2015')
% 
%             else
% 
%                 load('matlab_results/calibration_alt4_2015')
% 
%             end
% 
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         % Set up other params for 1970                 %
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%             other_params_alt4.delta = calibration_alt4_2015(1);
%             other_params_alt4.epsilon = calibration_alt4_2015(2);
%             other_params_alt4.mu = calibration_alt4_2015(5);
% 
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         % Set up moments to match for 1970             %
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
%             % Generate lower and upper bounds for minimization
% 
%             lb = [.05;2];
%             ub = [.5;20];
% 
%             param_guess = [.2014;5.1009];
% 
%             % Minimize 
% 
% 
%             if calibrate_model==1
% 
%                 calibration_alt4_1970 = fmincon(@(param) sec_stag_calibration_1970(param,other_params_alt4,moments_1970),param_guess,[],[],[],[],lb,ub);
% 
%                 save('matlab_results/calibration_alt4_1970','calibration_alt4_1970')
% 
%             else
% 
%                 load('matlab_results/calibration_alt4_1970')
% 
%             end
% 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % Step 8: Run Full Set of Results for alt4                     %
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       % Set up the calibrationparameters             %%
%       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%             params_alt4.gamma = other_params_alt4.gamma;
%             params_alt4.delta = calibration_alt4_2015(1); 
%             params_alt4.deprec = other_params_alt4.deprec;
%             params_alt4.epsilon = calibration_alt4_2015(2);
%             params_alt4.sigma = other_params_alt4.sigma;
%             params_alt4.mu = calibration_alt4_2015(5);
% 
%             params_alt4.dlim_pct_1970 = calibration_alt4_1970(1);
%             params_alt4.dlim_pct_2015 = calibration_alt4_2015(3);
% 
%             params_alt4.theta_1970 = calibration_alt4_1970(2);
%             params_alt4.theta_2015 = calibration_alt4_2015(4);
% 
%             params_alt4.delta_ss = .9895^(-1)-1;
%             params_alt4.gamma_as = .9195;
% 
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         % Create Table of Calibration Results                      %
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 
%             resultsfile = fopen('text_results/alt4_calib.txt','w');
% 
%             fprintf(resultsfile,['\n']);
%             fprintf(resultsfile,['2015 Alt4 Calibration Results']);
%             fprintf(resultsfile,['\n']);
%             fprintf(resultsfile,['\n']);
%             fprintf(resultsfile,['1.Beta..............,,,,,,,,,,,,,,,,,,..........' num2str(1/(1+params_alt4.delta),'%04.4f')  ]);
%             fprintf(resultsfile,['\n']);
%             fprintf(resultsfile,['2.Alpha.........................................' num2str(params_alt4.epsilon,'%04.4f')  ]);
%             fprintf(resultsfile,['\n']);
%             fprintf(resultsfile,['3.Mu............................................' num2str(params_alt4.mu,'%04.4f')  ]);
%             fprintf(resultsfile,['\n']);
%             fprintf(resultsfile,['4.Theta.........................................' num2str(params_alt4.theta_2015,'%04.4f')  ]);
%             fprintf(resultsfile,['\n']);
%             fprintf(resultsfile,['5.Debt limit....................................' num2str(params_alt4.dlim_pct_2015,'%04.4f')  ]);
% 
%             fprintf(resultsfile,['\n']);
%             fprintf(resultsfile,['\n']);
%             fprintf(resultsfile,['\n']);
%             fprintf(resultsfile,['1970 Alt4 Calibration Results']);
%             fprintf(resultsfile,['\n']);
%             fprintf(resultsfile,['\n']);
%             fprintf(resultsfile,['1.Theta.........................................' num2str(params_alt4.theta_1970,'%04.4f')  ]);
%             fprintf(resultsfile,['\n']);
%             fprintf(resultsfile,['2.Debt limit....................................' num2str(params_alt4.dlim_pct_1970,'%04.4f')  ]);
% 
%             fclose(resultsfile);
% 
%       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       % Set up controller, run it                    %%
%       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%             controller_alt4 = controller_main;
%             controller_alt4.sec_stag = 0;
% 
%             controller_alt4.savename = 'alt4_calibration';
% 
%             sec_stag_runmachine(params_alt4,controller_alt4)
% 
% end
