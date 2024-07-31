function [ps,gov,economy,prices,nipa_bank] = master_control(ps,prod,gov,policy,economy,run_schedule)

% This program controls which equilibria are calculated
% And the algorithm that calculates steady states, transition, etc. 

% Store original kldamp
kldamp_orig = run_schedule.kldamp;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 0 - Calculate Endogenous Quantities                   %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Do fertility as in AK Chapter 11
    
    if run_schedule.fertility_ak==1

        [ ps,gov,economy] = endog_fertility(ps,prod,gov,policy,economy,run_schedule,run_schedule.fertility_year);
        
    end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Schedule 1 - Iterative Algorithm for Initial Steady State  %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if run_schedule.iss.do==1
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Initialize Steady State Variables         %%
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
                
        [ps,gov,prices,economy] = initialize_iss(ps,prod,gov,policy,economy);

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Account for Unemployment                 %%
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        for type=1:economy.num_types

            ps(type).demog.hc = economy.employment * ps(type).demog.hc;

        end
        
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set Time Structure                       %%
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    tm.year_min = 1;
    tm.year_max = 1;
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Begin Algorithm                          %%
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if strcmp(run_schedule.iss.algo,'damp')==1 % If we are doing initial steady state

        keep_running = 1; % when this is set to zero, we stop
        iter = 1; % # of iterations

        while (keep_running==1) && iter < run_schedule.maxiter; % (keep_running==1) && 

            
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Partial Equilibrium Set Everything Up     %%
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

  
        if run_schedule.pe.run == 1

            % In partial equilibrum, we are given a vector of interest rates to
            % start with. Given the vector of interest rates, this
            % determines the capital stock. Given a fixed capital stock, we
            % simply run the program as normal. 

            % From this, we calculate the optimal decisions of individuals

            prices.r(1) = run_schedule.pe.r;

            economy.ag.A_adj = 1;
              
            % Now, given r, calculate capital stock

            guess = 2;
            
            % First, get capital capital labor ratio
            
                options = optimset('Display','off');
                kl_temp = fminunc(@(kl) kl_func(kl,prices.r(1),prod,economy),guess,options);

                kl_func(kl_temp,prices.r(1),prod,economy) ;
                % Now, get capital stock :)

                economy.ag.K(1) = kl_temp*economy.ag.L(1);
            
        end
            
            
            
            
            [ps,gov,economy,prices,keep_running] = equil_update(ps,prod,gov,policy,economy,prices,run_schedule,tm);

            % General equil
            if run_schedule.pe.run == 0
                
                disp(['The iteration is '  num2str(iter)])
                disp(['The aggregate capital stock is '  num2str(economy.ag.K(1))])
                disp(['The aggregate labor supply is '  num2str(economy.ag.L(1))])
            
            % Partial
            else
                
                disp(['The iteration is '  num2str(iter)])
                disp(['The aggregate capital stock is '  num2str(economy.ag.K_supply(1))])
                
            end
           % disp(['The shadow wages are is '])
            
            %p1.opt.sw(1,50:55)
            iter = iter + 1;
            
            % If there are too many iterations, decrease the damp
            if iter == 30
                
                
                run_schedule.kldamp = run_schedule.kldamp / 2;
                
            elseif iter==50
                
                
                run_schedule.kldamp = run_schedule.kldamp / 2;
                
            end
            
            if iter==(run_schedule.maxiter-1)
                
                warning('Nonconvergence!!!! Bad bad bad')
                warning('Nonconvergence!!!! Bad bad bad')
                warning('Nonconvergence!!!! Bad bad bad')
                warning('Nonconvergence!!!! Bad bad bad')
                
            end
            
          
        end
        
        % Double Check my Algorithm
        if any(strcmp('opt_check',fieldnames(run_schedule)))==1 % If the field exists
            
            if run_schedule.opt_check==1
            
                [opt_test] = opt_control_matlab(ps,gov,prod,policy,prices,economy,run_schedule,tm);
        
                check_ratio= (ps.opt.C_iss(1,:)' - opt_test.opt.C_iss(1,:)')./ps.opt.C_iss(1,:)';
                
                if max(abs(check_ratio)) > .01
                    
                    warning('Optimization Algorithm Failed')
                    warning('Optimization Algorithm Failed')
                    warning('Optimization Algorithm Failed')
                    warning('Optimization Algorithm Failed')

                else
                    
                     
                end
                
            end
        end
        % NIPA Statistics
        
            % Create the Matrix to Store Everything
            nipa_bank.a1 = zeros(policy.nt+2,36);
            nipa_bank.a2 = zeros(policy.nt+2,26);
            nipa_bank.personal_savings_rate = zeros(policy.nt+2,1);
            nipa_bank.KY = zeros(policy.nt+2,1);
            nipa_bank.con_pp = zeros(policy.nt+2,1);
            nipa_bank.con_pw = zeros(policy.nt+2,1);
            nipa_bank.Y_pp = zeros(policy.nt+2,1);
            nipa_bank.Y_pw = zeros(policy.nt+2,1);
            nipa_bank.KL = zeros(policy.nt+2,1);
            nipa_bank.KL_adj = zeros(policy.nt+2,1);
            nipa_bank.y_pw = zeros(policy.nt+2,1);
            nipa_bank.golden_rule = zeros(policy.nt+2,1);
            nipa_bank.support_ratio = zeros(policy.nt+2,1);
            nipa_bank.beq_cap = zeros(policy.nt+2,1);
            nipa_bank.beq_inc = zeros(policy.nt+2,1);
            nipa_bank.ss_inc = zeros(policy.nt+2,1);
            nipa_bank.debt_inc = zeros(policy.nt+2,1);
            nipa_bank.debt_cap = zeros(policy.nt+2,1);
            nipa_bank.income_gini = zeros(policy.nt+2,1);
            nipa_bank.wealth_gini = zeros(policy.nt+2,1);
            nipa_bank.KL = zeros(policy.nt+2,1);
            nipa_bank.average_age = zeros(policy.nt+2,1);
            nipa_bank.percent_120 = zeros(policy.nt+2,1);
            nipa_bank.percent_2140 = zeros(policy.nt+2,1);
            nipa_bank.percent_4156 = zeros(policy.nt+2,1);
            nipa_bank.pop_growth = zeros(policy.nt+2,1);
            nipa_bank.inc_share = zeros(policy.nt+2,3);
            nipa_bank.average_return = zeros(policy.nt+2,1);
            nipa_bank.IY = zeros(policy.nt+2,1);
            nipa_bank.IY_nom = zeros(policy.nt+2,1);

        % Calculate NIPA Statistic
        [nipa_bank] = nipa_display_control(nipa_bank,ps,prod,gov,policy,economy,prices,run_schedule,0,1,1);
        
        ak_display(nipa_bank,ps,prod,gov,policy,economy,prices,run_schedule,0,1);

    end % Algorithm
    
    
    % Reset KL Damp, in case it was changed
    run_schedule.kldamp = kldamp_orig;
    
end



% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Schedule 2 - Iterative Algorithm for Final SS              %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if run_schedule.fss.do==1
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Initialize Final Steady State Vars       %%
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
    [ ps,gov,prices,economy ] = initialize_fss(ps,gov,prices,policy,economy);

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set Time Structure                       %%
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    tm.year_min = policy.nt+2;
    tm.year_max = policy.nt+2;
    
    if strcmp(run_schedule.fss.algo,'damp')==1


        keep_running = 1;
        iter = 1;

        while ((keep_running==1) && iter < run_schedule.maxiter) || iter < 3;

            [ps,gov,economy,prices,keep_running] = equil_update(ps,prod,gov,policy,economy,prices,run_schedule,tm);

            disp(['The FSS iteration is '  num2str(iter)])
            disp(['The FSS aggregate capital stock is '  num2str(economy.ag.K(policy.nt+2))])
            disp(['The FSS aggregate labor supply is '  num2str(economy.ag.L(policy.nt+2))])
           % disp(['The shadow wages are is '])
            
            %p1.opt.sw(1,50:55)
            iter = iter + 1;
            
          
        end
        
        % Display
         [nipa_bank] = nipa_display_control(nipa_bank,ps,prod,gov,policy,economy,prices,run_schedule,0,policy.nt+2,policy.nt+2);
         
         ak_display(nipa_bank,ps,prod,gov,policy,economy,prices,run_schedule,151,1);

    end % Algorithm
    
end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Schedule 3 - Iterative Algorithm for Transition            %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if run_schedule.trans.do==1
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Initialize Transition Variables          %%
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
    [ ps,gov,prices,economy ] = initialize_trans(ps,gov,prices,policy,economy);
    
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set Time Structure                       %%
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    tm.year_min = 2;
    tm.year_max = policy.nt+2;
    
    if strcmp(run_schedule.trans.algo,'damp')==1

        keep_running = 1;
        iter = 1;
        
        for yy=2:152
            
            %ps.opt.sw(yy,1:55) = ps.opt.sw(1,1:55) ;
            %ps.opt.a(yy,1:55) = ps.opt.a(1,1:55);
            %ps.opt.l(yy,1:55) = ps.opt.l(1,1:55);
            
            %economy.ag.K(yy) = economy.ag.K(1);
            %economy.ag.L(yy) = economy.ag.L(1);

            %prices.r(yy) = prices.r(1);
            
            %prices.wages(yy) = (1.01)^(yy-1);
        end

        while ((keep_running==1) && iter < run_schedule.maxiter) || iter < 3;

            if iter==12
                
                test = 2;
                
            end
            
             [ps,gov,economy,prices,keep_running] = equil_update(ps,prod,gov,policy,economy,prices,run_schedule,tm);
            
             test = ps.opt.C(2:end,:) ./ ps.opt.C(1:end-1,:);
             test2 = economy.ag.K(2:end) ./ economy.ag.K(1:(end-1));
             test3 = ps.opt.a(2:end,:) ./ ps.opt.a(1:end-1,:);
             test4 = ps.opt.ben(2:end,:) ./ ps.opt.ben(1:end-1,:);
             
             %eagle = ps.opt.a(1,56)*(1+prices.r(1)) + ps.opt.ben(1,56);
            % keep_running = 1;
            disp(['The trans iteration is '  num2str(iter)])
            disp(['The trans aggregate capital stock year 2 is '  num2str(economy.ag.K(2)) ])
            disp(['The trans aggregate capital stock year 50 is '  num2str(economy.ag.K(50)) ])
            disp(['The trans aggregate capital stock year 100 is '  num2str(economy.ag.K(100)) ])
%             disp(['The trans aggregate labor supply is '  num2str(economy.ag.L(2)) ])
%             
%             disp(['The trans ubar is '  num2str(ps.opt.ubar) ])
%             
%             disp(['The trans vdebt is'  num2str(economy.ag.vdebt(50)) '  ' num2str(economy.ag.vdebt(100)) '  ' num2str(economy.ag.vdebt(152)) ])
%             
%             disp(['The V level is'  num2str(ps.opt.V(30)) '  ' num2str(ps.opt.V(60)) '  ' num2str(ps.opt.V(120)) '  ' num2str(ps.opt.V(152)) ])
%             
           % disp(['The shadow wages are is '])
            
            %p1.opt.sw(1,50:55)
            iter = iter + 1;
            
          
        end
        
        [nipa_bank] = nipa_display_control(nipa_bank,ps,prod,gov,policy,economy,prices,run_schedule,0,2,policy.nt+2);
        ak_display_control(nipa_bank,ps,prod,gov,policy,economy,prices,run_schedule,1);

    end % Algorithm
    
end

end

