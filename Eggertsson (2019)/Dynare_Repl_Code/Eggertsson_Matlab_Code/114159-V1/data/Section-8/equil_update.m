function [ps,gov,economy,prices,keep_running] = equil_update(ps,prod,gov,policy,economy,prices,run_schedule,tm)


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1: Update Data structures                             %%
% This summarizes data by creating new variables             %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [ps,gov,economy,prices] = create_profile(ps,prod,gov,policy,economy,prices,run_schedule,tm);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 2: Update Govenrment Variables: Taxes, Debt, Spending %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [ ps,gov] = gov_dec(ps,prod,gov,policy,economy,prices,run_schedule,tm);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 3: Calculate Tobin's Q and Interest Rate              %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     [ps,prices,economy] = calculate_Q(ps,prod,gov,policy,economy,prices,tm);
    
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 4: Calculate Social Security Benefits                 %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
     [ps] = update_ssben(ps,prod,gov,policy,economy,prices,run_schedule,tm);
    
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 5: Calculate Optimal Consumption Bundles              %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     if run_schedule.jr_method==1
         
        %prices.r(1) = 0.0057;
        % ps.opt.br(1,32) = 0;
        [ps] = opt_control_final(ps,gov,prod,policy,prices,economy,run_schedule,tm);
    
     elseif run_schedule.matlab_solve==1
         
         %prices.r(1) = .02;
         % ps.opt.br(1,32) = 0;
         %ps.opt.bgo = zeros(152,56);
         
         [ps] = opt_control_matlab(ps,gov,prod,policy,prices,economy,run_schedule,tm);
         
         %supertest = opt_control_matlab(ps,gov,prod,policy,prices,economy,run_schedule,tm);
         
         %[ps.opt.C_iss(1,:)',supertest.opt.C_iss(1,:)']
         % [ps.opt.a_iss(1,:)',supertest.opt.a_iss(1,:)',ps.opt.dlim(1,:)',[1:56]',ps.opt.C_iss(1,:)',supertest.opt.C_iss(1,:)']
         %  
         
%          for t=1:56
%              
%              ps.opt.dlim(1,t) = ps.opt.dlim(1,t)*(1+economy.ag.AL_growth_iss)^(t-1);
%              
%          end

     else
         %ps.prod_af = ones(ps.demog.lifespan,1);
         [ps_alt] = opt_control(ps,gov,prod,policy,prices,economy,run_schedule,tm.year_min,tm.year_max);
         
     end
     
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 6: Update Social Security Taxes                       %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     [ps,gov,economy] = update_sst(ps,prod,gov,policy,economy,prices,run_schedule,tm);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 7: Update Transfer Authority Amounts                  %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if tm.year_min~=tm.year_max && policy.lsra==1
        
        [ps,economy] = lsra(ps,gov,policy,economy,prices);
    
    end
    

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 8: Update aggregates and determine whether            %%
% Convergence has been reached                               %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

     [economy,keep_running] = repeatfunc(ps,prod,gov,policy,prices,economy,run_schedule,tm);
     
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 9: Update Aggregate Prices                            %%
% Given Aggregate Quantities                                 %%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %if year_min==year_max
        
    [prices,economy] = calculate_prices(prod,gov,economy,policy,prices,tm);
    
    %end

end