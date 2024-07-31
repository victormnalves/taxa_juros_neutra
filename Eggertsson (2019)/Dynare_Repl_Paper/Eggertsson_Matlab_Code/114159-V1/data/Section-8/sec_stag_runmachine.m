%function [ results ] = sec_stag_runmachine(params,controller)
function [ps_full,gov_full,economy_full,prices_full,nipa_bank_full,results] = sec_stag_runmachine(params,controller)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sec Stag Runmachine                                      %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Given a structure of parameters, calculates all of the tables we
% could possibly need or possibly want 

% Controller controles what results we actually want to run. Sometimes we
% want to do the transitions, sometimes we don't need to do them 

    results = 1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Open Results File                                %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Print results to a file
        resultsfile = fopen(strcat('text_results/',controller.savename,'_results.txt'),'w');
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set up the calibrationparameters             %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        gamma = params.gamma;
        delta = params.delta;
        deprec = params.deprec;
        epsilon = params.epsilon;
        sigma = params.sigma;
        gamma_as = params.gamma_as;
        
        % Debt limits change to match the moments
        dlim_pct_1970 = params.dlim_pct_1970;
        dlim_pct_2015 = params.dlim_pct_2015;

        % Bequest Parameters Change to match the moments
        
        mu = params.mu;
        
        % Theta parameter changes to match the moments
        theta_1970 = params.theta_1970;
        theta_2015 = params.theta_2015;

    % Number of transition periods
    
    policy.nt = 150;
    
    % Whether there is a lsra
    
    policy.lsra = 0;
    
    % year of transition
    
    policy.year_implemented = 2;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set up the Common Person Structure  %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % While there may be several people, they will share a common
    % structure
    
    % p1 stands for person 1 
    
        lifespan = 56;
        p1.demog.lifespan = lifespan;
        
        % Demographic characteristics. 
        % Need to set them for the initial steady state, the transition
        % years, and the final steady state
        
        p1.endogenous_labor = 0;
        
        % Dlim = 1 means that debt is a percentage of income
        % dlim = 2 means there is a set debt limit
        p1.dlim = 1;
               
        % The age of the parents when they have kids
        p1.demog.age_parent = 26*ones(policy.nt+2,1);
        
        % Set the number of kids that each parent has, to achieve a certain
        % level of population growth
        % In order for there to be a population growth rate of 1.01, each
        % parent must have 1.01^26 kids
        
        % Set years of retirement and human capital profile
        
        years_retirement = 16;
        hc = [18690.96 ; 19163.94 ; 19621.11 ; 20064.87 ; 20497.21 ; 20919.71 ; 21333.54 ; 21739.35 ; 22137.36 ; 22527.25 ; 22908.31 ; 23279.32 ; 23638.63 ; 23984.27 ; 24313.83 ; 24624.66 ; 24913.86 ; 25178.37 ; 25415.05 ; 25620.79 ; 25792.54 ; 25927.47 ; 26023.09 ; 26077.3 ; 26088.44 ; 26055.6 ; 25978.38 ; 25857.24 ; 25693.43 ; 25488.99 ; 25246.82 ; 24970.68 ; 24665.15 ; 24335.61 ; 23988.18 ; 23629.62 ; 23267.47 ; 22909.67 ; 22564.94 ; 22242.54];
        hc = [hc;0*ones(years_retirement,1)];
        hc = hc./18690.96;
        p1.demog.hc = repmat(hc',policy.nt + 2,1) .* ones(policy.nt+2,lifespan); % Human capital profile
        
        % Utility Parameters
        % Need to set them for the initial SS, transition years, and the
        % final steady state
       
        p1.util.gamma = gamma*ones(policy.nt+2,1); % Elasticity of intertemporal substitution
        p1.util.rho = .8*ones(policy.nt+2,1); % Intratemporal elasticity of substitution. Only necessary when there is endogenous labor supply.
        p1.util.delta = delta*ones(policy.nt+2,1); % Time discount rate
        p1.util.alpha = 1.5*ones(policy.nt+2,1); % Leisure preference parameter
        p1.util.a_bg = 56*ones(policy.nt+2,1); % year bequest is given
        
        % Kids always receive bequest at age a_bg (of their parents) - age_parent + 2
        % This is because the way that bequests work: when a parent gives
        % the bequest at 56, the child is of age 31. But the child doesn't
        % receive the bequest until the next year, so he is 32. 
        
        p1.util.a_br = 32*ones(policy.nt+2,1); % year bequest received. = 56 - 26 + 2
        
        p1.util.mu = mu*ones(policy.nt+2,1); % utility parameter for bequests
        
        % Debt limit vector will be created later
        
        p1.util.dlim = 0*ones(policy.nt + 2,lifespan); % gross debt limit, not as a percent of income. Only valid when p1.dlim = 2. Notice this stores the parameter in year age fashion, thus it can vary by age.
        
        % Annuity Markets
        % Participation = 1 means that all lending is done through annuity
        % markets
        % Participation < 1 means there are accidental bequests when
        % individuals die. 
        p1.util.annuity_participation = 1*ones(policy.nt + 2,lifespan);
        
        % Production Parameters
        prod.sigma = sigma*ones(policy.nt+2,1); % Production CES
        prod.epsilon = epsilon*ones(policy.nt+2,1); % Capital share
        prod.b = 0*ones(policy.nt+2,1); % Adjustment costs
        prod.deprec = deprec*ones(policy.nt+2,1); % Depreciaiton
        
        % Government Behavior
        % The government set taxes, social security replacement rate, etc
        
            gov.tax.yp = 0*ones(policy.nt+2,1); % proportional income tax
            gov.tax.yp(1) = 0;
            gov.tax.yg = 0*ones(policy.nt+2,1); % progressive income tax

            gov.tax.wp = 0*ones(policy.nt+2,1); % wage tax
            gov.tax.wp(1) = 0;
            gov.tax.wg = 0*ones(policy.nt+2,1);

            gov.tax.kp = 0*ones(policy.nt+2,1); % capital income tax
            gov.tax.kg = 0*ones(policy.nt+2,1);
            
         % Whether taxes are endogenous or not; 0 is exogenous, 1 is
         % endogenous
         
            gov.tax.ype = 0*ones(policy.nt+2,1); % proportional income tax
            gov.tax.ype(1) = 0;
            gov.tax.yge = 0*ones(policy.nt+2,1); % progressive income tax

            gov.tax.wpe = 1*ones(policy.nt+2,1); % wage tax
            gov.tax.wpe(1) = 1;
            gov.tax.wge = 0*ones(policy.nt+2,1);

            gov.tax.kpe = 0*ones(policy.nt+2,1); % capital income tax
            gov.tax.kge = 0*ones(policy.nt+2,1);

            gov.ss_age = 41*ones(policy.nt+2,1); % age individuals start receiving ss
            gov.rep = .0*ones(policy.nt+2,1); % social security replacement rate

            gov.z = 0*ones(policy.nt+2,1); % Investment tax credits

            % Government Spending
            
            gov.spend.ak = 0; % Do government spending as ak does
            
            % Government Spending as JR Intended 
                % We here set government spending as a % of GDP
                % We also set the Debt as % of GDP (will come later)
            
                gov.spend.jr = 1;
                
                gov.spend.amt = ones(policy.nt+2,1); % Government spending. Simply initializing the variable.
                
                gov.debt = 0*ones(policy.nt+2,1); % Debt as a percentage of the capital stock. 
                gov.debt_alt = 0*ones(policy.nt+2,1); % Gross debt. 
                gov.deficit = 0*ones(policy.nt+2,1); % Deficit. Since spend.jr=1, this doesn't matter
                
            % Proportions for Endogenous Taxes
            
            gov.tax.yrevpro = ones(policy.nt+2,1);
            gov.tax.rrevpro = ones(policy.nt+2,1);
            gov.tax.wrevpro = ones(policy.nt+2,1);
            
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set up the economy structure  %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        economy.num_types = 1; % Number of types for inequality
        economy.type_share = 1; % Types of 
        economy.initial_pop = 1; % Initial population
        economy.employment = 1; % Employment rate for secular stagnation
        economy.ff = .00*ones(policy.nt+2,1); % Financial Frictions
        
        economy.adj_prod = 1; % Whether we adjust productivity so that wages in the first year are zero
        
        % Add net foreign assets holdings, as a percent of the capital
        % stock
        % Later we will also add them as % of GDP
        
        economy.nfa = 0*ones(policy.nt + 2,1);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set up the economy structure  %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        % Hicks Neutral Growth
        economy.ag.A_growth = 0*ones(policy.nt+2,1);
        
        % Capital Productivity Growth
        economy.ag.AK_growth = 0*ones(policy.nt+2,1);
        
        % Set the capital prod, total prod
        
        economy.ag.A = ones(policy.nt+2,1);
        economy.ag.AK = ones(policy.nt+2,1);

        for t=2:(policy.nt + 2)
            
            economy.ag.A(t) = economy.ag.A(t-1)*(1+economy.ag.A_growth(t));
            
            economy.ag.AK(t) = economy.ag.AK(t-1)*(1+economy.ag.AK_growth(t));
            
        end
        
        % Steady State Productivity Growth
        
        economy.ag.A_growth_iss = 0;
        economy.ag.A_growth_fss = 0;
        
        economy.ag.AK_growth_iss = 0;
        economy.ag.AK_growth_fss = 0;
        
        economy.ag.A_growth_adj_iss = 0;
        economy.ag.A_growth_adj_fss = 0;

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set up the baby boom                                     %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Import UN Fertility data
            
                bb_data =  readtable('data/jr_tfr.xlsx');

                tfr_ma = table2array(bb_data(:,'tfr_ma'));

            % This sets the TFR for the dates
            
                bb_1945 = tfr_ma(1);
                bb_2015 = tfr_ma(end);

            % Set up transition vector
            
                num_bb_transition = bb_2015*ones(policy.nt+2,1);
                
                num_bb_transition(1:length(tfr_ma)) = tfr_ma;

            
                num_bb_1945 = bb_1945*ones(policy.nt+2,1);

                num_bb_2015 = bb_2015*ones(policy.nt+2,1);

            % Set up economy.n

                bb_transition = num_bb_transition.^(1./(p1.demog.age_parent-1)) - 1; % Population growth. Sometimes it is easier to specify it like this

            % Initial Stedey State

                bb_vec_1945 = (bb_1945*ones(policy.nt+2,1)).^(1./(p1.demog.age_parent-1)) - 1;
                bb_vec_2015 = (bb_2015*ones(policy.nt+2,1)).^(1./(p1.demog.age_parent-1)) - 1;

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set up Productivity Master Vector                        %
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Transition                                               %
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Import Fernald HP Trended Series
            
            productivity_data =  readtable('data/jr_tfp.xlsx');
            
            fernald_hp = table2array(productivity_data(:,'fernald_hp'));
            
            % Save Final Values
            tpg_1970 = fernald_hp(1);
            tpg_2015 = fernald_hp(end);

            % Create transition vector
            AL_growth_transition = tpg_2015*ones(policy.nt+2,1);
            AL_growth_transition(1:length(fernald_hp)) = fernald_hp;

            % Set total productivity

            AL_transition = ones(policy.nt+2,1);
            for t=2:(policy.nt + 2)

                AL_transition(t) = AL_transition(t-1)*(1+AL_growth_transition(t));

            end

            % Steady State Productivity Growth

            AL_growth_iss_transition = tpg_1970;
            AL_growth_fss_transition = tpg_2015;
            
           
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 1970                                                     %
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            
            AL_growth_1970 = tpg_1970*ones(policy.nt+2,1);

            % Set total productivity

            AL_1970 = ones(policy.nt+2,1);
            for t=2:(policy.nt + 2)


                AL_1970(t) = AL_1970(t-1)*(1+AL_growth_1970(t));

            end

            % Steady State Productivity Growth

            AL_growth_iss_1970 = tpg_1970;
            AL_growth_fss_1970 = tpg_1970;

        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % 2015                                                     %
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            AL_growth_2015 = tpg_2015*ones(policy.nt+2,1);

            % Set total productivity

            AL_2015 = ones(policy.nt+2,1);
            for t=2:(policy.nt + 2)


                AL_2015(t) = AL_2015(t-1)*(1+AL_growth_2015(t));

            end

            % Steady State Productivity Growth

            AL_growth_iss_2015 = tpg_2015;
            AL_growth_fss_2015 = tpg_2015;
            
            
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Optimistic Transition                                    %
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            AL_growth_otransition = tpg_1970*ones(policy.nt+2,1);
            AL_growth_otransition(1:length(fernald_hp)) = fernald_hp;
        
            % Set total productivity

            AL_otransition = ones(policy.nt+2,1);
            for t=2:(policy.nt + 2)

                AL_otransition(t) = AL_otransition(t-1)*(1+AL_growth_otransition(t));

            end

            % Steady State Productivity Growth

            AL_growth_iss_otransition = tpg_1970;
            AL_growth_fss_otransition = tpg_1970;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Relative Price of Investment Goods  and Capital Specific Productivity
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        
        relp_1970 = 1.3;
        
        relp_vec_1970 = relp_1970*ones(policy.nt+2,1);
        
        relp_2015 = 1;
        
        relp_vec_2015 = relp_2015*ones(policy.nt+2,1);
        
        % We can't have a constant decline, because otherwise the interest
        % rate will jump at the point at which there is no longer a change
        % in the RELP
        
        relative_price_decline1 = [0:-.0006:-.012]';
        relative_price_decline2 = [-.012:.0006:0]';
        
        relative_price_decline = [relative_price_decline1;relative_price_decline2;];
        
        relative_price_vector = relp_1970*ones(length(relative_price_decline),1);
        
        for t=2:length(relative_price_decline)
            
            relative_price_vector(t) = relative_price_vector(t-1)*(1+relative_price_decline(t-1));
            
            
        end
        
        relative_price_transition = relative_price_vector(end)*ones(policy.nt+2,1);
        
        relative_price_transition(1:length(relative_price_vector)) = relative_price_vector;
        
        % Set Capital Productivity for transition
        
        capital_prod = 1./relative_price_vector;
        
        AK_transition = capital_prod(end)*ones(policy.nt+2,1);

        AK_transition(1:length(capital_prod)) = capital_prod;
         
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Mortality Transition                %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Read in the vector of marginal survival
        ms_transition_vec = csvread('data/jr_survival.csv');

        % Create the MS transition variable. It's final steady state values
        % will be the 2015 mortality data
        ms_transition = ms_transition_vec(end,:);
        ms_transition = repmat(ms_transition,policy.nt + 2,1);

        % Replace the rest of the transition with the data from 1970-2015
        ms_transition(1:size(ms_transition_vec,1),:) = ms_transition_vec;

        % Create survival vec
        s_transition = gens(ms_transition);

        % Generate mhs vec
        mhs_transition = genmhs(ms_transition);

        ms_1970 = ms_transition(1,:);
        ms_1970 = repmat(ms_1970,policy.nt + 2,1);
        s_1970 = gens(ms_1970);
        mhs_1970 = genmhs(ms_1970);

        ms_2015 = ms_transition(end,:);
        ms_2015 = repmat(ms_2015,policy.nt + 2,1);
        s_2015 = gens(ms_2015);
        mhs_2015 = genmhs(ms_2015);

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set up Government Debt                    %%
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        % MA Government Debt
        
            debt_data =  readtable('data/jr_gdebt.xls');
            
            gov_debt_ma = table2array(debt_data(:,'gov_debt_ma_paper1'));
            
        % 1970 Debt
        
            gov_debt_gdp_1970 = gov_debt_ma(1)*ones(policy.nt+2,1);
            
        % 2015 Spending, Deficit, and Debt
        
            gov_debt_gdp_2015 = gov_debt_ma(end)*ones(policy.nt+2,1);
            
        % Transition
            
            gov_debt_transition = gov_debt_gdp_2015;
            
            gov_debt_transition(1:length(gov_debt_ma)) = gov_debt_ma;
            
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Gov Spending                              %%
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        % 1970 Spending, Deficit, and Debt
        
            gov_amt_gdp_1970 = .2128*ones(policy.nt + 2,1);
            
        % 2015 Spending, Deficit, and Debt
        
            gov_amt_gdp_2015 = .2128*ones(policy.nt + 2,1);
            
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set up Net Foreign Assets                 %%
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        % For this exercise, set NFA to zero :)
        nfa_1970 = 0 ;
        
        % We will simply
        nfa_2015 = 0 ;
            
        nfa_transition_vector = 0;
        
        % Initial vectors
        
        nfa_gdp_1970 = nfa_1970*ones(policy.nt+2,1);
        
        nfa_gdp_2015 = nfa_2015*ones(policy.nt+2,1);
        
        % Transition
        
        nfa_gdp_transition = nfa_2015*ones(policy.nt+2,1);
        
        nfa_gdp_transition(1:length(nfa_transition_vector)) = nfa_transition_vector;
        
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set Up Monopolistic Competition           %%
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        theta_vec_1970 = theta_1970*ones(policy.nt+2,1);
        
        theta_vec_2015 = theta_2015*ones(policy.nt+2,1);
        
        % Transition
    
        scale = (theta_2015 - theta_1970)/40;
        theta_transition_vector = [theta_1970:scale:theta_2015]';
        
        theta_transition = theta_vec_2015;
        
        theta_transition(1:length(theta_transition_vector)) = theta_transition_vector;
        
        
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set up Debt Limit                         %%
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        dlim_pct_vec_1970 = dlim_pct_1970*ones(policy.nt+2,1);
        
        dlim_pct_vec_2015 = dlim_pct_2015*ones(policy.nt+2,1);
        
        % Transition vector
        scale = (dlim_pct_2015 - dlim_pct_1970)/40;
        
        dlim_transition_vector = [dlim_pct_1970(1):scale:dlim_pct_2015(1)]';
        
        dlim_transition = dlim_pct_vec_2015;
        
        dlim_transition(1:length(dlim_transition_vector)) = dlim_transition_vector;
        
          
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set up the Run Schedule                   %%
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        run_schedule.iss.do = 1; % Do the ISS
        run_schedule.fss.do = 0; % Calculate FSS
        run_schedule.trans.do = 1; % Calculate transition
        run_schedule.pe.run=0;
        
        run_schedule.tol = 1e-5;%1e-3; % Tolerance
        run_schedule.maxiter = 200; % Max # of iterations
        run_schedule.iss.algo = 'damp'; % Type of algorithm for ISS
        run_schedule.fss.algo = 'damp'; % Type of algorithm for FSS
        run_schedule.trans.algo = 'damp'; % Type of algorithm for transition
        
        run_schedule.iss.iter = 1;
        run_schedule.swdamp = .5; % Dampening on shadow wages
        run_schedule.kldamp = .4; % Dampening on capital
        run_schedule.tdamp = .5; % Dampening on Taxes
        
        run_schedule.matlab_solve = 0;
        run_schedule.jr_method = 1; % Whether we use my method to calculate optimal consumption

        run_schedule.fertility_ak = 1; % Fertility is set endogenously as in aka
        run_schedule.fertility_year = 60; % year after which fertility is endogenous
   
        run_schedule.demog_shift = 1; % this means we specify fertility beginning in year 1, rather than 25 years ago
        
        run_schedule.nfa_gdp = 1;% NFA will be sepcified as a % of GDP
        
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Initially Set Everything to 1970 and do a test run                %%
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
            % Set num kids to 1970 version
            
                p1.demog.num_kids = num_bb_1945;

            % Set the pop growth
            
                economy.n = bb_vec_1945;
            
            % Set labor productivity to 1970 vectors

                economy.ag.AL_growth = AL_growth_1970;
                economy.ag.AL = AL_1970;

                economy.ag.AL_growth_iss = AL_growth_iss_1970;
                economy.ag.AL_growth_fss = AL_growth_fss_1970;

            % Set mortality / survival
            
                p1.demog.s = s_1970; % survival
                p1.demog.ms = ms_1970; % marginal survival; column i = prob of surviving from i to j
                p1.demog.mhs = mhs_1970;

            % Set government spending
            
                gov.spend.amt_gdp = gov_amt_gdp_1970;

            % Specify debt the way JR Intended

                gov.debt_gdp = gov_debt_gdp_1970;
                
            % Specify Net Foreign Assets
            
                economy.nfa_gdp = nfa_gdp_1970;
                
            % Specify Monopolistic Competition
            
                prod.theta = theta_vec_1970;
                
                prod.monop_comp= 1; 
                
            % Specify Relative Price
            
                economy.relp = relp_vec_1970;
                
            % Specify debt limit
            
                 p1.util.dlim_pct = dlim_pct_vec_1970;
                 
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Run 1970 Equilibrium                      %%
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % If this is set equal to 1, run the transition
    % Else, we only want to run it for 1970
    
    if controller.transition==1
        
%             [ps_test,gov_test,economy_test,prices_test,nipa_bank_test] = master_control(p1,prod,gov,policy,economy,run_schedule);   
%             
%             % Original r
%             
%             r_1970 = prices_test.r(1);
%                    
%             nipa_display(nipa_bank_test,ps_test,prod,gov_test,policy,economy_test,prices_test,run_schedule,1,1) ;  
%             
%             nipa_display(nipa_bank_test,ps_test,prod,gov_test,policy,economy_test,prices_test,run_schedule,152,1) ;  
%             
%             % Save Comparison Y_pp and Y_pw:) 
%             comparison.Y_pp = nipa_bank_test.Y_pp;
%             comparison.Y_pw = nipa_bank_test.Y_pw;
%             comparison.wages = prices_test.wages;
%             
%             graph_display(nipa_bank_test,ps_test,prod,gov_test,policy,economy_test,prices_test,run_schedule,[],strcat('test_',controller.savename),comparison)
%                         
%     % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % Do Baby Boom Transition                   %%
%     % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%         % Set num kids to transition
%         p1.demog.num_kids = num_bb_transition;
%         
%         % Set the pop growth
%         economy.n = bb_transition;
%                         
%         % RUNNIT
%         
%         [ps_bb,gov_bb,economy_bb,prices_bb,nipa_bank_bb] = master_control(p1,prod,gov,policy,economy,run_schedule);   
%         
%         graph_display(nipa_bank_bb,ps_bb,prod,gov_bb,policy,economy_bb,prices_bb,run_schedule,[],strcat('baby_boom_',controller.savename),comparison)
%                 
%         % Reset num kids to 1970 version
% 
%         p1.demog.num_kids = num_bb_1945;
% 
%         % Reset the pop growth
% 
%         economy.n = bb_vec_1945;
%         
%     % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % Do Optimistic Productivity Transition     %%
%     % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             
%         % Set productivity to transition vector
%         
%         economy.ag.AL_growth = AL_growth_otransition;
%         economy.ag.AL = AL_otransition;
%         
%         economy.ag.AL_growth_iss = AL_growth_iss_otransition;
%         economy.ag.AL_growth_fss = AL_growth_fss_otransition;
%          
%         % RUNNIT
%         
%         [ps_oprod,gov_oprod,economy_oprod,prices_oprod,nipa_bank_oprod] = master_control(p1,prod,gov,policy,economy,run_schedule);   
% 
%         graph_display(nipa_bank_oprod,ps_oprod,prod,gov_oprod,policy,economy_oprod,prices_oprod,run_schedule,[],strcat('oprod_',controller.savename),comparison)
%                       
%         % Reset labor productivity to 1970 vectors
% 
%         economy.ag.AL_growth = AL_growth_1970;
%         economy.ag.AL = AL_1970;
% 
%         economy.ag.AL_growth_iss = AL_growth_iss_1970;
%         economy.ag.AL_growth_fss = AL_growth_fss_1970;
%         
%         
% %     % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %     % Do Mortality Transition                   %%
% %     % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %    
% %         % Set mortality / survival
% %         p1.demog.s = s_transition; % survival
% %         p1.demog.ms = ms_transition; % marginal survival; column i = prob of surviving from i to j
% %         p1.demog.mhs = mhs_transition;
% %         
% %         % Runnit
% %         
% %         [ps_mort,gov_mort,economy_mort,prices_mort,nipa_bank_mort] = master_control(p1,prod,gov,policy,economy,run_schedule);   
% % 
% %         graph_display(nipa_bank_mort,ps_mort,prod,gov_mort,policy,economy_mort,prices_mort,run_schedule,[],'mort',comparison)
% % 
% %         % Reset mortality / survival
% % 
% %         p1.demog.s = s_1970; % survival
% %         p1.demog.ms = ms_1970; % marginal survival; column i = prob of surviving from i to j
% %         p1.demog.mhs = mhs_1970;
% 
% % plot(1:152,100*prices_mort.r-1)
% %         xlabel('Period')
% %         ylabel('Rates')
% %         title('Interest Rates')
% %         print(strcat('figures/','ineq','_transition','.png'),'-dpng')
% %         
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Do Full Transition                        %%
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
       % Set productivity to transition vector
        
        economy.ag.AL_growth = AL_growth_transition;
        economy.ag.AL = AL_transition;
        
        economy.ag.AL_growth_iss = AL_growth_iss_transition;
        economy.ag.AL_growth_fss = AL_growth_fss_transition;
        
        % Set num kids to transition
        p1.demog.num_kids = num_bb_transition;
        
        % Set the pop growth
        economy.n = bb_transition;
        
        % Set mortality / survival
        p1.demog.s = s_transition; % survival
        p1.demog.ms = ms_transition; % marginal survival; column i = prob of surviving from i to j
        p1.demog.mhs = mhs_transition;
        
        % Specify debt the way JR Intended

        gov.debt_gdp = gov_debt_transition;
        
        % Specify net financial assets

        economy.nfa_gdp = nfa_gdp_transition;
        
        % Set Theta
        
        prod.theta = theta_transition;
        
        % Set Relative Prices
        economy.relp = relative_price_transition;
        
        % Specify debt limit
            
        p1.util.dlim_pct = dlim_transition;
                
        
        % RUNNIT
        
        [ps_full,gov_full,economy_full,prices_full,nipa_bank_full] = master_control(p1,prod,gov,policy,economy,run_schedule);   

        %nipa_display(nipa_bank_full,ps_full,prod,gov_full,policy,economy_full,prices_full,run_schedule,1,1) ;  
           
        %nipa_display(nipa_bank_full,ps_full,prod,gov_full,policy,economy_full,prices_full,run_schedule,152,1) ;  
            
        %graph_display(nipa_bank_full,ps_full,prod,gov_full,policy,economy_full,prices_full,run_schedule,[],strcat('full_',controller.savename))

    else
        
        % Else don't do the transition, only run the 1970 numbers
        
        run_schedule.iss.do = 1; % Do the ISS
        run_schedule.fss.do = 0; % Calculate FSS
        run_schedule.trans.do = 0; % Calculate transition
        
        [ps_test,gov_test,economy_test,prices_test,nipa_bank_test] = master_control(p1,prod,gov,policy,economy,run_schedule);   

        % Save original r

        r_1970 = prices_test.r(1);    
        
        
        % Now, run the counterfactual with theta set to 2015
        
        prod.theta = theta_vec_2015;
        
        [ps_count,gov_count,economy_count,prices_count,nipa_bank_count] = master_control(p1,prod,gov,policy,economy,run_schedule);   

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Print 1970 Results                               %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        fprintf(resultsfile,['\n']);
        fprintf(resultsfile,['1970 Counterfactual Theta Statistics']);
        fprintf(resultsfile,['\n']);
        fprintf(resultsfile,['\n']);
        fprintf(resultsfile,['1.Natural Interest Rate........................' num2str(prices_count.r(1),'%04.4f')  ]);
        fprintf(resultsfile,['\n']);
        fprintf(resultsfile,['2.Nom Investment Output........................' num2str(nipa_bank_count.IY_nom(1),'%04.4f')  ]);
        fprintf(resultsfile,['\n']);
        fprintf(resultsfile,['3.Labor Share...................................' num2str(nipa_bank_count.inc_share(1,1),'%04.4f')  ]);
        fprintf(resultsfile,['\n']);
        fprintf(resultsfile,['4.Debt to income................................' num2str(nipa_bank_count.debt_inc(1),'%04.4f')  ]);
        fprintf(resultsfile,['\n']);
        fprintf(resultsfile,['5.Bequests to income............................' num2str(nipa_bank_count.beq_inc(1),'%04.4f')  ]);
        fprintf(resultsfile,['\n']);
        fprintf(resultsfile,['6.Average Return................................' num2str(nipa_bank_count.average_return(1),'%04.4f')  ]);
        fprintf(resultsfile,['\n']);
        fprintf(resultsfile,['7.MPK...........................................' num2str(prices_count.mpk(1),'%04.4f')  ]);
        fprintf(resultsfile,['\n']);
        
        prod.theta = theta_vec_1970;

    end
    
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % Print 1970 Results                               %
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         
%         fprintf(resultsfile,['\n']);
%         fprintf(resultsfile,['1970 Steady State Statistics']);
%         fprintf(resultsfile,['\n']);
%         fprintf(resultsfile,['\n']);
%         fprintf(resultsfile,['1.Natural Interest Rate........................' num2str(prices_test.r(1),'%04.4f')  ]);
%         fprintf(resultsfile,['\n']);
%         fprintf(resultsfile,['2.Nom Investment Output........................' num2str(nipa_bank_test.IY_nom(1),'%04.4f')  ]);
%         fprintf(resultsfile,['\n']);
%         fprintf(resultsfile,['3.Labor Share...................................' num2str(nipa_bank_test.inc_share(1,1),'%04.4f')  ]);
%         fprintf(resultsfile,['\n']);
%         fprintf(resultsfile,['4.Debt to income................................' num2str(nipa_bank_test.debt_inc(1),'%04.4f')  ]);
%         fprintf(resultsfile,['\n']);
%         fprintf(resultsfile,['5.Bequests to income............................' num2str(nipa_bank_test.beq_inc(1),'%04.4f')  ]);
%         fprintf(resultsfile,['\n']);
%         fprintf(resultsfile,['6.Average Return................................' num2str(nipa_bank_test.average_return(1),'%04.4f')  ]);
%         fprintf(resultsfile,['\n']);
%         fprintf(resultsfile,['7.MPK...........................................' num2str(prices_test.mpk(1),'%04.4f')  ]);
%         fprintf(resultsfile,['\n']);
%         
%         % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         % Reset run schedule                        %%
%         % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%         run_schedule.iss.do = 1; % Do the ISS
%         run_schedule.fss.do = 0; % Calculate FSS
%         run_schedule.trans.do = 0; % Calculate transition
%         
%         close all
%         
%         % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         % Reset to 2015 Values                      %%
%         % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%             % Set productivity to 2015 vectors
% 
%                 economy.ag.AL_growth = AL_growth_2015;
%                 economy.ag.AL = AL_2015;
% 
%                 economy.ag.AL_growth_iss = AL_growth_iss_2015;
%                 economy.ag.AL_growth_fss = AL_growth_fss_2015;
% 
%             % Set num kids to 2015
% 
%                 p1.demog.num_kids = num_bb_2015;
% 
%             % Set the pop growth to 2015
% 
%                 economy.n = bb_vec_2015;
% 
%             % Set mortality / survival
%             
%                 p1.demog.s = s_2015; % survival
%                 p1.demog.ms = ms_2015; % marginal survival; column i = prob of surviving from i to j
%                 p1.demog.mhs = mhs_2015;
%                 
%             % Specify debt the way JR Intended
% 
%                 gov.debt_gdp = gov_debt_gdp_2015;
%                 
%             % Set government spending
%             
%                 gov.spend.amt_gdp = gov_amt_gdp_2015;
%                 
%             % Specify Net Foreign Assets
%             
%                 economy.nfa_gdp = nfa_gdp_2015;
%                 
%             % Specify Monopolistic Competition
%             
%                 prod.theta = theta_vec_2015;
%                 
%             % Specify Relative Price
%                 
%                 % Let's change this up a little bit
%                 economy.relp = relp_vec_2015;
%                 
%             % Specify debt limit
%             
%                  p1.util.dlim_pct = dlim_pct_vec_2015;
                
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Financial Frictions                       %%
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if controller.ff==1 
            
            % RUNNIT
            % Financial Frictions
            
            economy.ff = .04*ones(policy.nt+2,1);
            
            run_schedule.matlab_solve = 1;
            run_schedule.jr_method = 0; % Whether we use my method to calculate optimal consumption

            [ps_ff,gov_ff,economy_ff,prices_ff,nipa_bank_ff] = master_control(p1,prod,gov,policy,economy,run_schedule);   

            % Reset
            
            economy.ff = 0*ones(policy.nt+2,1);
            
            run_schedule.matlab_solve = 0;
            run_schedule.jr_method = 1; % Whether we use my method to calculate optimal consumption

        end
               
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Show Elasticities                         %%
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
    if controller.elasticities==1
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Loop through Pop/EIS Params               %%
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            eis_vec = [.5;gamma;.99];
            n_vec = [0;.01;.02;.03;.04];

            elas_n_mat = zeros(length(n_vec),length(eis_vec));

            for i=1:length(eis_vec)

                % Replace the EIS

                p1.util.gamma = eis_vec(i)*ones(policy.nt+2,lifespan);

                for j=1:length(n_vec)

                     economy.n = n_vec(j)*ones(policy.nt+2,1);

                     p1.demog.num_kids = ((1+n_vec(j))*ones(policy.nt + 2,1)).^(p1.demog.age_parent-1);

                     [ps_temp,gov_test,economy_temp,prices_temp,nipa_bank_temp] = master_control(p1,prod,gov,policy,economy,run_schedule);   

                     elas_n_mat(j,i) = prices_temp.r(1);

                end

            end
            
            % Put things back
            
                % Set num kids to 1970 version

                    p1.demog.num_kids = num_bb_2015;

                % Set the pop growth

                    economy.n = bb_vec_2015;
                    
                % EIS
                
                    p1.util.gamma = gamma*ones(policy.nt+2,lifespan);

            plot(n_vec,elas_n_mat,'LineWidth',2)
            legend('EIS = .5',strcat('EIS = ',num2str(gamma)),'EIS=1','Location','northoutside')
            xlabel('Population Growth Rate (%)','FontSize',12)
            ylabel('Real Interest Rate (%)','FontSize',12)
            set(gca,'FontSize',12)
            print(strcat('figures/real_rate_by_n_',controller.savename),'-dpng') 
            savefig(strcat('figures/real_rate_by_n_',controller.savename))
            % Calculate level of pop necessary to get to get real interest
            % rates to 1
            
            n_required = interp1(elas_n_mat(:,2),n_vec,.01,'pchip');
            tfr_required = 2*(1+n_required).^(p1.demog.age_parent(1)-1);
            
            tfr_vec = 2*(1+n_vec).^(p1.demog.age_parent(1)-1);
            
            % Calculate derivative of interest rate with respect to TFR
            
            deriv_tfr = (elas_n_mat(end,2) - elas_n_mat(1,2)) / (tfr_vec(end) - tfr_vec(1));
            
            % Elasticity
            % Take this at 2015 steady state, change from tfr of 2.5 to 2
            elas_tfr = ( (elas_n_mat(2,2) - elas_n_mat(1,2))/abs(elas_n_mat(1,2)))/ ( (tfr_vec(2) - tfr_vec(1))/abs(tfr_vec(1)));
            
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Loop through Government Debt              %%
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            eis_vec = [.5;gamma;.99];
            gdebt_vec = [0;.4;.8;1.2;1.6;2];

            elas_gov_mat = zeros(length(gdebt_vec),length(eis_vec));

            for i=1:length(eis_vec)

                % Replace the EIS

                p1.util.gamma = eis_vec(i)*ones(policy.nt+2,lifespan);

                for j=1:length(gdebt_vec)

                     gov.debt_gdp = gdebt_vec(j)*ones(policy.nt+2,1);

                     [ps_temp,gov_temp,economy_temp,prices_temp,nipa_bank_temp] = master_control(p1,prod,gov,policy,economy,run_schedule);   

                     elas_gov_mat(j,i) = prices_temp.r(1);

                end

            end
            
            % Put things back
            
            gov.debt_gdp = gov_debt_gdp_2015;
            
            p1.util.gamma = gamma*ones(policy.nt+2,lifespan);
            
            % Plot

            plot(gdebt_vec,elas_gov_mat,'LineWidth',2)
            legend('EIS = .5',strcat('EIS = ',num2str(gamma)),'EIS=1','Location','northoutside')
            xlabel('Gov Debt / GDP (%)','FontSize',12)
            ylabel('Real Interest Rate (%)','FontSize',12)
            set(gca,'FontSize',12)
            print(strcat('figures/real_rate_by_gdebt_',controller.savename),'-dpng') 
            savefig(strcat('figures/real_rate_by_gdebt_',controller.savename))
            
            % Calculate level of gdebt necessary to get to get real interest
            % rates to 1
            
            gdebt_required = interp1(elas_gov_mat(:,2),gdebt_vec,.01,'pchip');
            
            % Derivative of gdebt
            deriv_gdebt = (elas_gov_mat(end,2) - elas_gov_mat(1,2)) / (gdebt_vec(end) - gdebt_vec(1));
            
            % Elasticity
            % Take this at 2015 steady state, approx gdebt of 100%
            elas_gdebt = ( (elas_gov_mat(4,2) - elas_gov_mat(3,2))/abs(elas_gov_mat(3,2)))/ ( (gdebt_vec(4) - gdebt_vec(3))/abs(gdebt_vec(3)));
            
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Loop through borrowing limit              %%
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            eis_vec = [.5;gamma;.99];
            blim_vec = [0;.2;.4;.6;.8;1;1.2];

            elas_blim_mat = zeros(length(blim_vec),length(eis_vec));
            
            
            elas_blimdebt_mat = zeros(length(blim_vec),length(eis_vec));

            for i=1:length(eis_vec)

                % Replace the EIS

                p1.util.gamma = eis_vec(i)*ones(policy.nt+2,lifespan);

                for j=1:length(blim_vec)

                     p1.util.dlim_pct = blim_vec(j)*ones(policy.nt + 2,1);

                     [ps_temp,gov_temp,economy_temp,prices_temp,nipa_bank_temp] = master_control(p1,prod,gov,policy,economy,run_schedule);   

                     elas_blim_mat(j,i) = prices_temp.r(1);
                     
                     elas_blimdebt_mat(j,i) = nipa_bank_temp.debt_inc(1);

                end

            end

            % Put things back
            
            p1.util.dlim_pct = dlim_pct_vec_2015;
            
            p1.util.gamma = gamma*ones(policy.nt+2,lifespan);
            
            % Plot

            plot(blim_vec,elas_blim_mat,'LineWidth',2)
            legend('EIS = .5',strcat('EIS = ',num2str(gamma)),'EIS=1','Location','northoutside')
            xlabel('Bor. Lim (% Income)','FontSize',12)
            ylabel('Real Interest Rate (%)','FontSize',12)
            set(gca,'FontSize',12)
            print(strcat('figures/real_rate_by_blim_',controller.savename),'-dpng') 
            
            % Plot debt as % of GDP

            plot(elas_blimdebt_mat,elas_blim_mat,'LineWidth',2)
            legend('EIS = .5',strcat('EIS = ',num2str(gamma)),'EIS=1','Location','northoutside')
            xlabel('Household debt (% GDP)','FontSize',12)
            ylabel('Real Interest Rate (%)','FontSize',12)
            set(gca,'FontSize',12)
            print(strcat('figures/real_rate_by_blimdebt_',controller.savename),'-dpng') 
            savefig(strcat('figures/real_rate_by_blimdebt_',controller.savename))
            
            % Calculate level of blim necessary to get to get real interest
            % rates to 1
            
            blim_required = interp1(elas_blim_mat(:,2),blim_vec,.01,'pchip');
            
            % Derivative of blim
            deriv_blim = (elas_blim_mat(end,2) - elas_blim_mat(1,2)) / (blim_vec(end) - blim_vec(1));
            
            % Elasticity
            % Take this at 2015 steady state, approx .4
            elas_blim = ( (elas_blim_mat(3,2) - elas_blim_mat(2,2))/abs(elas_blim_mat(2,2)))/ ( (blim_vec(3) - blim_vec(2))/abs(blim_vec(2)));
            
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Loop through Foreign Savings Glut         %%
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            eis_vec = [.5;gamma;.99];
            nfa_vec = [0;.2;.4;.6;.8;1;1.2];

            elas_nfa_mat = zeros(length(nfa_vec),length(eis_vec));

            for i=1:length(eis_vec)

                % Replace the EIS

                p1.util.gamma = eis_vec(i)*ones(policy.nt+2,lifespan);

                for j=1:length(nfa_vec)

                      economy.nfa_gdp = nfa_vec(j)*ones(policy.nt+2,1);

                     [ps_temp,gov_temp,economy_temp,prices_temp,nipa_bank_temp] = master_control(p1,prod,gov,policy,economy,run_schedule);   

                     elas_nfa_mat(j,i) = prices_temp.r(1);

                end

            end

            
            % Put things back
            
            economy.nfa_gdp = nfa_gdp_2015;
            
            p1.util.gamma = gamma*ones(policy.nt+2,lifespan);
            
            % Plot

            plot(nfa_vec,elas_nfa_mat,'LineWidth',2)
            legend('EIS = .5',strcat('EIS = ',num2str(gamma)),'EIS=1','Location','northoutside')
            xlabel('NFA / GDP (%)','FontSize',12)
            ylabel('Real Interest Rate (%)','FontSize',12)
            set(gca,'FontSize',12)
            print(strcat('figures/real_rate_by_nfa_',controller.savename),'-dpng') 
            savefig(strcat('figures/real_rate_by_nfa_',controller.savename))
            % Calculate level of nfa necessary to get to get real interest
            % rates to 1
            
            nfa_required = interp1(elas_nfa_mat(:,2),nfa_vec,.01,'pchip');
            
            % Derivative of nfa
            deriv_nfa = (elas_nfa_mat(end,2) - elas_nfa_mat(1,2)) / (nfa_vec(end) - nfa_vec(1));
            
            % Elasticity
            % Take this at 2015 steady state, approx 0
            elas_nfa = ( (elas_nfa_mat(3,2) - elas_nfa_mat(2,2))/abs(elas_nfa_mat(2,2)))/ ( (nfa_vec(3) - nfa_vec(2))/abs(nfa_vec(2)));
            
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Loop through Profit Share                 %%
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            eis_vec = [.5;gamma;.99];
            theta_vec = [5;6;7;8;9;10;11];

            elas_profit_mat = zeros(length(theta_vec),length(eis_vec));
            
            % Labor share 
            ls_profit_mat = zeros(length(theta_vec),length(eis_vec));

            for i=1:length(eis_vec)

                % Replace the EIS

                p1.util.gamma = eis_vec(i)*ones(policy.nt+2,lifespan);

                for j=1:length(theta_vec)

                     prod.theta = theta_vec(j)*ones(policy.nt+2,1);

                     [ps_temp,gov_temp,economy_temp,prices_temp,nipa_bank_temp] = master_control(p1,prod,gov,policy,economy,run_schedule);   

                     elas_profit_mat(j,i) = prices_temp.r(1);
                     
                     ls_profit_mat(j,i) = nipa_bank_temp.inc_share(1,1);

                end

            end

             
            % Put things back
            
            prod.theta = theta_vec_2015;
            
            p1.util.gamma = gamma*ones(policy.nt+2,lifespan);
            
            % Plot

            % Profit share vector
            prof_share_temp = 100./theta_vec;
            
            plot(prof_share_temp,elas_profit_mat,'LineWidth',2)
            legend('EIS = .5',strcat('EIS = ',num2str(gamma)),'EIS=1','Location','northoutside')
            xlabel('Profit / GDP (%)','FontSize',12)
            ylabel('Real Interest Rate (%)','FontSize',12)
            set(gca,'FontSize',12)
            print(strcat('figures/real_rate_by_pshare_',controller.savename),'-dpng') 
            
            % Labor share graph
            plot(ls_profit_mat,elas_profit_mat,'LineWidth',2)
            legend('EIS = .5',strcat('EIS = ',num2str(gamma)),'EIS=1','Location','northoutside')
            xlabel('Labor share','FontSize',12)
            ylabel('Real Interest Rate (%)','FontSize',12)
            set(gca,'FontSize',12)
            print(strcat('figures/real_rate_by_psharels_',controller.savename),'-dpng') 
            savefig(strcat('figures/real_rate_by_psharels_',controller.savename))
            
            % Calculate level of pshare necessary to get to get real interest
            % rates to 1
            
            profit_required = interp1(elas_profit_mat(:,2),prof_share_temp,.01,'pchip');
            
            % Derivative of profit
            deriv_profit = (elas_profit_mat(end,2) - elas_profit_mat(1,2)) / (prof_share_temp(end) - prof_share_temp(1));
            
            % Elasticity
            % Take this at 2015 steady state, approx 0
            elas_profit = ( (elas_profit_mat(2,2) - elas_profit_mat(1,2))/abs(elas_profit_mat(1,2)))/ ( (prof_share_temp(2) - prof_share_temp(1))/abs(prof_share_temp(1)));
            
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Loop through RELP                                    %%
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            prod_eis_vec = [1;sigma;1.25];
            relp_vec = [.6;.7;.8;1;1.2];

            elas_relp_mat = zeros(length(relp_vec),length(prod_eis_vec));

            for i=1:length(prod_eis_vec)

                % Replace the EIS

                prod.sigma = prod_eis_vec(i)*ones(policy.nt+2,1);

                for j=1:length(relp_vec)

                     % Set Relative Price
                     economy.relp = (relp_vec(j))*ones(policy.nt+2,1);
                     
                     [ps_temp,gov_temp,economy_temp,prices_temp,nipa_bank_temp] = master_control(p1,prod,gov,policy,economy,run_schedule);   

                     elas_relp_mat(j,i) = prices_temp.r(1);

                end

            end
            
            % Put things back
            
            economy.relp = relp_vec_2015;
            prod.sigma = sigma*ones(policy.nt+2,1);

            % Plot
            
            plot(relp_vec,elas_relp_mat,'LineWidth',2)
            legend('EIS = 1',strcat('EIS = ',num2str(sigma)),'EIS=1.25','Location','northoutside')
            xlabel('Relp. K. ','FontSize',12)
            ylabel('Real Interest Rate (%)','FontSize',12)
            set(gca,'FontSize',12)
            print(strcat('figures/real_rate_by_relp_',controller.savename),'-dpng') 
            savefig(strcat('figures/real_rate_by_relp_',controller.savename))
            
            % Calculate level of relp necessary to get to get real interest
            % rates to 1
            
            relp_required = interp1(elas_relp_mat(:,2),relp_vec,.01,'pchip');
            
            % Derivative of relp
            deriv_relp = (elas_relp_mat(end,1) - elas_relp_mat(1,1)) / (relp_vec(end) - relp_vec(1));
            
            % Elasticity
            % Take this at 2015 steady state, 1
            elas_relp = ( (elas_relp_mat(4,1) - elas_relp_mat(3,1))/abs(elas_relp_mat(3,1)))/ ( (relp_vec(4) - relp_vec(3))/abs(relp_vec(3)));
            
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Loop through Productivity/EIS             %%
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            eis_vec = [.5;gamma;.99];
            prod_vec = [.0;.005;.01;.015;.02;.03];

            elas_g_mat = zeros(length(prod_vec),length(eis_vec));

            for i=1:length(eis_vec)

                % Replace the EIS

                p1.util.gamma = eis_vec(i)*ones(policy.nt+2,lifespan);

                for j=1:length(prod_vec)

                    tpg = prod_vec(j);

                    % Hicks Neutral Productivity
                    economy.ag.AL_growth = tpg*ones(policy.nt+2,1);

                    % Set the labor productivity, capital prod, total prod

                    economy.ag.AL = ones(policy.nt+2,1);

                    for t=2:(policy.nt + 2)

                        economy.ag.AL(t) = economy.ag.AL(t-1)*(1+economy.ag.AL_growth(t));

                    end

                    % Steady State Productivity Growth

                    economy.ag.AL_growth_iss = tpg;
                    economy.ag.AL_growth_fss = tpg;

                    [ps_temp,gov_temp,economy_temp,prices_temp,nipa_bank_temp] = master_control(p1,prod,gov,policy,economy,run_schedule);   

                    elas_g_mat(j,i) = prices_temp.r(1);

                end

            end

            % Put Things Back
            
            p1.util.gamma = gamma*ones(policy.nt+2,lifespan);
            
            economy.ag.AL_growth = AL_growth_2015;
            economy.ag.AL = AL_2015;

            economy.ag.AL_growth_iss = AL_growth_iss_2015;
            economy.ag.AL_growth_fss = AL_growth_fss_2015;
            
            % plot

            plot(prod_vec,elas_g_mat,'LineWidth',2)
            legend('EIS = .5',strcat('EIS = ',num2str(gamma)),'EIS=1','Location','northoutside')
            xlabel('Productivity Growth Rate (%)','FontSize',12)
            ylabel('Real Interest Rate (%)','FontSize',12)
            set(gca,'FontSize',12)
            print(strcat('figures/real_rate_by_g_',controller.savename),'-dpng') 
            savefig(strcat('figures/real_rate_by_g_',controller.savename))
            
            % Calculate level of prod necessary to get to get real interest
            % rates to 1
            
            g_required = interp1(elas_g_mat(:,2),prod_vec,.01,'pchip');
            
            % Derivative of nfa
            deriv_prod = (elas_g_mat(end,2) - elas_g_mat(1,2)) / (prod_vec(end) - prod_vec(1));
            
            % Elasticity
            % Take this at 2015 steady state, approx .005
            elas_prod = ( (elas_g_mat(3,2) - elas_g_mat(2,2))/abs(elas_g_mat(2,2)))/ ( (prod_vec(3) - prod_vec(2))/abs(prod_vec(2)));
            
            fprintf(resultsfile,['\n']);
            fprintf(resultsfile,['Table: Necessary Changes to Achieve Real Rate of 1%']);
            fprintf(resultsfile,['\n']);
            fprintf(resultsfile,['\n']);
            fprintf(resultsfile,['1.Total Fertility...............................' num2str(tfr_required,'%04.4f')  ]);
            fprintf(resultsfile,['\n']);
            fprintf(resultsfile,['2.Government Debt...............................' num2str(gdebt_required,'%04.4f')  ]);
            fprintf(resultsfile,['\n']);
            fprintf(resultsfile,['3.Borrowing Limit...............................' num2str(blim_required,'%04.4f')  ]);
            fprintf(resultsfile,['\n']);
            fprintf(resultsfile,['4.Foreign Assets................................' num2str(nfa_required,'%04.4f')  ]);
            fprintf(resultsfile,['\n']);
            fprintf(resultsfile,['5.Relative Price................................' num2str(relp_required,'%04.4f')  ]);
            fprintf(resultsfile,['\n']);
            fprintf(resultsfile,['5.Productivity..................................' num2str(g_required,'%04.4f')  ]);
            fprintf(resultsfile,['\n']);
            
    end
    
        
%     % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % Original 2015 Interest Rate               %%
%     % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%             [ps_decomp,gov_decomp,economy_decomp,prices_decomp,nipa_bank_decomp] = master_control(p1,prod,gov,policy,economy,run_schedule);   
%             
%             r_2015 = prices_decomp.r(1);
%     
%         fprintf(resultsfile,['\n']);
%         fprintf(resultsfile,['2015 Steady State Statistics']);
%         fprintf(resultsfile,['\n']);
%         fprintf(resultsfile,['\n']);
%         fprintf(resultsfile,['1.Natural Interest Rate........................' num2str(prices_decomp.r(1),'%04.4f')  ]);
%         fprintf(resultsfile,['\n']);
%         fprintf(resultsfile,['2.Nom Investment Output.........................' num2str(nipa_bank_decomp.IY_nom(1),'%04.4f')  ]);
%         fprintf(resultsfile,['\n']);
%         fprintf(resultsfile,['3.Labor Share...................................' num2str(nipa_bank_decomp.inc_share(1,1),'%04.4f')  ]);
%         fprintf(resultsfile,['\n']);
%         fprintf(resultsfile,['4.Debt to income................................' num2str(nipa_bank_decomp.debt_inc(1),'%04.4f')  ]);
%         fprintf(resultsfile,['\n']);
%         fprintf(resultsfile,['5.Bequests to income............................' num2str(nipa_bank_decomp.beq_inc(1),'%04.4f')  ]);
%         fprintf(resultsfile,['\n']);
%         fprintf(resultsfile,['6.Average Return................................' num2str(nipa_bank_decomp.average_return(1),'%04.4f')  ]);
%         fprintf(resultsfile,['\n']);
%         fprintf(resultsfile,['7.MPK...........................................' num2str(prices_decomp.mpk(1),'%04.4f')  ]);
%         fprintf(resultsfile,['\n']);
        
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate Decomposition                   %%
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if controller.decomposition==1
            
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1970 Fertility                            %%
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        % Set num kids to 1970 version

            p1.demog.num_kids = num_bb_1945;

        % Set the pop growth

            economy.n = bb_vec_1945;
            
         % Get the interest rate
         
            [ps_decomp,gov_decomp,economy_decomp,prices_decomp,nipa_bank_decomp] = master_control(p1,prod,gov,policy,economy,run_schedule);   
            
            r_1970_n = prices_decomp.r(1);
            
        % Reset Values
        % Set num kids to 2015
        % Set the pop growth to 2015

            p1.demog.num_kids = num_bb_2015;
            economy.n = bb_vec_2015;
            
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1970 Productivity                         %%
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        % Set labor productivity to 1970 vectors

            economy.ag.AL_growth = AL_growth_1970;
            economy.ag.AL = AL_1970;

            economy.ag.AL_growth_iss = AL_growth_iss_1970;
            economy.ag.AL_growth_fss = AL_growth_fss_1970;
            
         % Get the interest rate
         
            [ps_decomp,gov_decomp,economy_decomp,prices_decomp,nipa_bank_decomp] = master_control(p1,prod,gov,policy,economy,run_schedule);   
            
            r_1970_g = prices_decomp.r(1);
                
        % Reset productivity to 2015 vectors

            economy.ag.AL_growth = AL_growth_2015;
            economy.ag.AL = AL_2015;

            economy.ag.AL_growth_iss = AL_growth_iss_2015;
            economy.ag.AL_growth_fss = AL_growth_fss_2015;
            
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1970 Mortality                            %%
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        % Set 1970 mortality / survival

            p1.demog.s = s_1970; % survival
            p1.demog.ms = ms_1970; % marginal survival; column i = prob of surviving from i to j
            p1.demog.mhs = mhs_1970;
    
         % Get the interest rate
         
            [ps_decomp,gov_decomp,economy_decomp,prices_decomp,nipa_bank_decomp] = master_control(p1,prod,gov,policy,economy,run_schedule);   
            
            r_1970_mort = prices_decomp.r(1);
            
        % Reset 2015 mortality / survival

            p1.demog.s = s_2015; % survival
            p1.demog.ms = ms_2015; % marginal survival; column i = prob of surviving from i to j
            p1.demog.mhs = mhs_2015;
            
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1970 Government debt                      %%
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        % 1970 debt

            gov.debt_gdp = gov_debt_gdp_1970;
            
         % Get the interest rate
         
            [ps_decomp,gov_decomp,economy_decomp,prices_decomp,nipa_bank_decomp] = master_control(p1,prod,gov,policy,economy,run_schedule);   
            
            r_1970_gdebt = prices_decomp.r(1);
                
        % Reset 2015 debt

            gov.debt_gdp = gov_debt_gdp_2015;
            
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1970 Profits                              %%
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        % Specify 1970 theta

            prod.theta = theta_vec_1970;
            
         % Get the interest rate
         
            [ps_decomp,gov_decomp,economy_decomp,prices_decomp,nipa_bank_decomp] = master_control(p1,prod,gov,policy,economy,run_schedule);   
            
            r_1970_pshare = prices_decomp.r(1);
            
        % Reset 2015 theta

            prod.theta = theta_vec_2015;
            
            
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1970 RELP                                 %%
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        % Specify Relative Price

            economy.relp = relp_vec_1970;
            
         % Get the interest rate
         
            [ps_decomp,gov_decomp,economy_decomp,prices_decomp,nipa_bank_decomp] = master_control(p1,prod,gov,policy,economy,run_schedule);   
            
            r_1970_relp = prices_decomp.r(1);
 
        % Replace relative price

            economy.relp = relp_vec_2015;
            
            
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 1970 Debt limit                           %%
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        % 1970 dlim

            p1.util.dlim_pct = dlim_pct_vec_1970;
            
         % Get the interest rate
         
            [ps_decomp,gov_decomp,economy_decomp,prices_decomp,nipa_bank_decomp] = master_control(p1,prod,gov,policy,economy,run_schedule);   
            
            r_1970_dlim = prices_decomp.r(1);
                
        % Reset 2015 dlim

            p1.util.dlim_pct = dlim_pct_vec_2015;

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Display Decomposition                     %%
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        r_change_prod = r_2015 - r_1970_g;

        r_change_bb = r_2015 - r_1970_n;
        
        r_change_mort = r_2015 - r_1970_mort;
        
        r_change_gdebt = r_2015 - r_1970_gdebt;
        
        r_change_pshare = r_2015 - r_1970_pshare;
        
        r_change_relp = r_2015 - r_1970_relp;
        
        r_change_dlim = r_2015 - r_1970_dlim;
        
        r_change_sum = r_change_prod + r_change_bb + r_change_mort + r_change_gdebt + r_change_pshare + r_change_relp + r_change_dlim;
        
        r_change_prod_pct = 100*r_change_prod/r_change_sum;
        r_change_bb_pct = 100*r_change_bb/r_change_sum;
        r_change_mort_pct = 100*r_change_mort/r_change_sum;
        r_change_gdebt_pct = 100*r_change_gdebt/r_change_sum;
        r_change_pshare_pct = 100*r_change_pshare/r_change_sum;
        r_change_relp_pct = 100*r_change_relp/r_change_sum;
        r_change_dlim_pct = 100*r_change_dlim/r_change_sum;
        
        r_change_total = r_2015 - r_1970;
        
        fprintf(resultsfile,['\n']);
        fprintf(resultsfile,['Decomposition of Interest Rate Change'])
        fprintf(resultsfile,['\n']);
        fprintf(resultsfile,['\n']);
        fprintf(resultsfile,['1.Total interest change.........................' num2str(r_change_total,'%04.4f') '.... 100'  ]);
        fprintf(resultsfile,['\n']);
        fprintf(resultsfile,['2.Productivity Interest Rate Change.............' num2str(r_change_prod,'%04.4f')  '.... ' num2str(r_change_prod_pct,'%03.0f')  ]);
        fprintf(resultsfile,['\n']);
        fprintf(resultsfile,['3.Population Interest Rate Change...............' num2str(r_change_bb,'%04.4f')  '.... ' num2str(r_change_bb_pct,'%03.0f')  ]);
        fprintf(resultsfile,['\n']);
        fprintf(resultsfile,['4.Mortality Change..............................' num2str(r_change_mort,'%04.4f')  '.... ' num2str(r_change_mort_pct,'%03.0f')  ]);
        fprintf(resultsfile,['\n']);
        fprintf(resultsfile,['5.Gov Debt Change................................' num2str(r_change_gdebt,'%04.4f')  '.... ' num2str(r_change_gdebt_pct,'%03.0f')  ]);
        fprintf(resultsfile,['\n']);
        fprintf(resultsfile,['6.Labor Share Change............................' num2str(r_change_pshare,'%04.4f')  '.... ' num2str(r_change_pshare_pct,'%03.0f')  ]);
        fprintf(resultsfile,['\n']);
        fprintf(resultsfile,['7.Relp Share Change.............................' num2str(r_change_relp,'%04.4f')  '.... ' num2str(r_change_relp_pct,'%03.0f')  ]);
        fprintf(resultsfile,['\n']);
        fprintf(resultsfile,['7.Dlim Share Change.............................' num2str(r_change_dlim,'%04.4f')  '.... ' num2str(r_change_dlim_pct,'%03.0f')  ]);
        fprintf(resultsfile,['\n']);
        fprintf(resultsfile,['\n']);
        
    end
    

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Secular Stagnation Exercise               %%
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
    if controller.sec_stag==1
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Step 1: Reset Beta and Save Statistics               %%
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
            % First, we need to choose a beta that corresponds to an output
            % gap (needs to be below -2% with a 2% inflation target)
            
            % After choosing beta, we run the economy to find the natural
            % rate of interest. Then, we save the debt limit / level of
            % government spending from this. 
    
            % Set Beta
            p1.util.delta = params.delta_ss;
            
            % Run the economy to save values from
            [ps_ss,gov_ss,economy_ss,prices_ss,nipa_bank_ss] = master_control(p1,prod,gov,policy,economy,run_schedule);   
            
            % Display some statistics from the flex price economy
            nipa_display(nipa_bank_ss,ps_ss,prod,gov_ss,policy,economy_ss,prices_ss,run_schedule,1,1) ;  
            
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Step 2: Change economy structure          %%
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
            % For the sec stag economy, government spending and borrowing
            % limits are simply a number. This code changes the economy's
            % structure to ensure this. 
    
            % Debt limit to be not a percent of income
            p1.dlim = 2;
            p1.util.dlim_pct = 0*ones(policy.nt + 2,1); % debt limit as a percent of income. Only valid when p1.dlim = 1
            p1.util.dlim = ps_ss.opt.dlim; %Set debt limit from where it was before
           
            % Get Productivity Adjustment -- in the new economy, we don't
            % adjust productivity so the wage is 1
            economy.ag.A = economy_ss.ag.A_adj;
            economy.adj_prod = 0;

            % Add Government Debt 
            
            gov.debt = 0*ones(policy.nt + 2,1); % Debt as a percentage of tte capital stock
            gov.debt_gdp = 0*ones(policy.nt+2,1);
            gov.debt_alt = gov_ss.debt.*economy_ss.ag.K.*ones(policy.nt + 2,1); % GRoss debt amount

        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Recalculate Equilibrium to Test           %%
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
            % This recalculates the equilibrium to make sure nothing has
            % changed
            
            [ps_ss,gov_ss,economy_ss,prices_ss,nipa_bank_ss] = master_control(p1,prod,gov,policy,economy,run_schedule);   

            % Save original r
            orig_r = prices_ss.r(1);
            
            % Save original y
            orig_y = economy_ss.ag.Y(1);
                  
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Step 3: Create Aggregate Demand Curve     %%
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
            % This generates an employment / inflation aggregate demand
            % curve. 

            % First, we set Monetary Policy
            
            % Pi Star is inflation target. Usually this will be 2%
            pi_star = .02;
            
            % i_star is the target nominal interest rate
            % This is set as (1+r_natural)*(1+pi_star)
            i_star = (1+orig_r)*(1+pi_star) - 1;
                        
            % Set monetary policy phi parameter 
            phi_pi = 2;

            % Kink is the portion at which the AD curve shifts. 
            % This can be calculated from equation 23: it is the level of
            % inflation such that you are constrained by the lower bound
            
            % If inflation is lower, then the nominal rate is 0, and the
            % real rate will simply be the 1 devided by the inflation rate
            kink = (1+i_star)^(-1/phi_pi)*(1+pi_star) - 1;

            % Remember -- our program takes in labor supply, and spits out
            % the equil real interest rate. This is exactly what we need
            % to do the IS curve. 
            
            % We will then convert the real interest rate to an inflation
            % rate. 
            
            LS_vec = [.80:.04:1]'; % labor supply vector
            
            Y_vec = zeros(length(LS_vec),1); % Output vector
            
            numl = length(LS_vec);
            r_vec = zeros(numl,1);
            adtop_vec = zeros(numl,1);

            for index=1:numl

                % Set the employment
                economy.employment = LS_vec(index);

                % Calculate equilibrium interest rate
                [ps_temp,gov_temp,economy_temp,prices_temp,nipa_bank_temp] = master_control(p1,prod,gov,policy,economy,run_schedule);   

                % Save interest rate
                r_vec(index) = prices_temp.r(1);
                
                % Save output as a fraction of full employment output
                Y_vec(index) = economy_temp.ag.Y(1) / orig_y;

            end

            % Top of AD
            % Remember -- for the top of the AD curve, higher inflation means
            % higher nomianl interest rates, thus higher real rates. 

            % All we have is real rates -- we need to convert these to
            % inflation rates. To do this, we use equation (23) in my version
            % of the sec stag paper. This calculates (1 + i). To calculate the
            % real interest rate, we must devide by inflation. 
            for index=1:numl

                % calculate inflation rates
                adtop_vec(index) = ((1+r_vec(index))*(1+pi_star)^(phi_pi)/(1+i_star))^(1/(phi_pi-1)) - 1;

            end

            % The top of the aggregate demand curve is the portion above the
            % kink
            real_atop = adtop_vec(find(adtop_vec >= kink));
            ls_atop = LS_vec(find(adtop_vec >= kink));
            y_atop = Y_vec(find(adtop_vec >= kink));

            % Bottom of the aggregate demand curve; now, we are constrained by
            % the zero lower bound. So what ever the real rate is, the
            % inflation is such that (1+i)=1=(1+r)*(1+pi)
            pi_vec = 1./(1+r_vec) - 1;

            % The aggregate demand is below the kink
            real_abottom = pi_vec(find(pi_vec < kink ));
            ls_abottom = LS_vec(find(pi_vec < kink ));
            y_abottom = Y_vec(find(pi_vec < kink ));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Step 4: Draw AS Curve                            %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Aggregate supply curve follows directly  from an equation in
        % Eggertsson & Mehrotra
        
        alpha = epsilon; % Capital Share
        pi_alt_rates = -[-pi_star:.00005:.005]'; % Inflation rates to loop through
        employment = zeros(length(pi_alt_rates),1); % Output as a percentage of full employment
        for index = 1:length(pi_alt_rates)

            % This is equation B23 from appendix of Eggertsson & mehrotra.
            employment(index) = ( (1-gamma_as*(1+pi_star)/(1+pi_alt_rates(index)))/(1-gamma_as) )^(1/alpha);

        end

        % Trim values to reasonable values for graphing purposes
        trim_values = find(employment > .1);
        employment = employment(trim_values);
        pi_alt_rates = pi_alt_rates(trim_values);

        % Now, get the vertical part of the aggregate supply

        as_ver_ls = ones(5,1);
        as_ver_pi = pi_star:.005:(pi_star + .005*4);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Step 5: Calculate intersection of AS and AD  Curve  %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        % Now get the intersection
        
            % Interpolate Aggregate Demand
        
            ls_d_interp = interp1(real_abottom,ls_abottom,pi_alt_rates,'pchip');
        
            % Find the intersection
            
            [value,intersection] = min(abs(employment - ls_d_interp));
            
            if abs(value) > .01
                
                error('WARNING NON CONVERGENCE OF SEC STAG EQUILIBRIUM')
                
            end
        
            % Employment
            
            employment_eq = employment(intersection);
            
            inflation_eq = pi_alt_rates(intersection);
            
            output_eq = interp1(LS_vec,Y_vec,employment_eq);
            
            output_supply = interp1(LS_vec,Y_vec,employment);
            
        % Now get the information on the secular stagnation equilibrium
        
            % Interpolate the relationship between employment and output
            
                economy.employment = employment_eq;

                %p1.demog.hc = economy.employment*repmat(hc',policy.nt + 2,1) .* ones(policy.nt+2,lifespan);

                [ps_ssss,gov_ssss,economy_ssss,prices_ssss,nipa_bank_ssss] = master_control(p1,prod,gov,policy,economy,run_schedule);   
                
                % Replace employment
                
                economy.employment = 1;
        % Plot the Employment Results 
        
        close all
        
        plot(ls_abottom,real_abottom,ls_atop,real_atop,employment,pi_alt_rates,as_ver_ls,as_ver_pi,'LineWidth',2)
        legend('Aggregate Demand 1','Aggregate Demand 2','Aggregate Supply 1','Aggregate Supply 2','Location','northoutside')
        xlabel('Employment Rate')
        ylabel('Inflation Rate (%)')
        axis([.2 1.2 0 inf])
        print(strcat('figures/as_ad_emp_',controller.savename),'-dpng') 
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Step 5: Plot Output Figure                       %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
        % Extend these things
        y_kink = interp1(real_abottom,y_abottom,kink,'pchip');
        
%         plot([y_abottom;y_kink],[real_abottom;kink],[y_atop;y_kink],[real_atop;kink],output_supply,pi_alt_rates,as_ver_ls,as_ver_pi,'LineWidth',2)
%         legend('Aggregate Demand 1','Aggregate Demand 2','Aggregate Supply 1','Aggregate Supply 2','Location','northoutside')
%         xlabel('Output % of FE')
%         ylabel('Inflation Rate (%)')
%         axis([.2 1.2 0 inf])
%         print(strcat('figures/as_ad_out_',controller.savename),'-dpng') 
%         
        plot([y_abottom;y_kink],[real_abottom;kink],'r',output_supply,pi_alt_rates,'b',[y_atop;y_kink],[real_atop;kink],'r',as_ver_ls,as_ver_pi,'b','LineWidth',2)
        legend('Aggregate Demand','Aggregate Supply','Location','northoutside')
        xlabel('Output % of FE')
        ylabel('Inflation Rate (%)')
        axis([.2 1.2 0 inf])
        print(strcat('figures/as_ad_out_',controller.savename),'-dpng') 
        
        % Export Table of These Things
        sec_table_output = table([y_abottom;y_kink;zeros(length(pi_alt_rates)-(length(y_abottom)+1),1)],'VariableNames',{'ad_output_bottom'});
        
        sec_table_output = [sec_table_output,table([real_abottom;kink;zeros(length(pi_alt_rates)-(length(y_abottom)+1),1)],'VariableNames',{'ad_inflation_bottom'})];

        sec_table_output = [sec_table_output,table(output_supply,'VariableNames',{'as_output_bottom'})];
        
        sec_table_output = [sec_table_output,table(pi_alt_rates,'VariableNames',{'as_inflation_bottom'})];
        
        sec_table_output = [sec_table_output,table([y_atop;y_kink;zeros(length(pi_alt_rates)-(length(y_atop)+1),1)],'VariableNames',{'ad_output_top'})];
        
        sec_table_output = [sec_table_output,table([real_atop;kink;zeros(length(pi_alt_rates)-(length(y_atop)+1),1)],'VariableNames',{'ad_inflation_top'})];

        writetable(sec_table_output,strcat('figures/',controller.savename,'_sec_tables','.csv'))
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Step 6: Table of Sec Stag Results                %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        fprintf(resultsfile,['Secular Stagnation Statistics']);
        fprintf(resultsfile,['\n']);
        fprintf(resultsfile,['\n']);
        fprintf(resultsfile,['1.Natural Interest Rate........................' num2str(prices_ss.r(1),'%04.4f')  ]);
        fprintf(resultsfile,['\n']);
        fprintf(resultsfile,['2.Real Interest Rate...........................' num2str(prices_ssss.r(1),'%04.4f')  ]);
        fprintf(resultsfile,['\n']);
        fprintf(resultsfile,['3.Inflation.....................................' num2str(inflation_eq,'%04.4f')  ]);
        fprintf(resultsfile,['\n']);
        fprintf(resultsfile,['4.Output........................................' num2str(output_eq,'%04.4f')  ]);
        fprintf(resultsfile,['\n']);
        fprintf(resultsfile,['5.Employment....................................' num2str(employment_eq,'%04.4f')  ]);
        fprintf(resultsfile,['\n']);
        fprintf(resultsfile,['6.Nominal Investment Output.....................' num2str(nipa_bank_ssss.IY_nom(1),'%04.4f')  ]);
        fprintf(resultsfile,['\n']);
        fprintf(resultsfile,['7.Labor Share...................................' num2str(nipa_bank_ssss.inc_share(1,1),'%04.4f')  ]);
        fprintf(resultsfile,['\n']);
        fprintf(resultsfile,['8.Debt to income................................' num2str(nipa_bank_ssss.debt_inc(1),'%04.4f')  ]);
        fprintf(resultsfile,['\n']);
        fprintf(resultsfile,['9.Bequests to income............................' num2str(nipa_bank_ssss.beq_inc(1),'%04.4f')  ]);
        fprintf(resultsfile,['\n']);
        fprintf(resultsfile,['10.Average Return...............................' num2str(nipa_bank_ssss.average_return(1),'%04.4f')  ]);
        fprintf(resultsfile,['\n']);
        
        
    end
    
     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Close result text file                           %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        fclose(resultsfile);
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Save Workspace                                   %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    save(strcat('matlab_results/results_',controller.savename))
    
end

