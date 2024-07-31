function [ obj_val] = sec_stag_calibration(param_vector,other_params,moments)

% This function will be minimized in order to choose the parameters to
% match the moments for the data

% It calibrates the model to the 2015 US economy

% Retrieve elements from the parameter vector

    delta = param_vector(1);
    epsilon = param_vector(2);
    dlim_pct = param_vector(3);
    theta = param_vector(4);
    mu = param_vector(5);
    
% Retrieve parameters which are not estimated

    gamma = other_params.gamma;
    deprec = other_params.deprec;
    sigma = other_params.sigma;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set up non calibrated parameters                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    % Number of transition periods
    
    policy.nt = 150;
    
    % Whether there is a lsra (lump sum redistribution authority -- see
    % some of the discussion in Auerbach & Kotlikoff)
    
    policy.lsra = 0;
    
    % year of transition
    
    policy.year_implemented = 2;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set up the Common Person Structure  %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % While there will be several people, they will share a common
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
        p1.util.dlim_pct = dlim_pct*ones(policy.nt + 2,1); % debt limit as a percent of income. Only valid when p1.dlim = 1. Note that this stores the dlim in the "year born" fashion. 
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
            
            % Government Spending as JR Intended :)
            
                % We here set government spending as a % of GDP
                % We also set the Debt as % of GDP (will come later)
            
                gov.spend.jr = 1;
                
                gov.spend.amt = ones(policy.nt+2,1); % Government spending. Simply initializing the variable.
                
                gov.debt = 0*ones(policy.nt+2,1); % Debt as a percentage of the capital stock. 
                gov.debt_alt = 0*ones(policy.nt+2,1); % Gross debt. 
                gov.deficit = 0*ones(policy.nt+2,1); % Deficit. Since spend.jr=1, this doesn't matter
                
            % Proportions for Endogenous TAxes
            
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

            % This sets the TFR for the dates -- note this is different
            % from the above
            
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
            
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Relative Price of Investment Goods  and Capital Specific Productivity
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        
        relp_1970 = 1.3;
        
        relp_2015 = 1;
        
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
    
        % For this exercise, set NFA to zero 
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
        
        theta_1970 = 7.15*ones(policy.nt+2,1);
        
        theta_2015 =  theta*ones(policy.nt+2,1);
        
        theta_transition_vector = [theta_1970(1):-.05:theta_2015(1)]';
        
        theta_transition = theta_2015;
        
        theta_transition(1:length(theta_transition_vector)) = theta_transition_vector;
        
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set up the Run Schedule                   %%
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        run_schedule.iss.do = 1; % Do the ISS
        run_schedule.fss.do = 0; % Calculate FSS
        run_schedule.trans.do = 0; % Calculate transition
        run_schedule.pe.run=0;
        
        run_schedule.tol = 1e-3; % Tolerance
        run_schedule.maxiter = 75; % Max # of iterations
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
    % Run it                                                            %%
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
            
                prod.theta = theta_1970;
                
                prod.monop_comp= 1; 
                
            % Specify Relative Price
            
                economy.relp = relp_1970*ones(policy.nt+2,1);
                
      
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Reset to 2015 Values                      %%
        % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Set productivity to 2015 vectors

                economy.ag.AL_growth = AL_growth_2015;
                economy.ag.AL = AL_2015;

                economy.ag.AL_growth_iss = AL_growth_iss_2015;
                economy.ag.AL_growth_fss = AL_growth_fss_2015;

            % Set num kids to 2015

                p1.demog.num_kids = num_bb_2015;

            % Set the pop growth to 2015

                economy.n = bb_vec_2015;

            % Set mortality / survival
            
                p1.demog.s = s_2015; % survival
                p1.demog.ms = ms_2015; % marginal survival; column i = prob of surviving from i to j
                p1.demog.mhs = mhs_2015;
                
            % Specify debt the way JR Intended

                gov.debt_gdp = gov_debt_gdp_2015;
                
            % Set government spending
            
                gov.spend.amt_gdp = gov_amt_gdp_2015;
                
            % Specify Net Foreign Assets
            
                economy.nfa_gdp = nfa_gdp_2015;
                
            % Specify Monopolistic Competition
            
                prod.theta = theta_2015;
                
            % Specify Relative Price
                
                % Let's change this up a little bit
                economy.relp = relp_2015*ones(policy.nt+2,1);
        
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Secular Stagnation Exercise               %%
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % We will take some values from this
            [ps_ss,gov_ss,economy_ss,prices_ss,nipa_bank_ss] = master_control(p1,prod,gov,policy,economy,run_schedule);   
            
            % Generate the objective function
            obj_val = 5000*(prices_ss.r(1) - moments.r)^2 + 1000*(nipa_bank_ss.IY(1) - moments.IY)^2 + 1000*(.01*(nipa_bank_ss.inc_share(1,1) - moments.ls))^2 + 1000*(nipa_bank_ss.debt_inc(1) - moments.debt_inc)^2 + 1000*(nipa_bank_ss.beq_inc(1) - moments.beq_inc)^2;

end

