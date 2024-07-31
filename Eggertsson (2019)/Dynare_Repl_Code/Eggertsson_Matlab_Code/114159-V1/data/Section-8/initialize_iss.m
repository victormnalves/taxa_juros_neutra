function [ ps,gov,prices,economy ] = initialize_iss(ps,prod,gov,policy,economy)

% Initializes variables

  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Economywide Variables    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

    economy.ag.K = 100 * ones(policy.nt + 2,1)  ; % aggregates. Initial guesses. 
    economy.ag.L = 20 * ones(policy.nt + 2,1);
    economy.ag.Y = 25*ones(policy.nt + 2,1);
    economy.ag.C = zeros(policy.nt+2,1);
    economy.ag.K_supply = zeros(policy.nt+2,1);
    
    economy.ag.profit = zeros(policy.nt+2,1);

    % Total Wage income in Economy
    economy.ag.total_wag = zeros(policy.nt+2,1);
    
    % Tax Bases
    economy.ag.rbase = zeros(policy.nt+2,1);
    economy.ag.rbasesq = zeros(policy.nt+2,1);
    economy.ag.wbase = zeros(policy.nt+2,1);
    economy.ag.wbasesq = zeros(policy.nt+2,1);
    economy.ag.ybase = zeros(policy.nt+2,1);
    economy.ag.ybasesq = zeros(policy.nt+2,1);
    
    % LSRA Vars
    economy.ag.vdebt = zeros(policy.nt+2,1);
    economy.ag.vc = zeros(policy.nt+2,1);
    economy.ag.xg = 0;
    
    % Transfer Related Payments for NIPA
    economy.ag.personal_transfer_payments = zeros(policy.nt+2,1);
    economy.ag.personal_transfer_receipts = zeros(policy.nt+2,1);
    
    
    % Bequests
    
    economy.ag.brt = zeros(policy.nt+2,1);
    economy.ag.brt_i = zeros(policy.nt+2,1);
    economy.ag.brt_u = zeros(policy.nt+2,1);
    economy.ag.bgt = zeros(policy.nt+2,1);
    economy.ag.bgt_i = zeros(policy.nt+2,1);
    economy.ag.bgt_u = zeros(policy.nt+2,1);
    
    
%%%%%%%%%%%%%%%%%%%%%
% Prices Initialization
%%%%%%%%%%%%%%%%%%%%%

    prices.mpk(1) = prod.epsilon(1)*economy.ag.A(1)*(economy.ag.K(1)/(economy.ag.L(1)))^(prod.epsilon(1)-1);
    prices.mpk = prices.mpk(1)*ones(policy.nt+2,1);
    prices.rentk = prices.mpk;
    
    prices.r = prices.mpk - prod.deprec;
    prices.wages = ones(policy.nt+2,1);
    
    prices.r(1) =  .02;
    prices.mpk(1) = .02 + prod.deprec(1);
    
    %prices.mpk(1) = -.0002 + prod.deprec(1);
    %prices.wages(1) = 1.0966;

for type=1:economy.num_types
    
        %%%%%%%%%%%%%%%%%%%%%
        % Optimal Prices
        %%%%%%%%%%%%%%%%%%%%%
        
        % Interest rates will vary by individual because of the gap between
        % the borrowing and the lending rate
        ps(type).prices.r = repmat(prices.r,1,ps(type).demog.lifespan);
        
        ps(type).prices.r_annuity = repmat(prices.r,1,ps(type).demog.lifespan);
        
        % Shadow prices -- this is a little bit of an anachronism now
        ps(type).opt.sp = zeros(policy.nt+2,ps(type).demog.lifespan);
        % Leisure Supply decisions
        ps(type).opt.l = .5*ones(policy.nt+2,ps(type).demog.lifespan); % Optimal Leisure
        
        % If endogenous labor, leisure is equal to zero
        if ps(type).endogenous_labor == 0 
            
            ps(type).opt.l = 0*ones(policy.nt+2,ps(type).demog.lifespan);
            
        end 
        
        
        ps(type).opt.sw = zeros(policy.nt+2,ps(type).demog.lifespan); % shadow wages
        ps(type).opt.C = zeros(policy.nt+2,ps(type).demog.lifespan); % Consumption
        ps(type).opt.a = 0*ones(policy.nt+2,ps(type).demog.lifespan); % assets at beginning of period
        ps(type).opt.inc = zeros(policy.nt+2,ps(type).demog.lifespan); % wage income
        ps(type).opt.cap_inc = zeros(policy.nt+2,ps(type).demog.lifespan); % capital income
        ps(type).opt.profit = zeros(policy.nt+2,ps(type).demog.lifespan); % profits
         
        ps(type).opt.br = zeros(policy.nt+2,ps(type).demog.lifespan); % bequest received
        ps(type).opt.bg = zeros(policy.nt+2,ps(type).demog.lifespan); % bequest given
        ps(type).opt.bgo = zeros(policy.nt+2,ps(type).demog.lifespan);
        ps(type).opt.bro = zeros(policy.nt+2,ps(type).demog.lifespan);
        ps(type).opt.bg_adj = zeros(policy.nt+2,ps(type).demog.lifespan);
        ps(type).opt.br_adj = zeros(policy.nt+2,ps(type).demog.lifespan);
        
        ps(type).opt.brt = zeros(policy.nt+2,1);
        ps(type).opt.brt_i = zeros(policy.nt+2,1);
        ps(type).opt.brt_u = zeros(policy.nt+2,1);
        ps(type).opt.bgt = zeros(policy.nt+2,1);
        ps(type).opt.bgt_i = zeros(policy.nt+2,1);
        ps(type).opt.bgt_u = zeros(policy.nt+2,1);

        
        % Social Security
        ps(type).opt.ben = zeros(policy.nt+2,ps(type).demog.lifespan); % social security benefits
        ps(type).opt.aime = zeros(policy.nt + 2,1); % annual index of monthly earnings, for social security
        
        % Taxes
        % Will have nt+2 rows; the first row is the initial SS
        % The final row is the final ss
        
        
        ps(type).tax.ya = 0*ones(policy.nt+2,ps(type).demog.lifespan); % average income
        ps(type).tax.ym = 0*ones(policy.nt+2,ps(type).demog.lifespan); % marginal income

        ps(type).tax.wa = zeros(policy.nt+2,ps(type).demog.lifespan); % average wage
        ps(type).tax.wm = zeros(policy.nt+2,ps(type).demog.lifespan); % marginal wage

        ps(type).tax.ka = zeros(policy.nt+2,ps(type).demog.lifespan); % capital average
        ps(type).tax.km = zeros(policy.nt+2,ps(type).demog.lifespan); % capital marginal

        ps(type).tax.sst = gov.rep(1)*(2/9)*ones(policy.nt+2,ps(type).demog.lifespan); % social security

        % Population
        
        ps(type).demog.pop = ones(policy.nt + 2,ps(type).demog.lifespan); % population 
        
        %%%%%%%%%%%%%%%%%%%%%
        % LSRA 
        %%%%%%%%%%%%%%%%%%%%%
        
            ps(type).opt.vcur = zeros(ps(type).demog.lifespan,1); % transfers to initial generations
            ps(type).opt.V = zeros(policy.nt+2,1); % transfers to future generations
            ps(type).opt.usave = zeros(policy.nt+2,1); % utility of initial steady state
            ps(type).opt.ubar = 0; % utility of future generations

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Government Initialization  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Total social security taxes
    gov.total_ss = zeros(policy.nt + 2,1);
    
    % Total social security payments
    gov.total_ss_pay = zeros(policy.nt+2,1);

    % Social Security Taxes
    gov.total_ss_tax = zeros(policy.nt+2,1);
    
    % Social Security Tax Rate
    gov.tax.sst = zeros(policy.nt + 2,1);
    
    % Tax Aggregates -- Fill in later
    gov.tax.resid_req = zeros(policy.nt + 2,1);
    gov.tax.exog_rev = zeros(policy.nt + 2,1);
    gov.tax.tax_rev_req = zeros(policy.nt + 2,1);
    gov.tax.revt = zeros(policy.nt + 2,1);
    
  

end

