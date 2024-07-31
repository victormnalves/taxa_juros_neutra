function [ ps,gov,prices,economy ] = initialize_ss(ps,gov,policy,economy)

% Initializes variables

%%%%%%%%%%%%%%%%%%%%%
% Prices Initialization
%%%%%%%%%%%%%%%%%%%%%

    prices.r = 0.07*ones(policy.nt+2,1);
    prices.mpk = prices.r;
    prices.wages = ones(policy.nt+2,1);
    

for type=1:economy.num_types
    

        ps(type).opt.l = .5*ones(policy.nt+2,ps(type).demog.lifespan); % Optimal Leisure
        ps(type).opt.sw = zeros(policy.nt+2,ps(type).demog.lifespan);
        ps(type).opt.C = zeros(policy.nt+2,ps(type).demog.lifespan);
        ps(type).opt.a = 0*ones(policy.nt+2,ps(type).demog.lifespan);
        ps(type).opt.inc = zeros(policy.nt+2,ps(type).demog.lifespan);
        ps(type).opt.cap_inc = zeros(policy.nt+2,ps(type).demog.lifespan);
        
        % Social Security
        
        ps(type).opt.ben = zeros(policy.nt+2,ps(type).demog.lifespan);
        
        ps(type).opt.aime = zeros(policy.nt + 2,1);
        
        % Taxes
        % Will have nt+2 rows; the first row is the initial SS
        % The final row is the final ss
        

        ps(type).tax.ya = .15*ones(policy.nt+2,ps(type).demog.lifespan);
        ps(type).tax.ym = .15*ones(policy.nt+2,ps(type).demog.lifespan);

        ps(type).tax.wa = zeros(policy.nt+2,ps(type).demog.lifespan);
        ps(type).tax.wm = zeros(policy.nt+2,ps(type).demog.lifespan);

        ps(type).tax.ka = zeros(policy.nt+2,ps(type).demog.lifespan);
        ps(type).tax.km = zeros(policy.nt+2,ps(type).demog.lifespan);

        ps(type).tax.sst = gov.rep(1)*(2/9)*ones(policy.nt+2,ps(type).demog.lifespan);

        % Population
        
        ps(type).demog.pop = ones(policy.nt + 2,ps(type).demog.lifespan);
        
        %%%%%%%%%%%%%%%%%%%%%
        % LSRA 
        %%%%%%%%%%%%%%%%%%%%%
        
        if policy.lsra==1
            
            ps(type).opt.vcur = zeros(ps(type).demog.lifespan,1);
            
            ps(type).opt.V = zeros(policy.nt+2,1);
            
            ps(type).opt.usave = zeros(policy.nt+2,1);
            
            ps(type).opt.ubar = 0;

        end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Government Initialization  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    gov.total_ss = zeros(policy.nt + 2,1);
    
    % Total social security payments
    gov.total_ss_pay = zeros(policy.nt+2,1);

    % Social Security Benefits
    gov.total_ss_tax = zeros(policy.nt+2,1);
    
    % Social Security Tax Rate
    gov.tax.sst = zeros(policy.nt + 2,1);
    
    % Tax Aggregates
    
    gov.tax.resid_req = zeros(policy.nt + 2,1);
    gov.tax.exog_rev = zeros(policy.nt + 2,1);
    gov.tax.tax_rev_req = zeros(policy.nt + 2,1);
    gov.tax.revt = zeros(policy.nt + 2,1);
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Economywide Variables    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Total Wage income in Economy
    economy.ag.total_wag = zeros(policy.nt+2,1);
    
    economy.ag.rbase = zeros(policy.nt+2,1);
    economy.ag.rbasesq = zeros(policy.nt+2,1);
    economy.ag.wbase = zeros(policy.nt+2,1);
    economy.ag.wbasesq = zeros(policy.nt+2,1);
    economy.ag.ybase = zeros(policy.nt+2,1);
    economy.ag.ybasesq = zeros(policy.nt+2,1);
    
    economy.ag.C = zeros(policy.nt+2,1);

    economy.ag.vdebt = zeros(policy.nt+2,1);
    
    economy.ag.vc = zeros(policy.nt+2,1);
    
    economy.ag.xg = 0;

end

