function [ ps,gov,prices,economy ] = initialize_fss(ps,gov,prices,policy,economy)

% Initializes variables

%%%%%%%%%%%%%%%%%%%%%
% Prices Initialization
%%%%%%%%%%%%%%%%%%%%%

    prices.r(policy.nt+2) = prices.r(1);
    prices.mpk(policy.nt+2) = prices.mpk(1);
    prices.rentk(policy.nt+2) = prices.rentk(1);
    
    prices.wages(policy.nt+2) = prices.wages(1)*(1+economy.ag.AL_growth_fss)^(policy.nt + 1)*(1+economy.ag.A_growth_adj_fss)^(policy.nt + 1);
    
%%%%%%%%%%%%%%%%%%%%%
% Individual Variable Initialization
%%%%%%%%%%%%%%%%%%%%%

for type=1:economy.num_types
    

        ps(type).opt.l(policy.nt+2,:) = ps(type).opt.l(1,:);
        ps(type).opt.sw(policy.nt+2,:) = ps(type).opt.sw(1,:);
        ps(type).opt.C(policy.nt+2,:) = ps(type).opt.C(1,:)*(1+economy.ag.AL_growth_fss)^(policy.nt + 1)*(1+economy.ag.A_growth_adj_fss)^(policy.nt + 1);
        ps(type).opt.a(policy.nt+2,:) = ps(type).opt.a(1,:)*(1+economy.ag.AL_growth_fss)^(policy.nt + 1)*(1+economy.ag.A_growth_adj_fss)^(policy.nt + 1);
        ps(type).opt.inc(policy.nt+2,:) = ps(type).opt.inc(1,:)*(1+economy.ag.AL_growth_fss)^(policy.nt + 1)*(1+economy.ag.A_growth_adj_fss)^(policy.nt + 1);
        ps(type).opt.cap_inc(policy.nt+2,:) = ps(type).opt.cap_inc(1,:)*(1+economy.ag.AL_growth_fss)^(policy.nt + 1)*(1+economy.ag.A_growth_adj_fss)^(policy.nt + 1);
        
        % Social Security
        
        ps(type).opt.ben(policy.nt+2,:) = ps(type).opt.ben(1,:)*(1+economy.ag.AL_growth_fss)^(policy.nt + 1)*(1+economy.ag.A_growth_adj_fss)^(policy.nt + 1);
        
        ps(type).opt.aime(policy.nt+2,:) = ps(type).opt.aime(1,:)*(1+economy.ag.AL_growth_fss)^(policy.nt + 1)*(1+economy.ag.A_growth_adj_fss)^(policy.nt + 1);
        
        % Taxes
        % Will have nt+2 rows; the first row is the initial SS
        % The final row is the final ss
        

        ps(type).tax.ya(policy.nt+2,:) = ps(type).tax.ya(1,:);
        ps(type).tax.ym(policy.nt+2,:) = ps(type).tax.ym(1,:);

        ps(type).tax.wa(policy.nt+2,:) = ps(type).tax.wa(1,:);
        ps(type).tax.wm(policy.nt+2,:) = ps(type).tax.wm(1,:);

        ps(type).tax.ka(policy.nt+2,:) = ps(type).tax.ka(1,:);
        ps(type).tax.km(policy.nt+2,:) = ps(type).tax.km(1,:);

        ps(type).tax.sst(policy.nt+2,:) = ps(type).tax.sst(1,:);

        

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Economywide Variables    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

    economy.ag.K(policy.nt+2) =  economy.ag.K(1)*(1+economy.n(policy.nt+2))^(policy.nt + 1)*(1+economy.ag.AL_growth_fss)^(policy.nt + 1)*(1+economy.ag.A_growth_adj_fss)^(policy.nt + 1)  ; 

    economy.ag.L(policy.nt+2) = economy.ag.L(1)*(1+economy.n(policy.nt+2))^(policy.nt + 1);


end

