function [ ps,gov,prices,economy ] = initialize_trans(ps,gov,prices,policy,economy)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize Variables for Transition %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for year=2:policy.nt+2

    %%%%%%%%%%%%%%%%%%%%%
    % Prices Initialization
    %%%%%%%%%%%%%%%%%%%%%

        prices.r(year) = prices.r(1)*(1+get_growth(prices.r))^(year-1);
        prices.mpk(year) = prices.mpk(1)*(1+get_growth(prices.mpk))^(year-1);
        prices.rentk(year) = prices.rentk(1)*(1+get_growth(prices.rentk))^(year-1);
        prices.wages(year) = prices.wages(1)*(1+get_growth(prices.wages))^(year-1);
        
    %%%%%%%%%%%%%%%%%%%%%
    % Economic Aggregates Initialization
    %%%%%%%%%%%%%%%%%%%%%

        economy.ag.K(year) = economy.ag.K(1)*(1+get_growth(economy.ag.K))^(year-1);
        economy.ag.L(year) = economy.ag.L(1)*(1+get_growth(economy.ag.L))^(year-1);
        
    %%%%%%%%%%%%%%%%%%%%%
    % Individual Variable Initialization
    %%%%%%%%%%%%%%%%%%%%%

    for type=1:economy.num_types


            ps(type).opt.l(year,:) = ps(type).opt.l(1,:) .* (1 + get_growth(ps(type).opt.l)) .^(year-1) ;
            ps(type).opt.sw(year,:) = ps(type).opt.sw(1,:) .* (1 + get_growth(ps(type).opt.sw)) .^(year-1);
            ps(type).opt.C(year,:) = ps(type).opt.C(1,:) .* (1 + get_growth(ps(type).opt.C)) .^(year-1);
            ps(type).opt.a(year,:) = ps(type).opt.a(1,:) .* (1 + get_growth(ps(type).opt.a)) .^(year-1);
            ps(type).opt.inc(year,:) = ps(type).opt.inc(1,:) .* (1 + get_growth(ps(type).opt.inc)) .^(year-1);
            ps(type).opt.cap_inc(year,:) = ps(type).opt.cap_inc(1,:) .* (1 + get_growth(ps(type).opt.cap_inc)) .^(year-1);
            
            ps(type).opt.br(year,:) = ps(type).opt.br(1,:) .* (1 + get_growth(ps(type).opt.br)) .^(year-1);

            % Social Security

            ps(type).opt.ben(year,:) = ps(type).opt.ben(1,:) .* (1 + get_growth(ps(type).opt.ben)) .^(year-1);

            ps(type).opt.aime(year,:) = ps(type).opt.aime(1,:) .* (1 + get_growth(ps(type).opt.aime)) .^(year-1);

            % Taxes
            % Will have nt+2 rows; the first row is the initial SS
            % The final row is the final ss


            ps(type).tax.ya(year,:) = ps(type).tax.ya(1,:) .* (1 + get_growth(ps(type).tax.ya)) .^(year-1);
            ps(type).tax.ym(year,:) = ps(type).tax.ym(1,:) .* (1 + get_growth(ps(type).tax.ym)) .^(year-1);

            ps(type).tax.wa(year,:) = ps(type).tax.wa(1,:) .* (1 + get_growth(ps(type).tax.wa)) .^(year-1);
            ps(type).tax.wm(year,:) = ps(type).tax.wm(1,:) .* (1 + get_growth(ps(type).tax.wm)) .^(year-1);

            ps(type).tax.ka(year,:) = ps(type).tax.ka(1,:) .* (1 + get_growth(ps(type).tax.ka)) .^(year-1);
            ps(type).tax.km(year,:) = ps(type).tax.km(1,:) .* (1 + get_growth(ps(type).tax.km)) .^(year-1);

            ps(type).tax.sst(year,:) = ps(type).tax.sst(1,:) .* (1 + get_growth(ps(type).tax.sst)) .^(year-1);

    end



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Economywide Variables    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%


        economy.ag.K(year) =   economy.ag.K(1)*(1+get_growth( economy.ag.K))^(year-1)  ; 

        economy.ag.L(year) =   economy.ag.L(1)*(1+get_growth( economy.ag.L))^(year-1);


end

end

