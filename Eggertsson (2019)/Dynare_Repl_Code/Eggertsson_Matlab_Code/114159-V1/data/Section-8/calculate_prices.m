function [ prices,economy ] = calculate_prices(prod,gov,economy,policy,prices,tm)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates the Prices in the Economy      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop through years from year_min to year_max, calculate interest rates
% and wages. Along the way we need to calculate Q as well. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1: Set Productivity in Year 1 so that wages are 1      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate what wage would be if no productivity adjustment, and use this
% to set the productivity adjustment to 1/wage_unadjusted

% Recall that this is simply a multiplicative factor! 
% It is not our main measure of productivity growth

if economy.adj_prod==0
    
    economy.ag.A_adj = 1;
    
end

if tm.year_min==1 && economy.adj_prod==1
    
    % Must also adjust wages for markups, so need this adjustment
    if prod.monop_comp==1
        
        markup_adj =  (prod.theta(1)-1)/prod.theta(1);
        
    else
        
        markup_adj = 1;
        
    end
    
    if prod.sigma(1)==1 % Cobb Douglas Production

        economy.ag.A_adj =  1/((1-prod.epsilon(1))*markup_adj*economy.ag.AL(1)*economy.ag.A(1)*(economy.ag.AK(1)*economy.ag.K(1))^(prod.epsilon(1))*(economy.ag.AL(1)*economy.ag.L(1))^(-prod.epsilon(1)));

        
    else
        
         omega_temp = ( prod.epsilon(1)*(economy.ag.AK(1)*economy.ag.K(1))^((prod.sigma(1)-1)/prod.sigma(1)) ...
                 + (1-prod.epsilon(1))*(economy.ag.AL(1)*economy.ag.L(1))^((prod.sigma(1)-1)/prod.sigma(1)) );
                          
        
        economy.ag.A_adj = ( omega_temp^(1/(prod.sigma(1)-1))...
                           * (1-prod.epsilon(1))*markup_adj*economy.ag.A(1)*economy.ag.AL(1)^((prod.sigma(1)-1)/prod.sigma(1))...
                           * economy.ag.L(1)^(-1/prod.sigma(1)) )^(-1);
                           
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 2: Loop through the years                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for year = tm.year_min:tm.year_max % Loop through the years
     
    
    % Aggregate Production
    
    if prod.sigma(year)==1 % Cobb douglas 
        
        economy.ag.Y(year) = economy.ag.A(year)*economy.ag.A_adj*(economy.ag.AK(year)*economy.ag.K(year))^(prod.epsilon(year))*(economy.ag.AL(year)*economy.ag.L(year))^(1-prod.epsilon(year));
    
    else % Generalized CES
        
         omega = ( prod.epsilon(year)*(economy.ag.AK(year)*economy.ag.K(year))^((prod.sigma(year)-1)/prod.sigma(year)) ...
                 + (1-prod.epsilon(year))*(economy.ag.AL(year)*economy.ag.L(year))^((prod.sigma(year)-1)/prod.sigma(year)) );
                          
         economy.ag.Y(year) = economy.ag.A(year)*economy.ag.A_adj*omega^(prod.sigma(year)/(prod.sigma(year)-1));
        
    end
    
    % Adjustment for Monopolistic Competition
    % Also Calculate Profits
    if prod.monop_comp==1
        
        prod.markup(year) = prod.theta(year)/(prod.theta(year)-1);
        
        economy.ag.profit(year) = economy.ag.Y(year) / prod.theta(year);
        
    else
        
        prod.markup(year) = 1;
        
        economy.ag.profit(year) = 0;
        
    end
    

    % Marginal Product of Capital
    % Remember, marginal product of capital has two components:
    % adjustment cots and production function
    
        if prod.sigma(year)==1
            
            prices.mpk(year) = prod.epsilon(year)*economy.ag.A(year)*economy.ag.A_adj*economy.ag.AK(year)*((economy.ag.AK(year)*economy.ag.K(year))/(economy.ag.AL(year)*economy.ag.L(year)))^(prod.epsilon(year)-1);
        
            prices.rentk(year) = (1/prod.markup(year))*prices.mpk(year);
            
        else

            
            prices.mpk(year) = omega^(1/(prod.sigma(year)-1))...
                               * prod.epsilon(year)*economy.ag.A(year)*economy.ag.A_adj*economy.ag.AK(year)^((prod.sigma(year)-1)/prod.sigma(year))...
                               * economy.ag.K(year)^(-1/prod.sigma(year));
                           
            prices.rentk(year) = (1/prod.markup(year))*prices.mpk(year);
        
        end
        
%     if tm.year_min==tm.year_max
%         
%         prices.mpk(year) = prices.mpk(year) + .5*prod.b(year)*economy.n(year)^2;
%     
%     end
%     
    % Interest Rate
    % Now includes relative price of capital goods
    
    % Year MInus 1 index
    if tm.year_min==tm.year_max

        ym1 = year;

    else

        ym1 = year - 1;

    end
    
    if (year == 1 || year==(policy.nt + 2))
        
        % This is equation 3.19 on pg 38. In the steady state, qtp1=qt
        % Need to think about depreciation
        prices.r(year) = (prices.rentk(year))/(economy.Q(year)*economy.relp(ym1)) +  (1-prod.deprec(year))*(economy.relp(year)/economy.relp(ym1)) - 1;
    
    else
        
        prices.r(year) = (prices.rentk(year))/(economy.Q(year)*economy.relp(ym1)) +  (1-prod.deprec(year))*(economy.relp(year)/economy.relp(ym1)) - 1;
    
    end
   

    % Wages
    if prod.sigma(year)==1
        
        prices.wages(year)= (1/prod.markup(year))*(1-prod.epsilon(year))*economy.ag.A(year)*economy.ag.A_adj*economy.ag.AL(year)*((economy.ag.AK(year)*economy.ag.K(year))/(economy.ag.AL(year)*economy.ag.L(year)))^(prod.epsilon(year));

    else
        
            prices.wages(year) = (1/prod.markup(year))*omega^(1/(prod.sigma(year)-1))...
                               * (1-prod.epsilon(year))*economy.ag.A(year)*economy.ag.A_adj*economy.ag.AL(year)^((prod.sigma(year)-1)/prod.sigma(year))...
                               * economy.ag.L(year)^(-1/prod.sigma(year));
        
    end
    
    % My Business
    %prices.r(year) = 0.093201;
    %prices.wages(year) = 0.786;
    
end



end

