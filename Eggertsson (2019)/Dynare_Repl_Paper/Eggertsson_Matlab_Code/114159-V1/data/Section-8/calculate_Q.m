function [ps,prices,economy ] = calculate_Q(ps,prod,gov,policy,economy,prices,tm)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculates Tobin's Q                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for year =tm.year_min:tm.year_max % Loop through the years
    
    % Calculate Tobin's Q; depends on whether in a steady state or
    % transition

    if (year==1 || year==(policy.nt+2)) % Initial or final steady stat
        
         % Essentially equation 3.17, but I need to work on the taxes...
         % I'm not sure if we can have Q if there are multiple types of
         % agents -- will have to do some derivations
         economy.Q(year) = (1-gov.z(year)*(gov.tax.yp(year) + gov.tax.kp(year))) + (1-(gov.tax.yp(year) + gov.tax.kp(year)))*prod.b(year)*( (economy.ag.K(year)*economy.n(year) )/economy.ag.K(year)) ; 
    
    else
        
        % Essentially equation 3.17, but note here we will use marginal taxes
        % Need to look into the taxes stuff
    
        % Will need to update this Q later
        % economy.Q(year) = (1-gov.z(year)*(gov.tax.yp(year) + gov.tax.kp(year))) + (1-(gov.tax.yp(year) + gov.tax.kp(year)))*prod.b(year)*(economy.ag.I(year)/economy.ag.K(year)) ; 
        economy.Q(year) = 1;
        
    end
    

    % Interest Rate for borrowers
%     if (year == 1 || year==(policy.nt + 2))
%         
%         % This is equation 3.19 on pg 38. In the steady state, qtp1=qt
%         % Need to think about depreciation
%         prices.r(year) = (prices.mpk(year))/economy.Q(year) - prod.deprec(year);
%     
%     else
%         
%         prices.r(year) = (prices.mpk(year))/economy.Q(year) - prod.deprec(year);
%         
%     end
    
 
end






end

