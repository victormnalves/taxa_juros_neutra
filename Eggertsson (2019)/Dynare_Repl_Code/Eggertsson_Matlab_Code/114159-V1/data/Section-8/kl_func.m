function [ diff ] = kl_func(kl,r,prod,economy)


% This will be minimized to find the optimal kl a firm demands, given the
% interest rate r

    year = 1;

    omega = ( prod.epsilon(year)*(economy.ag.AK(year)*kl)^((prod.sigma(year)-1)/prod.sigma(year)) ...
         + (1-prod.epsilon(year))*(economy.ag.AL(year))^((prod.sigma(year)-1)/prod.sigma(year)) );

     prod.markup(year) = prod.theta(year)/(prod.theta(year)-1);

     rho = (1/prod.markup(year))*omega^(1/(prod.sigma(year)-1))...
                               * prod.epsilon(year)*economy.ag.A(year)*economy.ag.A_adj*economy.ag.AK(year)^((prod.sigma(year)-1)/prod.sigma(year))...
                               * kl^(-1/prod.sigma(year));
                           
     diff = r - (rho - prod.deprec(year));
     
     diff = diff^2;



end

