function [U] = ucalc(ps,policy,age_min,age_max,year_init)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate Utility from age agemin to age agemax
% % for individuals who are agemin in year_init
%
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1: Calculate utility              %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

U = 0;

% Calculate Discount Rates
% Thus here we allow beta to vary throughout an individual's life span
discount = ones(ps.demog.lifespan,1); 

for age=2:ps.demog.lifespan
    
    par_index = year_init - age_min + 1 + age - 1;
    
    if year_init==1
        
        par_index = 1;
        
    end
    
    % If we are calculating people who are born in the early years, their
    % parameters are the steady state ones
    if par_index <= 1
        
        par_index = 1;
        
    end
    
    % Past the steady state, set it equal to the steady state
    
    if par_index>(policy.nt+2)
        
       par_index = (policy.nt+2);
       
    end
    
    discount(age) = discount(age-1)*(1 / (1+ ps.util.delta(par_index,age)));
    
end

% Loop Through ages
for age=age_min:age_max
    
    par_index = year_init + age - age_min ;
    
    if year_init==1
        
        par_index = 1;
        
    end
    
    % Past the steady state, set it equal to the steady state
    
    if par_index>(policy.nt+2)
        
       par_index = (policy.nt+2);
       
    end
    
    % Formula from AK 1987 pg 27
    
    u1 = ps.opt.C(par_index,age)^(1-1/ps.util.rho(par_index,age)) + ps.util.alpha(par_index,age)*ps.opt.l(par_index,age)^(1-1/ps.util.rho(par_index,age));
    
    u2 = u1^(1/(1-1/ps.util.rho(par_index,age)));
    
    % Very important the the index for discount is age!
    %u3 = (1/(1-1/ps.util.gamma(par_index,age)))*discount(age)*u2^(1-1/ps.util.gamma(par_index,age));
    
    % Use the AK definition of utility
    u3 = (1/25)*discount(age)*u2^(1-1/ps.util.gamma(par_index,age));

    U = U + u3;
    
end









end

