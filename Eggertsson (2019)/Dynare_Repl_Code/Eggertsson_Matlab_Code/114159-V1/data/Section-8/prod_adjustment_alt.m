function [ adjustment ] = prod_adjustment(initial_year,age,initial_age,year_min,year_max,economy,policy)

% Adjusts Productivity

% If year_min==year_max==1, we are calculating initial steady state, and we
% need to adjust productivity in the usual way

year_index = initial_year + age - initial_age;

if year_min==year_max && year_min==1
    
    adjustment = (1+economy.lprod_growth_iss)^(age - 1);
    
elseif year_min==year_max && year_min==(policy.nt + 2)
    
    adjustment = (1+economy.lprod_growth_fss)^(age-1);
    
elseif year_index > (policy.nt+2)
    
    adjustment = (1+economy.lprod_growth_fss)^(year_index - policy.nt+2);
    
elseif year_index < 1
    
    adjustment = (1+economy.lprod_growth_iss)^(year_index - 1);
    
else
    
    adjustment = 1;


end

