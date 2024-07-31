function [ adjustment ] = prod_adjustment(tm,age,economy,policy)

% Adjusts Productivity
% One of the most important programs
% In the year 1, sometimes we need to calculate variables from individuals
% who were born in year 0 and before
% Since we are in a steady state, we take the wage in year 1, and adjust
% for productivity

% Likewise, if the year is policy.nt +2, we need to calculate income of
% individuals who are born afterwards

% Finally, sometimes we calculate optimal decisions during the transition
% from period who are born close to policy.nt+2. When calculating wages etc
% for these people, need to adjust productivity. 

% If year_min==year_max==1, we are calculating initial steady state, and we
% need to adjust productivity in the usual way


% This may or may not be slightly confusing because it is different than
% the definition of year_index in the create_index program
year_index = tm.initial_year + age - tm.initial_age;

if tm.year_min==tm.year_max && tm.year_min==1
    
    % Note that if labor productiivty grows at rate g, then output and
    % capital also grow at g.
    
    % If Hicks Neutral grows at rate g, then output/wages/capital grow at
    % (1+g)^(1/(1-epsilon))
    adjustment = (1+economy.ag.AL_growth_iss)^(age - 1) * ...
                 (1+economy.ag.A_growth_adj_iss)^(age - 1) ;
    
elseif tm.year_min==tm.year_max && tm.year_min==(policy.nt + 2)
    
    adjustment = (1+economy.ag.AL_growth_fss)^(age-1) * ...
                 (1+economy.ag.A_growth_adj_fss)^(age - 1);
    
elseif year_index > (policy.nt+2)
    
    adjustment = (1+economy.ag.AL_growth_fss)^(year_index - (policy.nt+2)) * ...
                 (1+economy.ag.A_growth_adj_fss)^(year_index - (policy.nt+2));
    
    
elseif year_index < 1
    
    adjustment = (1+economy.ag.AL_growth_iss)^(year_index - 1) *...
                 (1+economy.ag.A_growth_adj_iss)^(year_index - 1);
else
    
    adjustment = 1;

end

