function [ ps ] = opt_control(ps,gov,prod,policy,prices,economy,run_schedule,init_year_min,init_year_max)

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop through Types     %
%%%%%%%%%%%%%%%%%%%%%%%%%%


for type = 1:economy.num_types
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Loop through Years     %
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % This indexes people born during this year,
    % or during the transition, people whose optimizaiton begins during
    % this year
    
    for year=init_year_min:init_year_max
        
       
        % For the first year of the transition, we loop through all ages
        % allive
        
        if year==2
            
            for age=1:ps(type).demog.lifespan
            
                ps = opt_lb(ps,gov,prod,policy,prices,economy,run_schedule,year,age,init_year_min,init_year_max);
            
                test=[0.556451914447317]*(1+prices.r(1)*(1-.15));
            end
        
        else % Else, we only do for people born at the beginning of the year
            
            ps = opt_lb(ps,gov,prod,policy,prices,economy,run_schedule,year,1,init_year_min,init_year_max);
            
        end
        
        
        
    end % Year
        
end % Type

end

