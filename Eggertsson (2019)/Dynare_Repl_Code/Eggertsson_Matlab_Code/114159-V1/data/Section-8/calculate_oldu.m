function [ps] = calculate_oldu(ps,policy,economy)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate Old Utility
%
% A lot of this code follows the original Auerbach & Kotlikoff codes
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1: Get Old Utility from individuals born before the change in policy %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Loop through types
for type = 1:economy.num_types 
    
    % Calculate utility of individuals in initial steady state
    ps(type).opt.bu = ucalc(ps(type),policy,1,ps(type).demog.lifespan,1);
       
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 2: Calculate Utilites of people alive at year 2     %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% We calculate their pre-intervention announcement utilies
% Their post intervention utilities
% and their total utilities

% Loop through types
for type = 1:economy.num_types 
    
    % Calculate utility of individuals in initial steady state
    u_pre = zeros(ps(type).demog.lifespan,1);
    u_post = zeros(ps(type).demog.lifespan,1);
    
    for age=2:ps(type).demog.lifespan
        
        % Calculate the utility of individuals pre announcement
        u_pre(age) = ucalc(ps(type),policy,1,age-1,1);
        
        % Calculate the utility of individuals post announcement
        u_post(age) = ucalc(ps(type),policy,age,ps(type).demog.lifespan,2);
        
    end
    
    u_pregen = u_pre + u_post;
    
    ps(type).opt.u_pre = u_pre;
    ps(type).opt.u_post = u_post;
    ps(type).opt.u_pregen = u_pregen;
    
end

end

