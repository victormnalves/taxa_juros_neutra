function [ s] = gens(ms)

grid_size = size(ms);

num_periods = grid_size(1);
lifespan = grid_size(2);

% First, Calculate Total Survival for First Period

s = zeros(num_periods,lifespan);
s(:,1) = ms(1,1);
for age=2:lifespan
    
    s(1,age) = s(1,age-1)*ms(1,age-1);
    
end

% Now Calculate Survival for Rest of Periods
for t=2:num_periods
    
    for age = 2:lifespan
        
      s(t,age) = s(t-1,age-1)*ms(t-1,age-1);
        
    end
    
    
    
end

% Devide everything so that survival is 1 at age 1
s = s/s(1,1);

end

