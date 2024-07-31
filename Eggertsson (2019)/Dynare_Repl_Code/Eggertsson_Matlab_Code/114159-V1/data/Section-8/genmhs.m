function [mhs] = genmhs(ms)

grid_size = size(ms);

num_periods = grid_size(1);
lifespan = grid_size(2);

% First, Calculate Total MHS for first period

mhs = zeros(num_periods,lifespan);
mhs(:,1) = 1;
for age=2:lifespan
    
    mhs(1,age) = ms(1,age-1);
    
end

% Now Calculate Survival for Rest of Periods
for t=2:num_periods
    
    for age = 2:lifespan
        
      mhs(t,age) = ms(t-1,age-1);
        
    end
    
end

end

