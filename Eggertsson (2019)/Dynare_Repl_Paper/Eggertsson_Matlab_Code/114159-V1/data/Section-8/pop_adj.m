function [adj] = pop_adj(tm,economy,policy)

% Creates a year index, year born, etc
adj = 1;
if tm.year_min==tm.year_max && tm.year_min==1 % In initial year, year_index is one
    
    adj = (1+economy.lprod_growth_iss);
    
elseif tm.year_min==tm.year_max && tm.year_max==(policy.nt+2) % Final Steady State
    
    %adj = (1+economy.lprod_growth_fss);
 
end

end

