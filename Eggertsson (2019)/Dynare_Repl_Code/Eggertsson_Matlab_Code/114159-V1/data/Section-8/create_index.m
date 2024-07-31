function [year_index,year_born,ym1_index,yp1_index,year_alt_index] = create_index(tm,policy)

% Creates a year index, year born, etc

if tm.year_min==tm.year_max && tm.year_min==1 % In initial year, year_index is one
    
    year_index = 1;
    year_born = 1;
    ym1_index = 1;
    yp1_index = 1;
    
    year_alt_index = 1;
    
elseif tm.year_min==tm.year_max && tm.year_max==(policy.nt+2) % Final Steady State
    
    year_index = policy.nt+2;
    year_born = policy.nt+2;
    ym1_index = policy.nt+2;
    yp1_index = policy.nt + 2;
    
    year_alt_index = policy.nt + 2;
    
else
    
    year_index = tm.initial_year + tm.age - tm.initial_age;
    
    % Year alt index will tell us if we should record it or not! 
    year_alt_index = year_index; 
    
    year_index = min(year_index,policy.nt + 2);
    year_index = max(year_index,1);
    
    ym1_index = year_index - 1;
    ym1_index = max(ym1_index,1);
    
    year_born = tm.initial_year - tm.initial_age + 1;
    year_born = max(year_born,1);
    
    yp1_index = year_index + 1;
    yp1_index = min(yp1_index,policy.nt + 2);
    
    

end

end

