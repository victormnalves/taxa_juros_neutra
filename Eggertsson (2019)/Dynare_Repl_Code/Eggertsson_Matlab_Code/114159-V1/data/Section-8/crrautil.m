function [ util] = crrautil(CB,ps,tm,policy)

% CRAA Utility Function, with Bequests, No Labor

    util = 0;

    % First ,Consumption
    
    for age=1:ps.demog.lifespan

        tm.age = age;
        [year_index,year_born] = create_index(tm,policy);
        
        util = util + ps.demog.s(year_index,age)*(1+ps.util.delta(year_born))^(-(age-1))*CB(age)^(1-1/ps.util.gamma(year_born));

    end
    
    % Now, bequests
    tm.age = ps.util.a_bg(year_born);
    [year_index,year_born] = create_index(tm,policy);
    util = util + ps.demog.s(year_index,ps.util.a_bg(year_born))*(1+ps.util.delta(year_born))^(-(ps.util.a_bg(year_born)-1))*ps.util.mu(year_born)*CB(ps.demog.lifespan + 1)^(1-1/ps.util.gamma(year_born));
    
    % Multiply
    
    util = util*1/(1-1/ps.util.gamma(year_born));

    % We take negative since we will be minimizing
    
    util = -util;

end

