function [ mult] = r_mult_fin(ps,policy,initial_age,initial_year,final_age)

mult = 1;

for age=(initial_age + 1):final_age
    
    year_index = initial_year + age - initial_age;

    if initial_year==1

        year_index = 1;

    end

    if year_index > ( policy.nt+2 )
        
        year_index = policy.nt+2;
        
    end
    
    mult = mult*(1+ps.prices.r(year_index,age)*(1-ps.tax.ya(year_index,age) - ps.tax.ka(year_index,age)) );
    
end



end
