function [ pv] = getpv1(ps,economy,prod,gov,policy,prices,initial_year,initial_age)

% Gets the present value of income for a particular individual
% who is an initial age 'age' in year 'year' 


% Initialize present value and the discount
pv = 0;

if initial_age > 1
    
    pv = pv + ps.opt.a(1,initial_age)*(1+prices.r(initial_year)*(1-ps.tax.ya(initial_year,initial_age) - ps.tax.ka(initial_year,initial_age))) ;
    
    
    pv = pv + ps.opt.vcur(initial_age);
    
end

if initial_age==1 && initial_year > 1
    
    pv = pv + ps.opt.V(initial_year);
    
end

discount = 1;

% Loop thorugh ages, beginning at initial age
for age=initial_age:ps.demog.lifespan
    
    % 'year_index' is an index corresponding to the year at which
    % an individual who is age 'initial_age' at 'initial_yearyear' is 
    % age 'age'
    
    
    year_index = initial_year + age - initial_age;
        
    % If it is the initial steady state, the year index is 1 only!
    if initial_year==1
        
        % Now, we need an adjustment factor as well, for the fact that
        % wages are growing throughout your life.
        
        year_index = 1;
        
    end
    
    % If beyond final ss, set year index to fss
    if year_index > (policy.nt+2)
        
        year_index = min(year_index,policy.nt+2);
    
    end
    
    % This isn't in the original KA code
    year_born = year_index - age + 1;

    if year_born < 1
        year_born = 1;
    end
    
     % Add in bequests -- note that we use the discount from last year, as
     % in AK. 
     if age== ps.util.a_br(year_born)
         
         pv = pv + ps.prod_af(age)*ps.opt.br(year_index,age)/discount;
         
     end
    
    % We don't discount the initial year
    if age > initial_age
        
        % Update Discount; we use average tax rates for the discount

        discount = discount*(1+prices.r(year_index) * ...
                  (1-ps.tax.ka(year_index,age)-ps.tax.ya(year_index,age)));
              
                      
    end
       
    % Update the present value of earnings
    pv = pv + ps.prod_af(age) * ps.opt.wstar(year_index,age) * ...
              (1-ps.tax.ya(year_index,age) - ps.tax.wa(year_index,age) )...
              /discount;
          
          
     % Add in social security income
     pv = pv + ps.opt.ben(year_index,age)/discount; 
     
    

    
     % Subtract social security taxes
     pv = pv - ps.tax.sst(year_index,age)*ps.prod_af(age)*ps.opt.inc(year_index,age)/discount;
     

end




end

