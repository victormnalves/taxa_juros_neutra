function [ Y, r_vec, dlim, delta, gamma, mu,initial_assets_end ] = getinc_vec(ps,economy,prod,gov,policy,prices,initial_year,initial_age)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Utility Parameters        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

year_born = initial_year - initial_age + 1;
year_born = max(year_born,1);

delta = ps.util.delta(year_born);
gamma = ps.util.gamma(year_born);
mu = ps.util.mu(year_born);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial Assets            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Initial Assets, in end of period value. Thus must multiply through by
    % (1+r). 
    
    initial_assets_end = 0;

    if initial_age > 1

        initial_assets_end = initial_assets_end + ps.opt.a(1,initial_age)*(1+prices.r(initial_year)*(1-ps.tax.ya(initial_year,initial_age) - ps.tax.ka(initial_year,initial_age))) ;

        initial_assets_end = initial_assets_end + ps.opt.vcur(initial_age);

    end

    if initial_age==1 && initial_year > 1

        initial_assets_end = initial_assets_end + ps.opt.V(initial_year);

    end


Y = zeros(ps.demog.lifespan - initial_age + 1);

r_vec = ones(ps.demog.lifespan - initial_age + 1);

dlim = ones(ps.demog.lifespan - initial_age + 1);

% Loop thorugh ages, beginning at initial age
for age=initial_age:ps.demog.lifespan
    
    % Index corresponds to the place in the income vector
    index = age - initial_age + 1;
    
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
    
    
    % Tabulate Interest Rate
    if age > initial_age
        
        % Update Discount; we use average tax rates for the discount
        r_vec(index) = (1+prices.r(year_index) * ...
                       (1-ps.tax.ka(year_index,age)-ps.tax.ya(year_index,age)));
                      
    end
    
     % Add in bequests -- note that we use the discount from last year, as
     % in AK. 
     if age== ps.util.a_br(year_born)
         
         Y(index) = Y(index) + ps.prod_af(age)*ps.opt.br(year_index,age);

     end
    
       
    % Add in Earnings
    Y(index) = Y(index) +  ps.prod_af(age) * ps.opt.wstar(year_index,age) * ...
              (1-ps.tax.ya(year_index,age) - ps.tax.wa(year_index,age) );
              

     % Add in social security income  
     Y(index) = Y(index) + ps.opt.ben(year_index,age);
    
     % Social Security Taxes
     Y(index) = Y(index) - ps.tax.sst(year_index,age)*ps.prod_af(age)*ps.opt.inc(year_index,age);
     
     % Get Debt limit
     
     dlim(index) = ps.opt.dlim(year_index,age);
     
end



end

