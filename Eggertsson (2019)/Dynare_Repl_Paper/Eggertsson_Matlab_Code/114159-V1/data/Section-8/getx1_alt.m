function [ zxi, x1, nu,H] = getx1_alt(ps,policy,prices,economy,tm,rec)

% Calculates the percentage of the present value of income at age beginning_age
%, GETX1. 
% Also, calculate consumption in period as a ratio to consumption in period
% initial age, ZXI


% Even though we may not start from age 1, we still want these vectors
% to be the size of the lifespan

nu = zeros(ps.demog.lifespan,1);
H = zeros(ps.demog.lifespan,1);
zxi = zeros(ps.demog.lifespan,1);
x1 = 0;

% Discount rates
s = 1; % marginal taxes discount
sa = 1; % average tax rates discount

bequest_discount = 1; % Sometimes we need to keep track of the discount during the year of the bequest

for age=rec.beginning_age:rec.final_age % from initial age to final_age

    % update year index
    tm.age = age;
    [year_index,year_born] = create_index(tm,policy);

    % Create year that indexes the beginning year
    if age==rec.beginning_age
        
        beginning_year = year_index;
    end

    % The 'nu' equation 3.11 pg 32 AK
    %nu(age) = (1+ps.util.alpha(year_index) * ps.util.rho(year_index) * ( ps.opt.wstar(year_index,age) * (1 - ps.tax.ym(year_index,age) - ps.tax.wm(year_index,age) ) )^(1-ps.util.rho(year_index)))^((ps.util.rho(year_index)-ps.util.gamma(year_index))/(1-ps.util.rho(year_index)));
    
    % Alt Calculation done by Auerbach & Kotlikoff's old fortran code. It
    % turns out that the book has a small mistake in it -- it should be
    % alpha to the rho, not alpha times rho.
    
    if ps.endogenous_labor ==1
        
        nu(age) = (1+ps.util.alpha(year_born)*(prod_adjustment(tm,age,economy,policy)*ps.opt.wstar(year_index,age)*(1-ps.tax.ym(year_index,age) - ps.tax.wm(year_index,age) )/ps.util.alpha(year_born))^(1-ps.util.rho(year_born)))^((ps.util.rho(year_born)-ps.util.gamma(year_born))/(1-ps.util.rho(year_born)));
        H(age) = nu(age)^(1/ps.util.gamma(year_born));

    else
        
        nu(age) = 1;
        H(age) = 1;
    end
    % Calculate Discount rates; we don't discount first year
    if age~=rec.beginning_age
        
         % Marginal Discount rate
         s = s*(1+ps.prices.r(year_index,age)*(1-ps.tax.ym(year_index,age) - ps.tax.km(year_index,age) ) );
         
         % Average discount rate
         sa = sa*(1+ps.prices.r(year_index,age)*(1-ps.tax.ya(year_index,age) - ps.tax.ka(year_index,age)) );
         
         if age==ps.util.a_bg(year_born)
             
             bequest_discount = sa;
             
         end
    end
    
    if age==rec.beginning_age
        
        zxi(age) = 1;
        
    else
        
        % Note the age-initial age for the discount
        
        % ZXI depends on whether or not there is exogenous labor
        
        % Adjust the survival adjustment
        
        
        if ps.endogenous_labor ==1
            
            zxi(age) = (s*(ps.demog.s(year_index,age)/ps.demog.s(beginning_year,rec.beginning_age))/(1+ps.util.delta(year_born))^(age-rec.beginning_age))^(ps.util.gamma(year_born))*(nu(age)/nu(rec.beginning_age));
            
        else
            
            %zxi(age) = (s*(ps.demog.s(year_index,age)/ps.demog.s(tm.initial_year,tm.initial_age))/(1+ps.util.delta(year_born))^(age-rec.beginning_age))^(ps.util.gamma(year_born));
            
            % Perhaps it should be beginning age, not initial age
            % In addition, we should devide by beginning year, I think :)
            
            zxi(age) = (s*(ps.demog.s(year_index,age)/ps.demog.s(beginning_year,rec.beginning_age))/(1+ps.util.delta(year_born))^(age-rec.beginning_age))^(ps.util.gamma(year_born));

        end
    end
    
    if ps.endogenous_labor ==1
        
        % make a warning if there is wages are zero
        if ps.opt.wstar(year_index,age)==0
            
            warning('Are you sure you want endogenous labor? We have that wstar is zero')
 
        end
        
        x1 = x1 + zxi(age)*(1+(prod_adjustment(tm,age,economy,policy)*ps.opt.wstar(year_index,age)*(1-ps.tax.ym(year_index,age) - ps.tax.wm(year_index,age))/ps.util.alpha(year_born))^(-ps.util.rho(year_born))*(1-ps.tax.ya(year_index,age) - ps.tax.wa(year_index,age)  )*prod_adjustment(tm,age,economy,policy)*ps.opt.wstar(year_index,age))/sa ;
    
    else
        
        x1 = x1 + zxi(age)*(1)/sa ;

    end
    
end

% Adjust for bequests -- and remember we insert back into original budget
% constraint


% Adjust for bequests, only if we need to
if rec.final_age >= ps.util.a_bg(year_born)
x1 = x1 + zxi(ps.util.a_bg(year_born))*ps.util.mu(year_born)^(ps.util.gamma(year_born)) ...
          * ps.demog.num_kids(year_born)^(1-ps.util.gamma(year_born)) ...
          * H(ps.util.a_bg(year_born))^(-ps.util.gamma(year_born))...
          *(1/bequest_discount) ;

    
end

x1 = 1/x1;

end

