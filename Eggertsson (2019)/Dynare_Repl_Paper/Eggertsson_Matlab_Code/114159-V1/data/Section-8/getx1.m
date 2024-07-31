function [ zxi, x1, nu,H,af] = getx1(ps,policy,prices,initial_year,initial_age)

% Calculates the percentage of the present value of income at age initial_age
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

% Calculate Adjustment Discount for Adjustment Factor
% adj_disc = ones(ps.demog.lifespan,1);
% 
% for age=initial_age:ps.demog.lifespan
%     
%     % The year at which people who are 'initial_age' at year 'initial_year'
%     % are age 'age'
%     year_index = initial_year + age - initial_age;
%     
%     % If beyond final ss, set year index to fss
%     year_index = min(year_index,policy.nt+2);
%     
%     if initial_year==1
%        
%         year_index = 1;
% 
%     end
%     
%     
%     if age~=initial_age
%         
%          % Average discount rate
%          adj_disc(age) = adj_disc(age-1)*(1+prices.r(year_index)*(1-ps.tax.ya(year_index,age) - ps.tax.ka(year_index,age)) );
% 
%     end
%     
%     
% end

% Calculate Adjustment Factor
af = zeros(ps.demog.lifespan,1);
for iage=initial_age:ps.demog.lifespan
    

    af(iage) = 1;
    
    % Formula I derived
    for age = iage:ps.demog.lifespan
        
        % The year at which people who are 'initial_age' at year 'initial_year'
        % are age 'age'
        year_index = initial_year + age - initial_age;

        % If beyond final ss, set year index to fss
        year_index = min(year_index,policy.nt+2);

        if initial_year==1

            year_index = 1;

        end
        
        af(iage) = af(iage) + r_mult_fin(ps,policy,initial_age,initial_year,age)*ps.opt.sp(year_index,age);

    end
    
end

for age=initial_age:ps.demog.lifespan % from initial age to lifespan

    % The year at which people who are 'initial_age' at year 'initial_year'
    % are age 'age'
    year_index = initial_year + age - initial_age;
    
    % If beyond final ss, set year index to fss
    year_index = min(year_index,policy.nt+2);
    
    year_born = initial_year - initial_age + 1;
    year_born = max(year_born,1);
    
    if initial_year==1
       
        year_index = 1;

    end
    
    % The 'nu' equation 3.11 pg 32 AK
    %nu(age) = (1+ps.util.alpha(year_index) * ps.util.rho(year_index) * ( ps.opt.wstar(year_index,age) * (1 - ps.tax.ym(year_index,age) - ps.tax.wm(year_index,age) ) )^(1-ps.util.rho(year_index)))^((ps.util.rho(year_index)-ps.util.gamma(year_index))/(1-ps.util.rho(year_index)));
    
    % Alt Calculation done by Auerbach & Kotlikoff's old fortran code. It
    % turns out that the book has a small mistake in it -- it should be
    % alpha to the rho, not alpha times rho.
    
    if ps.endogenous_labor ==1
        
        nu(age) = (1+ps.util.alpha(year_index,age)*(ps.prod_af(age)*ps.opt.wstar(year_index,age)*(1-ps.tax.ym(year_index,age) - ps.tax.wm(year_index,age) )/ps.util.alpha(year_index,age))^(1-ps.util.rho(year_index,age)))^((ps.util.rho(year_index,age)-ps.util.gamma(year_index,age))/(1-ps.util.rho(year_index,age)));
        H(age) = nu(age)^(1/ps.util.gamma(year_index,age));

    else
        
        nu(age) = 1;
        H(age) = 1;
    end
    % Calculate Discount rates; we don't discount first year
    if age~=initial_age
        
         % Marginal Discount rate
         s = s*(1+prices.r(year_index)*(1-ps.tax.ym(year_index,age) - ps.tax.km(year_index,age) ) );
         
         % Average discount rate
         sa = sa*(1+prices.r(year_index)*(1-ps.tax.ya(year_index,age) - ps.tax.ka(year_index,age)) );
         
         if age==ps.util.a_bg(year_born)
             
             bequest_discount = sa;
             
         end
    end
    
    if age==initial_age
        
        zxi(age) = 1;
        
    else
        
        % Note the age-initial age for the discount
        
        % ZXI depends on whether or not there is exogenous labor
        
        if ps.endogenous_labor ==1
            
            zxi(age) = (s/(1+ps.util.delta(year_index,age))^(age-initial_age))^(ps.util.gamma(year_index,age))*(nu(age)/nu(initial_age));
            
            % Now, add adjustment factors
            zxi(age) = zxi(age)*(af(initial_age)/af(age))^(ps.util.gamma(year_index,age));
        
        else
            
            zxi(age) = (s/(1+ps.util.delta(year_index,age))^(age-initial_age))^(ps.util.gamma(year_index,age));
            
            % Now, add adjustment factors
            zxi(age) = zxi(age)*(af(initial_age)/af(age))^(ps.util.gamma(year_index,age));
            
        end
    end
    
    if ps.endogenous_labor ==1
        
        % make a warning if there is wages are zero
        if ps.opt.wstar(year_index,age)==0
            
            % If wstar is zero, the below will be "not a number"
            warning('Are you sure you want endogenous labor? We have that wstar is zero')
 
        end
        
        x1 = x1 + zxi(age)*(1+(ps.prod_af(age)*ps.opt.wstar(year_index,age)*(1-ps.tax.ym(year_index,age) - ps.tax.wm(year_index,age))/ps.util.alpha(year_index,age))^(-ps.util.rho(year_index,age))*(1-ps.tax.ya(year_index,age) - ps.tax.wa(year_index,age)  )*ps.prod_af(age)*ps.opt.wstar(year_index,age))/sa ;
    
    else
        
        x1 = x1 + zxi(age)*(1)/sa ;

    end
    
end

% Adjust for bequests -- and remember we insert back into original budget
% constraint


year_born = initial_year - initial_age + 1;
year_born = max(year_born,1);

x1 = x1 + zxi(ps.util.a_bg(year_born))*ps.util.mu(year_born)^(ps.util.gamma(initial_year,initial_age)) ...
          * ps.demog.num_kids(year_born)^(1-ps.util.gamma(initial_year,initial_age)) ...
          * H(ps.util.a_bg(year_born))^(-ps.util.gamma(initial_year,initial_age))...
          *(1/bequest_discount) ;

x1 = 1/x1;

end

