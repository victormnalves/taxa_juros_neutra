function [ pv] = getpv1_alt(ps,economy,prod,gov,policy,prices,tm,rec)

% Gets the present value of income for a particular individual
% Formulas mainly taken from Auerbach & Kotlikoff


% Here we make a distinction between initial age and beginning age
    % initial age is the age an individual starts to to optimize from after
    % there is a change in the policy, and the age at which the individual
    % receives his or her transfers

    % beginning age is the age we optimize from in the recursive loop

% Sometimes we come into the period with initial debt

% Sometimes we want to end the optimization with debt, so we add on final
% debt to the present value of income

% Beginning year is the year when the individual is 'beginning_age'
beginning_year = tm.initial_year + rec.beginning_age - tm.initial_age;
if tm.initial_year==1
    
    beginning_year = 1;
    
elseif beginning_year > (policy.nt + 2)
    
    beginning_year = policy.nt + 2;
    
end

% We need to multiply by (1+r) because the present value of income is taken
% at the end of period 1
pv = -rec.initial_debt*(1+ps.prices.r(beginning_year,rec.beginning_age)*(1-ps.tax.ya(beginning_year,rec.beginning_age) - ps.tax.ka(beginning_year,rec.beginning_age)));

% IF we are optimizing from beginning_age==initial_age, then there is a
% change in policy and we need to take into account the assets we have
% before, as well as transfers
if tm.initial_age > 1 && tm.initial_age==rec.beginning_age
    
    % This needs to be adjusted for productivity growth in initial steady
    % state
    % The reason is that in the steady state, optimal assets is growing at
    % a constant rate
    % Thus to get assets in period 2, update assets in period 1 by growth
    % rate
    % We also adjust by interest rate because present value is take at the
    % end of the period
    pv = pv + (1+economy.ag.AL_growth_iss)*(1+economy.ag.A_growth_adj_iss)*ps.opt.a(1,tm.initial_age)*(1+ps.prices.r(tm.initial_year,tm.initial_age)*(1-ps.tax.ya(tm.initial_year,tm.initial_age) - ps.tax.ka(tm.initial_year,tm.initial_age))) ;
    
    pv = pv + ps.opt.vcur(tm.initial_age);
    
end

% If we are at beginning_age = 1, and year > 1, we need to take into
% transfers from the authority
if tm.initial_age==1 && tm.initial_year > 1 && rec.beginning_age==1
    
    pv = pv + ps.opt.V(tm.initial_year);
    
end


discount = 1;

% Loop thorugh ages, beginning at beginning_age, and ending at final_age
for age=rec.beginning_age:rec.final_age
    
    % update year index
    tm.age = age;
    [year_index,year_born] = create_index(tm,policy);


     % Add in bequests -- note that we use the discount from last year, as
     % in AK, since the bequests are received in the beginning of the
     % period
     
     % Also, we don't add in bequests if initial_age >= bequest year! This would
     % double count bequests, since the bequests are already included in
     % the initial assets of the person above
     if age== ps.util.a_br(year_born) && tm.initial_age < ps.util.a_br(year_born)
         
         % Note that we do adjust for productivity here. What is stored in
         % br is simply the bequests received by people who are currently
         % that age. In the calculation of the ISS, bequests for the
         % individual born in period 1 will thus be significantly higher
         % than the bequests stored in ps.opt.br
         pv = pv + prod_adjustment(tm,age,economy,policy)*ps.opt.br(year_index,age)/discount;
         
     end
    
    % We don't discount the initial year
    if age > rec.beginning_age
        
        % Update Discount; we use average tax rates for the discount

        discount = discount*(1+ps.prices.r(year_index,age) * ...
                  (1-ps.tax.ka(year_index,age)-ps.tax.ya(year_index,age)));
              
                      
    end
       
    % Update the present value of earnings
    pv = pv + prod_adjustment(tm,age,economy,policy) * ps.opt.wstar(year_index,age) * ...
              (1-ps.tax.ya(year_index,age) - ps.tax.wa(year_index,age) )...
              /discount;
          
    % Add in profits
    pv = pv + prod_adjustment(tm,age,economy,policy) * ps.opt.profit(year_index,age) ...
              /discount;
               
          
     % Add in social security income
     pv = pv + prod_adjustment(tm,age,economy,policy)*ps.opt.ben(year_index,age)/discount; 
     
     % Subtract social security taxes
     pv = pv - ps.tax.sst(year_index,age)*prod_adjustment(tm,age,economy,policy)*ps.opt.inc(year_index,age)/discount;
     
     % Add in final debt to PV
     if age==rec.final_age
         
         % Final debt is already adjusted for productivity
         pv = pv + rec.final_debt/discount;
         
     end
     
end

end

