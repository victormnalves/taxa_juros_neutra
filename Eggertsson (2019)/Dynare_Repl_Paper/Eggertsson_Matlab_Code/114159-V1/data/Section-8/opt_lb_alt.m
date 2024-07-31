function [ C, a, neg_flag] = opt_lb_alt(ps,gov,prod,policy,prices,economy,run_schedule,tm,rec)

% Calculates optimal consumption from rec.beginning_age to rec.final_age
% assuming no borrowing limit
% Does this through formulas mainly taken from Auerbach & Kotlikoff. We
% first calculate the present value of income, then calculte consumption in
% period 1 as a function of this value. Afterwards, we calculate optimal
% assets. We don't have any long term storage of this information, we
% simply pass them back as vectors C and a. 

% Finally, if the optimal consumption bundle violates the borrowing limit,
% we will send back neg_flag = 1

    % In this program we make a distinction between initial age 
    % and beginning age

    % initial age is the age an individual starts to to optimize from after
    % there is a change in the policy, and the age at which the individual
    % receives his or her transfers from the authority

    % beginning age is the age we optimize from in the recursive loop, starting
    % out with initial_debt


    % Construct present value of earnings, including SS and all taxes
    [pv1] = getpv1_alt(ps,economy,prod,gov,policy,prices,tm,rec);
    

    
    % Get variables that describe optimal consumption as a percentage
    % of the present value of income (x1), and zxi which describes
    % consumption at each age relative to initial consumption
    [zxi,x1] = getx1_alt(ps,policy,prices,economy,tm,rec);
               
    % Consumption for initial age
    C = zeros(ps.demog.lifespan,1);
    
    C(rec.beginning_age) = x1*pv1;

    % Consumption for other ages
    for age=(rec.beginning_age + 1):rec.final_age

        C(age) = zxi(age) * C(rec.beginning_age);

    end
    
%     % Fill In Consumption for age before initial_age
%     % This is necessary to calculate assets at initial age
%     if tm.initial_age > 1 && rec.beginning_age==tm.initial_age
%         
%         C(tm.initial_age - 1) = ps.opt.C(tm.initial_year - 1,tm.initial_age - 1);
%         
%     end



    % Now get leisure supply
    l = zeros(ps.demog.lifespan,1);
    if ps.endogenous_labor ==1

        for age = rec.beginning_age:rec.final_age

            % update year index
            tm.age = age;
            [year_index,year_born] = create_index(tm,policy);

            l(age) = (prod_adjustment(tm,age,economy,policy)*ps.opt.wstar(year_index,age)*(1-ps.tax.ym(year_index,age) - ps.tax.wm(year_index,age) )/ps.util.alpha(year_born))^(-ps.util.rho(year_born))*C(age); 


        end
        
    else
        
         l = zeros(ps.demog.lifespan,1);
    
    end
    
    % Update assets; remember, a(t,age) are the assets held at the
    % beginning of period 't' at age 'age'
    
    a = zeros(ps.demog.lifespan,1);
    
    neg_flag = 0; % Flag whether there is too much borrowing
    
    for age = rec.beginning_age:rec.final_age % Dec 2017 note that this is different than the new code in the profits folder, opt_lb_alt....not sure what's going on. In the new code, it's + 1

        % update year index
        tm.age = age;
        [year_index,year_born,ym1_index] = create_index(tm,policy);
        
        agem1_index = age-1;
           
        % Afterwards, it is a steady state
            
            if age==1

                a(age) = 0;
                
                % After initial steady state, we add in the transfers from
                % the authority
                
                if tm.initial_year > 1
                    
                    a(age) = a(age) + ps.opt.V(year_index)/(1+ps.prices.r(year_index,age)*(1- ps.tax.ya(year_index,age) - ps.tax.ka(year_index,age) )); 
                                                     
                    
                end
                                

            else

                  % If beginning age > initial_age, we neglect whatever assets you may have before
                  % and say you simply come into the world with
                  % initial_debt. Thus we have a(beginning_age) =
                  % -initial_debt
                  if age==rec.beginning_age && rec.beginning_age > tm.initial_age
                      
                      a(age) = a(age) - rec.initial_debt;
                      
                  % If we are at initial age and initial_age > 1, assets
                  % come from decisions made last period.
                  % REVISION -- could we simply take assets from period 1,
                  % updated by productivity growth, as we do in getpv1_alt?
                  % Since age==initial_age always occurs when year==2...
                  
                  elseif tm.initial_age > 1 && rec.beginning_age == tm.initial_age && age==tm.initial_age
                      
                     a(age) = (1 + ps.prices.r(ym1_index,agem1_index)*( 1 - ps.tax.ya(ym1_index,agem1_index) - ps.tax.ka(ym1_index,agem1_index))) * ps.opt.a(ym1_index,agem1_index) + ...
                                                      (1-ps.opt.l(ym1_index,agem1_index)) * prod_adjustment(tm,agem1_index,economy,policy) * ps.opt.wprod(ym1_index,agem1_index) * ... % I'mt not convinced this productivity adjustment does anything -- it may be superfluous. Remain for now. 
                                                      (1-ps.tax.ya(ym1_index,agem1_index) - ps.tax.wa(ym1_index,agem1_index) - ps.tax.sst(ym1_index,agem1_index)) + ...
                                                      + prod_adjustment(tm,agem1_index,economy,policy) * ps.opt.profit(ym1_index,agem1_index) + ...
                                                      ps.opt.ben(ym1_index,agem1_index) - ... 
                                                      ps.opt.C(ym1_index,agem1_index);
                                            
                  else
                      
                     % prod_adjustment: for ISS, need to calculate assets
                     % at age 5, which depends on wage at age 4. Need to
                     % update what the wage will be. 
                     a(age) = (1 + ps.prices.r(ym1_index,agem1_index)*( 1 - ps.tax.ya(ym1_index,agem1_index) - ps.tax.ka(ym1_index,agem1_index))) * a(age-1) + ...
                                                      (1-l(age-1)) * prod_adjustment(tm,agem1_index,economy,policy) * ps.opt.wprod(ym1_index,agem1_index) * ...
                                                      (1-ps.tax.ya(ym1_index,agem1_index) - ps.tax.wa(ym1_index,agem1_index) - ps.tax.sst(ym1_index,agem1_index)) + ...
                                                      prod_adjustment(tm,agem1_index,economy,policy) * ps.opt.profit(ym1_index,agem1_index) + ...
                                                      prod_adjustment(tm,agem1_index,economy,policy)*ps.opt.ben(ym1_index,agem1_index) - ... % Since there is a possiblity this can be someone from year 1, we do adjust for productivity
                                                      C(age-1);
                  end
                  
                  % Giving transfers to the peeps who are born when the
                  % policy changes
                  if tm.initial_age > 1 && age==tm.initial_age
                      
                      a(age) = a(age) + ...
                               ps.opt.vcur(age)/(1+ps.prices.r(year_index,age)*(1- ps.tax.ya(year_index,age) - ps.tax.ka(year_index,age) ));
                                                   
                      
                  end
                         

                    
                 % Add in bequests
                 if age== ps.util.a_br(year_born)

                     % Note that we do adjust for productivity here. What is stored in
                     % br is simply the bequests received by people who are currently
                     % that age. But as productivity increases, so do bequests

                     a(age) = a(age) + prod_adjustment(tm,age,economy,policy)*ps.opt.br(year_index,age);

                 end 
                                              
            end
            
            % Implement Flag Procedure
            if ps.dlim > 0
                if (a(age) + prod_adjustment(tm,age,economy,policy)*ps.opt.dlim(year_index,age)) < 0 && (age > rec.beginning_age)

                    neg_flag = 1;

                end
            end
            
    end
    
    if rec.beginning_age==rec.final_age 
        
        neg_flag = 0;
        
    end
    
    if ps.dlim==0
        
        neg_flag = 0;
        
    end
    
    % Finally, restrict the vectors sent
    C = C(rec.beginning_age:rec.final_age);
    a = a(rec.beginning_age:rec.final_age);

end

