function [ ps,gov,economy,prices] = create_profile(ps,prod,gov,policy,economy,prices,run_schedule,tm)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create Profile                    
% This function takes data already calculated during the last iteration,
% and combines it to create new variables or to reshape old variables in 
% matrices. Most of the data will be in "year by age" matrices;
% I.e., C(I,J) will store data from the consumption of individuals
% in year I who are age J. Thus, to follow the consumption of one
% particular individual you will have to look on a diagonal!

% Arguments:
% tm. This stores information about what years we are calculating this for.  
 

% Initialize Variables
for year=tm.year_min:tm.year_max
       
    % Tax Variables
    economy.ag.wbase(year) = 0;
    economy.ag.wbasesq(year) = 0;
    economy.ag.ybase(year) = 0;
    economy.ag.ybasesq(year) = 0;
    economy.ag.rbase(year) = 0;
    economy.ag.rbasesq(year) = 0;
    
    % Compute Aggregate Interest Income       
    economy.ag.personal_interest_income(year) = 0;
    economy.ag.personal_interest_payments(year) = 0;
    
    % Compute Pure Interest Income
    % (not including annuity interest/payments)
    economy.ag.personal_pure_interest_income(year) = 0;
    economy.ag.personal_pure_interest_payments(year) = 0;
    
    % Annuity Payments / Income
    
    %economy.ag.personal_annuity_income(year) = 0;
    %economy.ag.personal_annuity_payments(year) = 0;
    
    
end

% Loop through types
for type = 1:economy.num_types 
                
    % loop through years
    for year=tm.year_min:tm.year_max
                
        % Loop through ages
        for age=1:ps(type).demog.lifespan
            
            % Now, calculate annuity interest rate
            % Because annuity participation is not 100 percent, the annuity
            % interest rate is a convex combination of the annuity and the non
            % annuity interest rate
            
            % Since ms(t) is the probability surviving from age t to t+1,
            % we need the survival rate ms(t-1)
            
            % Annuity interest rate is (1+r(t))/ms(t-1) - 1
            
            % If age is 1, the annuity rate is simply r
            if age==1
                
                ps(type).prices.r(year,age) = ps(type).util.annuity_participation(year,age)*( (1+prices.r(year))/1 - 1) + (1-ps(type).util.annuity_participation(year,age))*prices.r(year);
                ps(type).prices.r_annuity(year,age) = ps(type).prices.r(year,age);
                
            else % else, annuity rate depends on survival probability of last probability in this year. Thus we need the ym1 index
                 % REVISION couldn't I just used marginal probability of
                 % having survived? 
                 % I'll try it for now
                 
%                 % Need the index of year - 1
%                 if year==1
%                     ym1_index = 1;
%                 else
%                     ym1_index = year - 1;
%                 end
                
                 ps(type).prices.r(year,age) = ps(type).util.annuity_participation(year,age)*( (1+prices.r(year))/ps(type).demog.mhs(year,age) - 1) + (1-ps(type).util.annuity_participation(year,age))*prices.r(year);
                 
                 ps(type).prices.r_annuity(year,age) = ( (1+prices.r(year))/ps(type).demog.mhs(year,age) - 1);
            end
            
            % Set Interest Rates for individual, with financial frictions
            
            if ps(type).opt.a(year,age) > 0
                
                ps(type).prices.r_ff(year,age) = ps(type).prices.r(year,age) - economy.ff(year);
                
            else
                
                ps(type).prices.r_ff(year,age) = ps(type).prices.r(year,age);
                
            end
            
            % Set Interest Rates for Government; this is always the low
            % rate
            
            prices.r_gov(year) = prices.r(year) - economy.ff(year);
            
            % Pre Tax Labor income for individual of age 'age' in year 'year' of type 'type'
            % This is equal to wages X human capital X labor supply
            ps(type).opt.inc(year,age) = prices.wages(year) * ...
                                            ps(type).demog.hc(year,age) * ...
                                            (1-ps(type).opt.l(year,age));
                                        
            % Distribute Profits
            
            if prod.monop_comp==1
                
                % Profits are distributed with respect to labor income :)
                ps(type).opt.profit(year,age) = ( (ps(type).demog.hc(year,age)*(1-ps(type).opt.l(year,age))) ...
                                                / economy.ag.L(year) ) * ...
                                                economy.ag.profit(year);
                                                
            else
                
                ps(type).opt.profit(year,age) = 0;
                                            
            end
                                        
            % Computes human-capital adjusted wages plus shadow wages                     
            ps(type).opt.wstar(year,age) = prices.wages(year) * ...
                                               ps(type).demog.hc(year,age) + ...
                                               ps(type).opt.sw(year,age);
                                           
                                           
            % Computes human-capital adjusted wages                   
            ps(type).opt.wprod(year,age) = prices.wages(year) * ...
                                           ps(type).demog.hc(year,age) ;
                                       
                                       
            % Compute Capital Income
            ps(type).opt.cap_inc(year,age) = ps(type).prices.r_ff(year,age)*ps(type).opt.a(year,age);
                                              
            
            year_born = year - age + 1;
            year_born = max(year_born,1);
            
            % Sometimes debt limit is specified as a percent; otherwise it
            % is specified as a level
            
            % The variable opt.dlim(year,age) specifies the numeric debt
            % limit individuals face. The input can either be as a
            % percentage of income (ps.dlim==1) or as a % of income
            % (ps.dlim==2)
            
            if ps(type).dlim==1 % i.e., specified as a percent
                
                % REVISION -- low priority -- at some point it may be a
                % good idea to specific opt.dlim_pct as a year by age
                % variable, rather a year_born variable. 
                ps(type).opt.dlim(year,age) = ps(type).util.dlim_pct(year_born)*ps(type).opt.inc(year,age);
            
            elseif ps(type).dlim==2 % i.e., specified as a level
                
                % Continue to specify dlim with lenders interest rate
                ps(type).opt.dlim(year,age) = ps(type).util.dlim(year_born,age)/(1+prices.r(year));
            
            end
            
            % New Initialize Population
             
             % In the first year, we set the population of age 1 to be 1
             % All other ages have the same relative relaiton to the age 1
             % population
             
             if year==1
                 
                 % Set Initial populaiton to 1
                 if age==1
                    
                     ps(type).demog.pop(year,age) = 1;
                     
                 % Everyone else, create a constant ratio to the kids at
                 % age 1. We want the relative sizes between ages to
                 % be constant at the growth of the population size;
                 % The population is growign at (# kids / parent)^(1/(age parents have kids))
                 
                 % Finally, we must also adjust for survival
                 else
                     
                     % ps(type).demog.pop(year,age) = ps(type).demog.pop(year,age-1) * ps(type).demog.ms(year,age-1) / ((ps(type).demog.num_kids(year))^(1/(ps(type).demog.age_parent(year)-1)));
                     ps(type).demog.pop(year,age) = ps(type).demog.pop(year,age-1) * ps(type).demog.ms(year,age-1) / (1 + economy.n(1));
                     
                 end
                 
             % Final Steady State Population 
             elseif ((tm.year_min==tm.year_max) && tm.year_min==(policy.nt+2))
                    
                 % Set Initial populaiton to 1
                 if age==1
                    
                     % Note that this is just the rough population! In
                     % the actual transition, the population size will
                     % depend upon actual birthrates, mortality rates, etc
                     ps(type).demog.pop(year,age) = (1+economy.n(year))^(policy.nt+1);
                     
                 else
                     
                     % ps(type).demog.pop(year,age) = ps(type).demog.pop(year,age-1) * ps(type).demog.ms(year,age-1) / (ps(type).demog.num_kids(year))^(1/(ps(type).demog.age_parent(year)-1));
                     
                     ps(type).demog.pop(year,age) = ps(type).demog.pop(year,age-1) * ps(type).demog.ms(year,age-1) / (1+economy.n(policy.nt+2));
                     
                 end
                 
             % Now, not in the first year    
             else
                                  
                 % The number of people who are born in year `year', is
                 % simply the population of their parent generation when
                 % they were age 1, times the number of kids.
                 
                 % In order to find the age of the parent generation, we
                 % simply note that if economic maturity is age 26, when
                 % parents are 26, children are 1. Thus parents are 25
                 % years older than the children
                 
                 if age==1
                     
                     % Assuming the baby boom starts in 1945, this won't
                     % show up in the data until 1970. Thus we have to
                     % simulate all of the additional years from 1945-1970
                     % if wewant to do things right.
                     
                     % Alternatively, we can assume the baby boom happens
                     % in 1945 but nobody "notices" it until 1970. 
                     
                     if run_schedule.demog_shift == 0
                         
                         % I made changes jere July 2016. 
                         % Note that age of parent is constant here, thus using
                         % the index (year) should not be problematic for now. 
                         year_parent_born = max(year - (ps(type).demog.age_parent(year)-1),1);

                         num_kids_year = year_parent_born;
                         
                     elseif run_schedule.demog_shift==1
                         
                         num_kids_year = year;
                         
                     end
                     
                     % Note that the (year-1) here is necessary beacuse we
                     % don't know the population in the current year!!
                     ps(type).demog.pop(year,age) = (1/ps(type).demog.s(year-1,ps(type).demog.age_parent(year)-1)) * ... % Find the parent population size when they are age 1
                                                    ps(type).demog.pop(year-1,(ps(type).demog.age_parent(year)-1)) * ... % population of the parent last year
                                                    ps(type).demog.num_kids(num_kids_year); % Number of kids each parent has
                 
                                                
                 % Everyone else, the population is the population last year times probability of surviving                          
                 else
                     
                     ps(type).demog.pop(year,age) = ps(type).demog.pop(year-1,age-1)*ps(type).demog.ms(year-1,age-1);
                     
                 end
                                  
             end
             
        end % end ages loop
        


    end % end year loop
    
    % Multiply the population by the type share

    ps(type).demog.pop = ps(type).demog.pop*economy.type_share(type);

end % and type loop

% New Series of Loops for the Tax BAse

% Loop through types
for type = 1:economy.num_types 
                
    % loop through years
    for year=tm.year_min:tm.year_max
                
        % Loop through ages
        for age=1:ps(type).demog.lifespan
            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Compute Tax Bases
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % Wage

                economy.ag.wbase(year) = economy.ag.wbase(year) + ps(type).opt.inc(year,age)*ps(type).demog.pop(year,age);
            
            % Wage Squared
            
                economy.ag.wbasesq(year) = economy.ag.wbasesq(year) + ps(type).opt.inc(year,age)^2*ps(type).demog.pop(year,age);
            
            % Capital Base
            
                economy.ag.rbase(year) = economy.ag.rbase(year) + ps(type).opt.cap_inc(year,age)*ps(type).demog.pop(year,age);
        
            % Capital Base SQ
            
                economy.ag.rbasesq(year) = economy.ag.rbasesq(year) + ps(type).opt.cap_inc(year,age)^2*ps(type).demog.pop(year,age);
     
            % Interest Income
                
                % "Interest income" comes from capital only, not on
                % payments to and from annuity firms. This tries to capture
                % this. 
                
                % However, we will still calculate "personal interest
                % income" as an individual's total interest income from the
                % annuity markets
                
                economy.ag.personal_interest_income(year) = economy.ag.personal_interest_income(year) +  (ps(type).opt.cap_inc(year,age)>0)*ps(type).opt.cap_inc(year,age)*ps(type).demog.pop(year,age);
                
                economy.ag.personal_interest_payments(year) = economy.ag.personal_interest_payments(year) +  (ps(type).opt.cap_inc(year,age)<0)*ps(type).opt.cap_inc(year,age)*ps(type).demog.pop(year,age);

            % Pure Interest Income
            % This does not include annuity payments
            % Also, note that it is necessary to include the interest
            % income of people who died already
            
                economy.ag.personal_pure_interest_income(year) = economy.ag.personal_pure_interest_income(year) +  (ps(type).opt.a(year,age)*prices.r(year))*ps(type).demog.pop(year,age)*(1/ps(type).demog.mhs(year,age))*((ps(type).opt.a(year,age)*prices.r(year))>0);
                
                economy.ag.personal_pure_interest_payments(year) = economy.ag.personal_pure_interest_payments(year) +  (ps(type).opt.a(year,age)*prices.r(year))*ps(type).demog.pop(year,age)*(1/ps(type).demog.mhs(year,age))*((ps(type).opt.a(year,age)*prices.r(year))<0);

            % Personal Annuity Income and Payments
            
                % When individuals participate in the annuity markets, they
                % receive annuity income from people who die. 
                
                % The annuity income is the difference between the interest
                % rate they receive and the risk free interest rate
                
                % When they borrow from the annuity market, they pay extra,
                % and these are the annuity payments
            
                % Note we only consider this income if it is greater than
                % zero
                %economy.ag.personal_annuity_income(year) = economy.ag.personal_annuity_income(year)  + ps(type).opt.a(year,age)*(ps(type).prices.r(year,age) - prices.r(year)) * ( (ps(type).opt.a(year,age)*(ps(type).prices.r(year,age) - prices.r(year))) > 0 )*ps(type).demog.pop(year,age)  ;
                
               % economy.ag.personal_annuity_payments(year) = economy.ag.personal_annuity_payments(year)  + ps(type).opt.a(year,age)*(ps(type).prices.r(year,age) - prices.r(year)) * ( (ps(type).opt.a(year,age)*(ps(type).prices.r(year,age) - prices.r(year))) < 0 )*ps(type).demog.pop(year,age)  ;
                
        end % End Age
        
        % Income Tax Base

        economy.ag.ybase(year) = economy.ag.wbase(year) + economy.ag.rbase(year);
        
        economy.ag.ybasesq(year) = economy.ag.wbasesq(year) + economy.ag.rbasesq(year);
        
        
    end % End Year
    
end % End Type


end

