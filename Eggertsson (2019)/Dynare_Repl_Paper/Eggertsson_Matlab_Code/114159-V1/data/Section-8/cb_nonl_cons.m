function [c,ceq,a] = cb_nonl_cons(CB,ps,prices,economy,policy,tm)

% This function displays whether or not an individual satisfies the budget
% constraint, given a choice of consumption and income. 

% It would be a simple budget constraint, but there is the wrinkle of the
% financial friction 

% Thus this function proceeds in two steps:
    % (1) Given Y and C, calculate assets; if assets are negative, the
    % interest rate for that year is set to be the borrower's rate; else it
    % is set to be the lender's rate
    
    % (2) Given the interest rates calculated above, calculate the pdv of
    % income and the pdv of consumption

    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Step 1: Determine whether there is a borrowers or lender's interest
    % rate
    % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    r = zeros(ps.demog.lifespan,1); % Vector of interest rates faced
    
    a = zeros(ps.demog.lifespan,1); % Assets of individual
    
    % In order to do so, need to calculate assets
    
    for age = tm.initial_age:ps.demog.lifespan

        % update year index
        tm.age = age;
        [year_index,year_born,ym1_index] = create_index(tm,policy);
        
        agem1_index = age-1;
           
        % If age=1, special case
        if age==1

            a(age) = 0;

            % After initial steady state, we add in the transfers from
            % the authority

            if tm.initial_year > 1

                if ps.opt.V(year_index) > 0 % Lending rate
                    
                    r(1) = ps.prices.r(year_index,age) - economy.ff(year_index);
                    
                else
                    
                    r(1) = ps.prices.r(year_index,age);
                    
                end
                
                a(age) = a(age) + ps.opt.V(year_index)/(1+r(age)*(1- ps.tax.ya(year_index,age) - ps.tax.ka(year_index,age) )); 


            end


        else

              % If we are at initial age and initial_age > 1, assets
              % come from decisions made last period.
              if tm.initial_age > 1 && age==tm.initial_age
                 
                 % Need to know if you were a borrower or a lender last
                 % period :)
                 
                 if  ps.opt.a(ym1_index,agem1_index) > 0 % Lending Rate
                     
                     r_prev = ps.prices.r(ym1_index,agem1_index) - economy.ff(ym1_index);
                     
                 else
                     
                     r_prev = ps.prices.r(ym1_index,agem1_index);
                     
                 end
                 
                 a(age) = (1 + r_prev*( 1 - ps.tax.ya(ym1_index,agem1_index) - ps.tax.ka(ym1_index,agem1_index))) * ps.opt.a(ym1_index,agem1_index) + ...
                                                  prod_adjustment(tm,agem1_index,economy,policy) * ps.opt.wprod(ym1_index,agem1_index) * ... % I'm pretty sure this prod adjustment does nothing
                                                  (1-ps.tax.ya(ym1_index,agem1_index) - ps.tax.wa(ym1_index,agem1_index) - ps.tax.sst(ym1_index,agem1_index)) + ...
                                                  prod_adjustment(tm,agem1_index,economy,policy) * ps.opt.profit(ym1_index,agem1_index)+ ...
                                                  ps.opt.ben(ym1_index,agem1_index) - ...
                                                  ps.opt.C(ym1_index,agem1_index);

              else

                 a(age) = (1 + r(agem1_index)*( 1 - ps.tax.ya(ym1_index,agem1_index) - ps.tax.ka(ym1_index,agem1_index))) * a(age-1) + ...
                                                  prod_adjustment(tm,agem1_index,economy,policy) * ps.opt.wprod(ym1_index,agem1_index) * ...
                                                  (1-ps.tax.ya(ym1_index,agem1_index) - ps.tax.wa(ym1_index,agem1_index) - ps.tax.sst(ym1_index,agem1_index)) + ...
                                                  prod_adjustment(tm,agem1_index,economy,policy) * ps.opt.profit(ym1_index,agem1_index) + ...
                                                  ps.opt.ben(ym1_index,agem1_index) - ...
                                                  CB(age-1);
              end


             % Add in bequests
             if age== ps.util.a_br(year_born)
                 
                 % Note that we do adjust for productivity here. What is stored in
                 % br is simply the bequests received by people who are currently
                 % that age. But as productivity increases, so do bequests
             
                 a(age) = a(age) + prod_adjustment(tm,age,economy,policy)*ps.opt.br(year_index,age);

             end 
             
              % Giving transfers to the peeps who are born when the
              % policy changes
              if tm.initial_age > 1 && age==tm.initial_age

                  % We need to know if, given the transfers from the
                  % authority, the individual has positive or negative
                  % assets
                  if a(age) + ps.opt.vcur(age) > 0
                      
                      r_current = ps.prices.r(year_index,age) - economy.ff(year_index);
                      
                  else
                      
                      r_current = ps.prices.r(year_index,age);
                      
                  end
                  
                  a(age) = a(age) + ...
                           ps.opt.vcur(age)/(1+r_current*(1- ps.tax.ya(year_index,age) - ps.tax.ka(year_index,age) ));


              end

        end
            
            
        % Update interest rate
        if a(age) > 0 % if you are a lender, low interest rate
            
            r(age) = ps.prices.r(year_index,age) - economy.ff(year_index);
            
        else  % you are a borrower, high interest rate
            
            r(age) = ps.prices.r(year_index,age);
            
        end
            

            
    end
    
    % Now, construct "c" 
    
    c = zeros(ps.demog.lifespan,1);
    for age=1:ps.demog.lifespan
        
        tm.age = age;
        [year_index,year_born,ym1_index] = create_index(tm,policy);
        
        c(age) = -a(age) - prod_adjustment(tm,age,economy,policy)*ps.opt.dlim(year_index,age);
        
    end
    

% Now, we simply have the budget constraint; pdvc should equal pdvy

    pv = 0;
    pdvcb = 0;
    
    discount = 1;

    % Loop thorugh ages, beginning at beginning_age, and ending at final_age
    for age=tm.initial_age:ps.demog.lifespan

        % update year index
        tm.age = age;
        [year_index,year_born] = create_index(tm,policy);

         % Add in bequests -- note that we use the discount from last year, as
         % in AK, since the bequests are received in the beginning of the
         % period


         if age== ps.util.a_br(year_born) 

             % Note that we do adjust for productivity here. What is stored in
             % br is simply the bequests received by people who are currently
             % that age. But as productivity increases, so do bequests
             
             pv = pv + prod_adjustment(tm,age,economy,policy)*ps.opt.br(year_index,age)/discount;

         end

        % We don't discount the initial year
        if age > tm.initial_age

            % Update Discount; we use average tax rates for the discount

            discount = discount*(1+r(age) * ...
                      (1-ps.tax.ka(year_index,age)-ps.tax.ya(year_index,age)));


        end

        % Update the present value of earnings
        pv = pv + prod_adjustment(tm,age,economy,policy) * ps.opt.wstar(year_index,age) * ...
                  (1-ps.tax.ya(year_index,age) - ps.tax.wa(year_index,age) )...
                  /discount;
              
        % Add Profits
        pv = pv + prod_adjustment(tm,age,economy,policy) * ps.opt.profit(year_index,age) ...
                  /discount;

         % Add in social security income
         pv = pv + ps.opt.ben(year_index,age)/discount; 

         % Subtract social security taxes
         pv = pv - ps.tax.sst(year_index,age)*prod_adjustment(tm,age,economy,policy)*ps.opt.inc(year_index,age)/discount;


         % Now, calculate pdvcb
         
            % Consumption
            
            pdvcb = pdvcb + CB(age)/discount;
            
            % Bequest Given
            
            if age== ps.util.a_bg(year_born) 

                pdvcb = pdvcb + ps.demog.num_kids(year_born)*CB(ps.demog.lifespan + 1)/discount;

            end
            
    end

   
    % Finally, the budget constraint
    
    ceq = pv - pdvcb;
    
%     temp = 0;
%     for age=1:56
%         
%         temp = temp + CB(age)/(1+prices.r(1))^(age-1);
%         
%     end
%     
%     temp = temp + CB(57)/(1+prices.r(1))^(55);
%     
    
end

