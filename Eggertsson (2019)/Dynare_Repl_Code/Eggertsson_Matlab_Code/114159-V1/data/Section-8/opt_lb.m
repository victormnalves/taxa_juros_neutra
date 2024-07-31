function [ ps ] = opt_lb(ps,gov,prod,policy,prices,economy,run_schedule,year,initial_age,year_min,year_max)

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop through Types     %
%%%%%%%%%%%%%%%%%%%%%%%%%%

for type = 1:economy.num_types

    % Construct present value of earnings, including SS and all taxes
    
    [pv1] = getpv1(ps,economy,prod,gov,policy,prices,year,initial_age);

    % Get variables that describe optimal consumption as a percentage
    % of the present value of income (x1), and zxi which describes
    % consumption as a fraction of consumption at the initial age
    [zxi,x1,nu,H,af] = getx1(ps,policy,prices,year,initial_age);

    % Consumption for initial age
    ps(type).opt.C(year,initial_age) = x1*pv1;

    % Consumption for other ages
    for age=(initial_age + 1):ps(type).demog.lifespan

        year_index = year + age - initial_age;

        if year==1

            year_index = 1;

        end
        
        % I.e., only in final steady state
        if (year_min==year_max) && (year_min==(policy.nt+2))
            
            year_index = policy.nt+2;
            
        end

        % Afterwards, it is a steady state; no need to fill in
        if year_index <= (policy.nt + 2)
            
            ps(type).opt.C(year_index,age) = zxi(age)* ps(type).opt.C(year,initial_age);
        
        end

    end

    % Now get leisure supply
    
    if ps(type).endogenous_labor ==1

        for age = initial_age:ps(type).demog.lifespan

           year_index = year + age - initial_age;

            if year==1

                year_index = 1;

            end

            if (year_min==year_max) && (year_min==(policy.nt+2))

                year_index = policy.nt+2;

            end

            % Afterwards, it is a steady state; no need to fill in
            if year_index <= (policy.nt + 2)

                ps(type).opt.l(year_index,age) = (ps(type).prod_af(age)*ps(type).opt.wstar(year_index,age)*(1-ps(type).tax.ym(year_index,age) - ps(type).tax.wm(year_index,age) )/ps(type).util.alpha(year_index,age))^(-ps(type).util.rho(year_index,age))*ps(type).opt.C(year_index,age); 


            end
          % l_alt(age) = (wstar_alt(age)*(1-ps(type).tax.ym(year_index,type) - ps(type).tax.wm(year_index,type) )/ps(type).util.alpha(year_index))^(-ps(type).util.rho(year_index))*ps(type).opt.C(year_index,age); 

        end
        
    else
        
         ps(type).opt.l = zeros(policy.nt+2,ps(type).demog.lifespan);
    
    end

    % Update assets; remember, a(t,age) are the assets held at the
    % beginning of period 't' at age 'age'
    a_save = zeros(ps(type).demog.lifespan,1);
    for age = initial_age:ps(type).demog.lifespan

       % Get year Born
       year_born = year - initial_age + 1;
       year_born = max(year_born,1);
       
       % this year's index for an individual of age 'age'
       year_index = year + age - initial_age;
             
       % Last year's index for an individual of age 'age'
       ym1_index = year_index - 1;

       agem1_index = age-1;
       
       
        if year==1

            year_index = 1;

            ym1_index = 1;

        end
        
        if (year_min==year_max) && (year_min==(policy.nt+2))
            
            year_index = policy.nt+2;
            
            ym1_index = policy.nt+2;
            
        end
        

        % Afterwards, it is a steady state
        if year_index <= (policy.nt + 2)
            
            if age==1

                ps(type).opt.a(year_index,age) = 0;
                
                % After initial steady state, we add in the transfers from
                % the authority
                
                if year > 1
                    
                    ps(type).opt.a(year_index,age) = ps(type).opt.a(year_index,age) + ...
                                                     ps(type).opt.V(year_index)/(1+prices.r(year_index)*(1- ps(type).tax.ya(year_index,age) - ps(type).tax.ka(year_index,age) )); 
                    
                end
                                
            else

                 ps(type).opt.a(year_index,age) = (1 + prices.r(ym1_index)*( 1 - ps(type).tax.ya(ym1_index,agem1_index) - ps(type).tax.ka(ym1_index,agem1_index))) * ps(type).opt.a(ym1_index,agem1_index) + ...
                                                  (1-ps(type).opt.l(ym1_index,agem1_index)) * ps(type).prod_af(agem1_index) * ps(type).opt.wprod(ym1_index,agem1_index) * ...
                                                  (1-ps(type).tax.ya(ym1_index,agem1_index) - ps(type).tax.wa(ym1_index,agem1_index) - ps(type).tax.sst(ym1_index,agem1_index)) + ...
                                                  ps(type).opt.ben(ym1_index,agem1_index) - ...
                                                  ps(type).opt.C(ym1_index,agem1_index);
                                              
                                              

                  % Giving transfers to the peeps who are born when the
                  % policy changes
                  if initial_age > 1 && age==initial_age
                      
                      ps(type).opt.a(year_index,age) = ps(type).opt.a(year_index,age) + ...
                                                       ps(type).opt.vcur(age)/(1+prices.r(year_index)*(1- ps(type).tax.ya(year_index,age) - ps(type).tax.ka(year_index,age) ));
                                                   
                      
                  end
                         

                    
                 % Add in bequests
                 if age== ps(type).util.a_br(year_born)

                     ps(type).opt.a(year_index,age) = ps(type).opt.a(year_index,age) + ps(type).opt.br(year_index,age);

                 end 
                                              
            end
            
            
            a_save(age) = ps(type).opt.a(year_index,age);

        end

    end

    % Update Shadow Wages

    if ps(type).endogenous_labor ==1
        
        shadow_wages = zeros(ps(type).demog.lifespan,1);

        for age = initial_age:ps(type).demog.lifespan

            year_index = year + age - initial_age;

            if year==1

                year_index = 1;

            end

            if (year_min==year_max) && (year_min==(policy.nt+2))

                year_index = policy.nt+2;

            end

            % Afterwards, it is a steady state
            if year_index <= (policy.nt + 2)

                shadow_wages(age) = (ps(type).util.alpha(year_index,age)/(1-ps(type).tax.ym(year_index,age) - ps(type).tax.wm(year_index,age) ))*(ps(type).opt.C(year_index,age))^(1/ps(type).util.rho(year_index,age)) - ps(type).prod_af(age) *ps(type).opt.wprod(year_index,age);

                if shadow_wages(age) < 0

                    shadow_wages(age) = 0;

                end


                % Update shadow wages with a dampening
                ps(type).opt.sw(year_index,age) =  run_schedule.swdamp*shadow_wages(age) + (1-run_schedule.swdamp)*ps(type).opt.sw(year_index,age);

            end

        end % age for shadow wage

    end

    % Calculate Optimal Bequests given (this is the quantity of bequests
    % per child)
    

    year_born = year - initial_age + 1;
    year_born = max(year_born,1);
    
    bequest_year = year + ps(type).util.a_bg(year_born) - initial_age;
    bequest_year_received = bequest_year + 1;
    
    if year==1 
        
        year_born = 1;
        bequest_year = 1;
        bequest_year_received = 1;
        
    end
    
    if (year_min==year_max) && (year_min==(policy.nt+2))

        year_born = policy.nt+2;
        bequest_year = policy.nt+2;
        
        bequest_year_received = policy.nt+2;

    end

    if bequest_year <= (policy.nt + 2) && bequest_year >=1
    
        ps(type).opt.bg(bequest_year,ps(type).util.a_bg(year_born)) = ps(type).opt.C(bequest_year,ps(type).util.a_bg(year_born)) ...
                                                           * ps(type).util.mu(year_born)^(ps(type).util.gamma(bequest_year,ps(type).util.a_bg(year_born))) ...
                                                           * ps(type).demog.num_kids(year_born)^(-ps(type).util.gamma(bequest_year,ps(type).util.a_bg(year_born))) ...
                                                           * H(ps(type).util.a_bg(year_born))^(-ps(type).util.gamma(year,initial_age));
                                                        
        % Now, set the level of bequests that the kids receive
        % The age of the kids at bequest_year + 1 is 
        % (a_bg + 1) - (age of parents - 1)
        
        child_bequest_age = (ps(type).util.a_bg(year_born)+1) - (ps(type).demog.age_parent(year_born) - 1);
        
        if bequest_year_received <= (policy.nt + 2)
            
            ps(type).opt.br(bequest_year_received,child_bequest_age) = ps(type).opt.bg(bequest_year,ps(type).util.a_bg(year_born)) ; 
                                                                       % / ps(type).demog.pop(bequest_year_received,child_bequest_age);
        
        end
    end
    
    
    % Now, update shadow prices!
    % Only if there is a debt limit
    if ps(type).dlim > 0
        
        spnew = zeros(ps(type).demog.lifespan,1); % Will store the updated shadow prices

        for age=initial_age:(ps(type).demog.lifespan-1) % No need to update the final shadow prices, for it is zero; you can't die in debt

            year_index = year + age - initial_age;
            yp1_index = year_index + 1;

            if year==1

                year_index = 1;
                yp1_index = 1;

            end

            if (year_min==year_max) && (year_min==(policy.nt+2))

                year_index = policy.nt+2;
                yp1_index = policy.nt+2;

            end

            % When updating in the final year, afterwards we are in a steady
            % state
            if year_index==policy.nt+2

                yp1_index = policy.nt+2;

            end


            % Afterwards, it is a steady state
            if year_index <= (policy.nt + 2)

                % Level of consumption this period which leaves assets next
                % period at the borrowing limit

                C_res = ps(type).opt.C(year_index,age) + ps(type).opt.a(yp1_index,age+1) + ps(type).opt.dlim(yp1_index,age + 1);

                % This comes from the formula I derived

                spnew(age) = af(age+1)*ps(type).opt.C(yp1_index,age+1)^(1/ps(type).util.gamma(year_index,age)) ...
                             * (1+ps(type).util.delta(year_index,age)) ...
                             / ( C_res^(1/ps(type).util.gamma(year_index,age)) ...
                             * (1+ps(type).prices.r(yp1_index,age+1)*(1-ps(type).tax.ym(yp1_index,age+1) - ps(type).tax.km(yp1_index,age+1) ) ) ) ...
                             * (nu(age)/nu(age+1))^(1/ps(type).util.gamma(year_index,age))...
                             - 1;

                % The whole problem with this is that we need to look forward
                % to get where the old shadow prices are stored
                % Potentially in the future I coudl change this so in every
                % row, we store all of the shadow prices for people born in
                % year 'rownum'
                for newage= (age+1):(ps(type).demog.lifespan)

                    alt_year_index = year_index + newage - age;

                    if year==1

                        alt_year_index = 1;

                    end

                    if (year_min==year_max) && (year_min==(policy.nt+2))

                        alt_year_index = policy.nt+2;

                    end

                    % When updating in the final year, afterwards we are in a steady
                    % state
                    if alt_year_index>policy.nt+2

                        alt_year_index = policy.nt+2;

                    end

                    spnew(age) = spnew(age) - r_mult_fin(ps(type),policy,initial_age,year,newage)*ps(type).opt.sp(alt_year_index,newage);

                end

                spnew(age) = spnew(age) / r_mult_fin(ps(type),policy,initial_age,year,age) ;

                if spnew(age) < 0

                    spnew(age) = 0;

                end

                ps(type).opt.sp(year_index,age) = .05*spnew(age) + .95*ps(type).opt.sp(year_index,age);

            end

        end

    end
    
    % Make updates to Save Year
    
    % Get year Born
    year_born = year - initial_age + 1;
    year_born = max(year_born,1);
    
    savnew =  find(a_save > 0.001 ,1) ;
    
    if isempty(savnew)
        
        savnew = ps(type).demog.lifespan + 1;
        
    end
    
    ps.savyear(year_born) = savnew;
    
    % Make adjustments because of productivity!!!!!
    % We have already made adjustments to the PVI. 
    % But now, we need to realize that people born in previous generations
    % had different productivities, and thus we need to adjust their
    % optimal quantities down by a factor
    
    
    if year==1 || ((year_min==year_max) && (year_min==(policy.nt+2)))
        ps(type).opt.C_iss = ps(type).opt.C(1,:);
        ps(type).opt.a_iss = ps(type).opt.a(1,:);
        
        for age=initial_age:ps(type).demog.lifespan
            
            %prod_af = 1/( (1+economy.lprod_growth_iss)^(1/(1-prod.epsilon(1))) )^(age - 1);

            ps(type).opt.C(year,age) = ps(type).opt.C(year,age)*(1/ps(type).prod_af(age));
            
            ps(type).opt.a(year,age) = ps(type).opt.a(year,age)*(1/ps(type).prod_af(age));
            
        end
        
        % Adjust bequests received for productivity
        % It is always your parents who give you bequests; so the
        % productivity difference is the age of the parents
        % Note that it may or may not be necessary to adjust bequests
        % given, but I am doing it anyway for now.

        ps(type).opt.br(year,ps(type).util.a_br(year)) = ps(type).opt.br(year,ps(type).util.a_br(year))*(1/ps(type).prod_af(ps(type).util.a_bg(year))) * (1/ps(type).prod_af(2));
        
        ps(type).opt.bg(year,ps(type).util.a_bg(year)) = ps(type).opt.bg(year,ps(type).util.a_bg(year))*(1/ps(type).prod_af(ps(type).util.a_bg(year)));

    end
    
    
end % Type

end

