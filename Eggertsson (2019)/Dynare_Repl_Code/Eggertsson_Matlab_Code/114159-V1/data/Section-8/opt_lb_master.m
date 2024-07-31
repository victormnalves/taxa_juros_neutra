function [ ps ] = opt_lb_master(ps,gov,prod,policy,prices,economy,run_schedule,tm,rec,C,a)

% This program mainly stores and compiles the information computed about
% the optimal consumption bundles already computed in the vectors C and a


    % The only things we really need are the H from this    
    rec.beginning_age = tm.initial_age;
    rec.final_age = ps.demog.lifespan;
    [zxi,x1,nu,H] = getx1_alt(ps,policy,prices,economy,tm,rec);

    % Put Consumption away
    % First, extend it so that the indexes align
    if tm.initial_age > 1
        
       C = [zeros(tm.initial_age-1,1);C];
       
    end
    
    % [C,ps.opt.C(1,:)']

    
    for age=(tm.initial_age):ps.demog.lifespan

        % update year index
        tm.age = age;
        [year_index,~,~,~,year_alt_index] = create_index(tm,policy);
        % Afterwards, it is a steady state; no need to fill in
        if year_alt_index <= (policy.nt + 2)
            
            ps.opt.C(year_index,age) = C(age);
        
        end

    end

    % Now get leisure supply
    
    if ps.endogenous_labor ==1

        for age = tm.initial_age:ps.demog.lifespan

            % update year index
            tm.age = age;
            [year_index,year_born,~,~,year_alt_index]= create_index(tm,policy);

            % Afterwards, it is a steady state; no need to fill in
            if year_alt_index <= (policy.nt + 2)

                ps.opt.l(year_index,age) = (prod_adjustment(tm,age,economy,policy)*ps.opt.wstar(year_index,age)*(1-ps.tax.ym(year_index,age) - ps.tax.wm(year_index,age) )/ps.util.alpha(year_born))^(-ps.util.rho(year_born))*ps.opt.C(year_index,age); 

            end

        end
        
    else
        
         ps.opt.l = zeros(policy.nt+2,ps.demog.lifespan);
    
    end

   % Put Assets away
    
    % First, extend it so that the indexes align
    if tm.initial_age > 1
        
       a = [zeros(tm.initial_age-1,1);a];
       
    end
    
    for age=(tm.initial_age):ps.demog.lifespan

        % update year index
        tm.age = age;
        [year_index,~,~,~,year_alt_index] = create_index(tm,policy);
        
        % Afterwards, it is a steady state; no need to fill in
        if year_alt_index <= (policy.nt + 2)
            
            ps.opt.a(year_index,age) = a(age);
            
        
        end

    end

    % Update Shadow Wages

    if ps.endogenous_labor ==1
        
        shadow_wages = zeros(ps.demog.lifespan,1);

        for age = tm.initial_age:ps.demog.lifespan

            % update year index
            tm.age = age;
            [year_index,year_born,~,~,year_alt_index]= create_index(tm,policy);

            % Afterwards, it is a steady state
            if year_alt_index <= (policy.nt + 2)

                shadow_wages(age) = (ps.util.alpha(year_born)/(1-ps.tax.ym(year_index,age) - ps.tax.wm(year_index,age) ))*(ps.opt.C(year_index,age))^(1/ps.util.rho(year_born)) - prod_adjustment(tm,age,economy,policy)*ps.opt.wprod(year_index,age);

                if shadow_wages(age) < 0

                    shadow_wages(age) = 0;

                end

                % Update shadow wages with a dampening
                ps.opt.sw(year_index,age) =  run_schedule.swdamp*shadow_wages(age) + (1-run_schedule.swdamp)*ps.opt.sw(year_index,age);

            end

        end % age for shadow wage

    end

    % Calculate Optimal Bequests given (this is the quantity of bequests
    % per child)
    

    [~,year_born] = create_index(tm,policy);
    
    bequest_year = tm.initial_year + ps.util.a_bg(year_born) - tm.initial_age;
    bequest_year_received = bequest_year + 1;
    
    if tm.initial_year==1 
        
        year_born = 1;
        bequest_year = 1;
        bequest_year_received = 1;
        
    end
    
    if (tm.year_min==tm.year_max) && (tm.year_min==(policy.nt+2))

        year_born = policy.nt+2;
        bequest_year = policy.nt+2;
        
        bequest_year_received = policy.nt+2;

    end

    if bequest_year <= (policy.nt + 2) && bequest_year >=1
    
        ps.opt.bgo(bequest_year,ps.util.a_bg(year_born)) = ps.opt.C(bequest_year,ps.util.a_bg(year_born)) ...
                                                           * ps.util.mu(year_born)^(ps.util.gamma(year_born)) ...
                                                           * ps.demog.num_kids(year_born)^(-ps.util.gamma(year_born)) ...
                                                           * H(ps.util.a_bg(year_born))^(-ps.util.gamma(year_born));
                                  
        % Total bequests given = bequest per kid * Num Kids * population                                               
        ps.opt.bg_adj(bequest_year,ps.util.a_bg(year_born)) = ps.opt.bgo(bequest_year,ps.util.a_bg(year_born))*ps.demog.pop(bequest_year,ps.util.a_bg(year_born))*ps.demog.num_kids(year_born);                                            
        

        % Now, set the level of bequests that the kids receive
        % The age of the kids at bequest_year + 1 is 
        % (a_bg + 1) - (age of parents - 1)
        
        child_bequest_age = (ps.util.a_bg(year_born)+1) - (ps.demog.age_parent(year_born) - 1);
        
        if bequest_year_received <= (policy.nt + 2)
            
            ps.opt.bro(bequest_year_received,child_bequest_age) = ps.opt.bgo(bequest_year,ps.util.a_bg(year_born)) ; 
                                                                       % / ps.demog.pop(bequest_year_received,child_bequest_age);
             
            % Total bequests received                                                           
            ps.opt.br_adj(bequest_year_received,child_bequest_age) = ps.opt.bg_adj(bequest_year,ps.util.a_bg(year_born));
            
            % ps.opt.brt(bequest_year_received) = ps.opt.brt(bequest_year_received) + ps.opt.br_adj(bequest_year_received,child_bequest_age);
            
            
        end
    end
    

    
    
    % If this is year 2, we need to get the bequest received from year 1 :)
    % Since this is the first year of the transition
    if tm.initial_year==2 && tm.initial_age==1;
        
        child_bequest_age = (ps.util.a_bg(1)+1) - (ps.demog.age_parent(1) - 1);
         
        ps.opt.bro(2,child_bequest_age) = ps.opt.bgo(1,ps.util.a_bg(1));
        
        ps.opt.br_adj(2,child_bequest_age) = ps.opt.bg_adj(1,ps.util.a_bg(1));
        
        % ps.opt.brt(2) = ps.opt.brt(2) + ps.opt.br_adj(2,child_bequest_age);

    end
    
    
   
    % Make updates to Save Year
    
%     ps.savyear =  find(a > 0.001 ,1) ;
%     
%     if isempty(ps.savyear)
%         
%         ps.savyear = ps.demog.lifespan + 1;
%         
%     end
    
    % ps.savyear(year_born) = savnew;
    
    % Make adjustments because of productivity!!!!!
    % We have already made adjustments to the PVI. 
    % But now, we need to realize that people born in previous generations
    % had different productivities, and thus we need to adjust their
    % optimal quantities down by a factor
    
    
    if tm.initial_year==1 || ((tm.year_min==tm.year_max) && (tm.year_min==(policy.nt+2)))
        
        ps.opt.C_iss = ps.opt.C(1,:);
        ps.opt.a_iss = ps.opt.a(1,:);
        
        ps.opt.u_new(1) = -crrautil([ps.opt.C_iss';.00001],ps,tm,policy);
        
        
        for age=tm.initial_age:ps.demog.lifespan
            
            %prod_af = 1/( (1+economy.lprod_growth_iss)^(1/(1-prod.epsilon(1))) )^(age - 1);
            
            ps.opt.C(tm.initial_year,age) = ps.opt.C(tm.initial_year,age)*(1/prod_adjustment(tm,age,economy,policy));
            
            ps.opt.a(tm.initial_year,age) = ps.opt.a(tm.initial_year,age)*(1/prod_adjustment(tm,age,economy,policy));
            
        end
        
        % Adjust bequests received for productivity
        % It is always your parents who give you bequests; so the
        % productivity difference is the age of the parents
        % Note that it may or may not be necessary to adjust bequests
        % given, but I am doing it anyway for now.

%         if tm.initial_year==1
%             
%             temp_adjustment = (1+economy.lprod_growth_iss);
%             
%         elseif tm.initial_year==(policy.nt + 2)
%             
%             temp_adjustment = (1+economy.lprod_growth_fss);
%         end
        
        % Since the person who gave hte bequest to the person receiving it
        % at year 1 is age 56, we need to compare his optimal bequest
        % choice to hte person born at age 1.
        
        ps.opt.bro(tm.initial_year,ps.util.a_br(tm.initial_year)) = ps.opt.bro(tm.initial_year,ps.util.a_br(tm.initial_year))*(1/prod_adjustment(tm,ps.util.a_bg(tm.initial_year) + 1,economy,policy)) ;
        
        % Adjust bequests received for population
        % It was given by the people in year 0

        pop_adj = (1 + economy.n(1));
            
        
        ps.opt.br_adj(tm.initial_year,ps.util.a_br(tm.initial_year)) = ps.opt.br_adj(tm.initial_year,ps.util.a_br(tm.initial_year))*(1/prod_adjustment(tm,ps.util.a_bg(tm.initial_year) + 1,economy,policy))*(1/pop_adj) ;
        
        % ps.opt.brt(tm.initial_year,ps.util.a_br(tm.initial_year)) = ps.opt.brt(tm.initial_year,ps.util.a_br(tm.initial_year))*(1/prod_adjustment(tm,ps.util.a_bg(tm.initial_year) + 1,economy,policy))*(1/pop_adj(tm,economy,policy)) ;
        
        % The person giving the bequest is only 55 -- thus there is less of
        % a discount
        ps.opt.bgo(tm.initial_year,ps.util.a_bg(tm.initial_year)) = ps.opt.bgo(tm.initial_year,ps.util.a_bg(tm.initial_year))*(1/prod_adjustment(tm,ps.util.a_bg(tm.initial_year),economy,policy));

        
        ps.opt.bg_adj(tm.initial_year,ps.util.a_bg(tm.initial_year)) = ps.opt.bg_adj(tm.initial_year,ps.util.a_bg(tm.initial_year))*(1/prod_adjustment(tm,ps.util.a_bg(tm.initial_year),economy,policy));

        % ps.opt.bgt(tm.initial_year,ps.util.a_bg(tm.initial_year)) = ps.opt.bgt(tm.initial_year,ps.util.a_bg(tm.initial_year))*(1/prod_adjustment(tm,ps.util.a_bg(tm.initial_year),economy,policy));
    
    end
    
    % %%%%%%%%%%%%%%%%%%%%
    % Calculate Savings  %
    % %%%%%%%%%%%%%%%%%%%%
    
    % This is a little bit complicated; we can't simply do income minus
    % consumption, because there is shit like bequests 
    
    % Instead, we calculate savings as assets next period minus assets this
    % period
    
    % We will also have to account for productivity growth in the initial
    % steady states
    
    for age = tm.initial_age:(ps.demog.lifespan)

    % update year index
    tm.age = age;
    [year_index,year_born,ym1_index,yp1_index,year_alt_index] = create_index(tm,policy);

    % Afterwards, it is a steady state
    if year_alt_index <= (policy.nt + 2)
        
        % If it is the final age, set savings to negative of assets
        if age==ps.demog.lifespan
            
            ps.opt.savings(year_index,age) = -ps.opt.a(year_index,age);
            
        else

            % Adjust for productivity in the initial steady state
            if tm.year_min==tm.year_max && tm.year_min==1

                ps.opt.savings(year_index,age) = ps.opt.a(yp1_index,age+1)*(1+economy.ag.AL_growth_iss)*(1+economy.ag.A_growth_adj_iss) - ps.opt.a(year_index,age);

            elseif year_alt_index == (policy.nt + 2) % account for productivity growth in the final steady state

                ps.opt.savings(year_index,age) = ps.opt.a(yp1_index,age+1)*(1+economy.ag.AL_growth_fss)*(1+economy.ag.A_growth_adj_fss) - ps.opt.a(year_index,age);

            else % otherwise, we do not adjust for productivity, because it is already adjusted!

                ps.opt.savings(year_index,age) = ps.opt.a(yp1_index,age+1) - ps.opt.a(year_index,age);

            end

        end
        
    end
    
    end
    
    % Final year savings are zero :)
    
    
    
%     pvalt = 0;
%     for age=(tm.initial_age):ps.demog.lifespan
%         
%         tm.age = age;
%         [year_index] = create_index(tm,policy);
%         
%         pvalt = pvalt + C(age)/(1+prices.r(1))^(age-tm.initial_age);
%         
%     end
%     
%     pvalt = pvalt + ps.opt.bgo(bequest_year,ps.util.a_bg(year_born))/(1+prices.r(1))^(ps.util.a_bg(year_born)-tm.initial_age);
%     
%     
    % Now, add to the bequests from people that die! 
%     
%     for age=(tm.initial_age):ps.demog.lifespan
% 
%         % update year index
%         tm.age = age;
%         [year_index,year_born,ym1_index] = create_index(tm,policy);
%         
%         prod_adj = 1;
%         pop_adj = 1;
%         
%         if year_index==1
%             
%             pop_adj = (1 + economy.n(1));
%             
%         elseif year_index==(policy.nt+2)
%             
%             prod_adj = (1+ economy.lprod_growth_fss);
%             
%         end
%         
%         % Afterwards, it is a steady state; no need to fill in
%         if year_index <= (policy.nt + 2)
%             
%             if age > 1
%                 
%                 ps.opt.brt(year_index) = ps.opt.brt(year_index) + ps.opt.a(year_index,age)*ps.demog.pop(ym1_index,age-1)*(1-ps.demog.ms(ym1_index,age-1))*(1/pop_adj);
%                 
%                 ps.opt.bgt(year_index) = ps.opt.bgt(year_index) + ps.opt.a(yp1_index,age)*ps.demog.pop(year_index,age-1)*(1-ps.demog.ms(year_index,age-1))*prod_adj;
%             
%             end
%         
%         end
% 
%     end
    


end

