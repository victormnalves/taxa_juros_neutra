function [ C, l] = opt_lb_final_calc(ps,gov,prod,policy,prices,economy,run_schedule,year,initial_age,year_min,year_max,C)

% Calculates the Final Information we need

  
    % Now get leisure supply
    
    l = zeros(ps.demog.lifespan,1);
    if ps.endogenous_labor ==1

        for age = initial_age:ps(type).demog.lifespan

           year_index = year + age - initial_age;

            if year==1

                year_index = 1;

            end

            if (year_min==year_max) && (year_min==(policy.nt+2))

                year_index = policy.nt+2;

            end

                l(age) = (ps(type).prod_af(age)*ps(type).opt.wstar(year_index,age)*(1-ps(type).tax.ym(year_index,age) - ps(type).tax.wm(year_index,age) )/ps(type).util.alpha(year_index,age))^(-ps(type).util.rho(year_index,age))*C(age); 

          % l_alt(age) = (wstar_alt(age)*(1-ps(type).tax.ym(year_index,type) - ps(type).tax.wm(year_index,type) )/ps(type).util.alpha(year_index))^(-ps(type).util.rho(year_index))*ps(type).opt.C(year_index,age); 

        end
        
    else
        
         l = zeros(ps.demog.lifespan,1);
    
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



            shadow_wages(age) = (ps(type).util.alpha(year_index,age)/(1-ps(type).tax.ym(year_index,age) - ps(type).tax.wm(year_index,age) ))*(C(age))^(1/ps(type).util.rho(year_index,age)) - ps(type).prod_af(age) *ps(type).opt.wprod(year_index,age);

            if shadow_wages(age) < 0

                shadow_wages(age) = 0;

            end


            % Update shadow wages with a dampening
            shadow_wages(age) =  run_schedule.swdamp*shadow_wages(age) + (1-run_schedule.swdamp)*ps(type).opt.sw(year_index,age);

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
    
    bequest_age = ps(type).util.a_bg(year_born);
    % In order to calculate bequests, we need to know the H
    if ps.endogenous_labor ==1
        
        nu_beq = (1+ps.util.alpha(bequest_year,bequest_age)*(ps.prod_af(bequest_age)*ps.opt.wstar(bequest_year,bequest_age)*(1-ps.tax.ym(bequest_year,bequest_age) - ps.tax.wm(bequest_year,bequest_age) )/ps.util.alpha(bequest_year,bequest_age))^(1-ps.util.rho(bequest_year,bequest_age)))^((ps.util.rho(bequest_year,bequest_age)-ps.util.gamma(bequest_year,bequest_age))/(1-ps.util.rho(bequest_year,bequest_age)));
        H_beq = nu_beq^(1/ps.util.gamma(bequest_year,bequest_age));

    else
        
        H_beq = 1;
        
    end


    bg = C(ps(type).util.a_bg(year_born)) ...
         * ps(type).util.mu(year_born)^(ps(type).util.gamma(bequest_year,ps(type).util.a_bg(year_born))) ...
         * ps(type).demog.num_kids(year_born)^(-ps(type).util.gamma(bequest_year,ps(type).util.a_bg(year_born))) ...
         * H_beq^(-ps(type).util.gamma(year,initial_age));

    br = bg;
 
    
    % Make adjustments because of productivity!!!!!
    % We have already made adjustments to the PVI. 
    % But now, we need to realize that people born in previous generations
    % had different productivities, and thus we need to adjust their
    % optimal quantities down by a factor
    
    
    if year==1 || ((year_min==year_max) && (year_min==(policy.nt+2)))

        for age=initial_age:ps(type).demog.lifespan
            
            %prod_af = 1/( (1+economy.lprod_growth_iss)^(1/(1-prod.epsilon(1))) )^(age - 1);

            C(age) =C(age)*(1/ps(type).prod_af(age));
            
            a(age) = a(age)*(1/ps(type).prod_af(age));
            
        end
        
        % Adjust bequests received for productivity
        % It is always your parents who give you bequests; so the
        % productivity difference is the age of the parents
        % Note that it may or may not be necessary to adjust bequests
        % given, but I am doing it anyway for now.

        % There is probably the additional division by the adjustment
        % factor b/c of the difference between generations
        
        ps(type).opt.br(year,ps(type).util.a_br(year)) = ps(type).opt.br(year,ps(type).util.a_br(year))*(1/ps(type).prod_af(ps(type).util.a_bg(year))) * (1/ps(type).prod_af(2));
        
        ps(type).opt.bg(year,ps(type).util.a_bg(year)) = ps(type).opt.bg(year,ps(type).util.a_bg(year))*(1/ps(type).prod_af(ps(type).util.a_bg(year)));

    end
    
end

