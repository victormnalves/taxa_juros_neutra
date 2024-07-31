function [ ps,gov] = gov_dec(ps,prod,gov,policy,economy,prices,run_schedule,tm)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Gov Dec
% 
% The government is characterized by three decisions: taxes, spending, and
% the proportion of government revenue that is financed by the deficit
%
% Individuals specify two of the above as exogenous; the third is
% endogenously determined. 
%
% As it happens, government spending will never be endogenous -- the
% program so far is not set up to do this


% Loop Through Years

for year=tm.year_min:tm.year_max
    
    % JR Government spending. This sets government spending in a particular
    % way:
        % Government:
            % Simply Set Government Spending as a % of GDP. Accepts a
            % vector
        % Debt
            % During year 1, set debt to be a given % of GDP. Accepts a
            % vector. Alternatively, can be set to be a level of GDP
            
            % Year 2 should be set the same as year 1. This is because we
            % are initially in a steady state, and thus the debt cannot be
            % changed in year 2
            % 
            % Given the vector of debt, the deficits are chosen
            % endogenously to make sure debt equals the vector
            
    if gov.spend.jr==1
        
        % Set Government spending to a % of GDP
        gov.spend.amt(year) = gov.spend.amt_gdp(year)*economy.ag.Y(year);
        
        % Debt specified as level, not % of GDP
        % Specified for each year
        if gov.debt_alt(1) ~= 0
            
            gov.debt_gdp = gov.debt_alt./economy.ag.Y;
            
        end
        
        if year==1
            
            % First, calculate initial debt to set this at a level of GDP
            % given
            % Since gov.debt is debt as % of capital stock, need to do this
            gov.debt(1) = gov.debt_gdp(1)*economy.ag.Y(1)/economy.ag.K(1);
            
            % Now, in the initial steady state, the debt must grow at the
            % same proportion as GDP, which grows at the rate of
            % productivity and population growth
            
            % Let D2 be debt in period 2, and D1 be debt in period 1
            % The government's budget constraint is given as 
            
            %D2 = D1(1+r) + G - T
            
            % Now, as in the program below, we define REVT as total
            % governemtn revenue (which is equal to Taxes plus change in
            % debt). Since total revenue equals total spending (G + rD1),
            % we can rewrite
            
            % D2 = D1 + REVT - T
            
            % Now, from below, we have T = REVT*(1-DEFICIT), where deficit
            % is the proportion of total revenue funded by the deficit. 
            
            % Thus we have
            %D2 = D1 + REVT - REVT + REVT*DEF = D1 + REVT*DEF
            
            % In the state, D2/D1 = (1+g)(1+n). Thus, dividng by D1, we
            % have
            
            %(1+g)(1+n) = 1 + REVT*DEF/D1
            % Thus DEF = ( (1+g)(1+n) - 1)*D1/REVT
            
            % We thus set this below. 
            if gov.tax.revt(1)~=0 % on starting value, this will be 0. We don't want to create an error

                gov.deficit(1) = ( (1+economy.ag.AL_growth_iss)*(1+economy.ag.A_growth_adj_iss)*(1+economy.n(1)) - 1)...
                                 *(gov.debt(1)*economy.ag.K(1))/gov.tax.revt(1);

            end
            
        elseif year==(policy.nt+2) && tm.year_min==tm.year_max
            
            % First, calculate initial debt to set this at a level of GDP
            % given
            gov.debt(policy.nt+2) = gov.debt_gdp(policy.nt+2)*economy.ag.Y(policy.nt+2)/economy.ag.K(policy.nt+2);
            
            % Now, set the deficit, so that debt grows at the rate of
            % population growth + productivity growth
            if gov.tax.revt(policy.nt+2)~=0 % on starting value, this will be 0

                gov.deficit(policy.nt+2) = ( (1+economy.ag.AL_growth_fss)*(1+economy.ag.A_growth_adj_fss)*(1+economy.n(policy.nt+2)) - 1)...
                                 *(gov.debt(policy.nt+2)*economy.ag.K(policy.nt+2))/gov.tax.revt(policy.nt+2);

            end
            
        elseif year==(policy.nt+2) % We must do this as well along the transition
            
            % Now, set the deficit, so that debt grows at the rate of
            % population growth + productivity growth
            if gov.tax.revt(policy.nt+2)~=0 % on starting value, this will be 0

                gov.deficit(policy.nt+2) = ( (1+economy.ag.AL_growth_fss)*(1+economy.ag.A_growth_adj_fss)*(1+economy.n(policy.nt+2)) - 1)...
                                 *(gov.debt(policy.nt+2)*economy.ag.K(policy.nt+2))/gov.tax.revt(policy.nt+2);

            end
            
        else % Beyond initial year, calculate deficit so as to match the debt to GDP ratio next period
                        
            yp1 = year + 1;
            
            if gov.tax.revt(year)~=0
                
                % Recall that the government budget constraint is given by 
                % D2 = D1 + REVT*DEF
                
                % Dividing by Y2, we have
                % D2/Y2 = D1/Y2 + REVT*DEF / Y2
                
                % We are given a series for D2/Y2, thus we must set DEF to
                % match this
                
            if year==128 || year==129 || year==130
                
                
                eagle = 55;
                
            end
                
                % Thus DEF = (d2/Y2 - D1/Y2)*Y2/REVT
                gov.deficit(year) = ( gov.debt_gdp(yp1) - gov.debt(year)*economy.ag.K(year)/economy.ag.Y(yp1) )...
                                    * economy.ag.Y(yp1)/gov.tax.revt(year);
                                
                % My little corrections
                % Sometimes because of initial starting values this can get
                % a little out of wack, so we need this correction
                if gov.deficit(year) > .7 || gov.deficit(year) < -.7
                    
                    gov.deficit(year) = gov.deficit(1);
                    
                end

            end
            

        end
                
    end
    
    % AK Government spending == set spending equal to tax revenue collected
    % in initial steady state minus payments on debt
    if gov.spend.ak==1
        
        if year==1
            
            gov.spend.amt(year) = gov.tax.exog_rev(year) - ...
                            gov.debt(year)*economy.ag.K(year)*prices.r_gov(year);
                        
            % Add in non-proportional debt
            %gov.spend.amt(year) = gov.spend.amt(year) - gov.debt_alt(year)*prices.r_gov(year);

                        
        elseif (tm.year_min==tm.year_max) && tm.year_min==(policy.nt+2) % Final Steady State
            
            gov.spend.amt(year) = gov.spend.amt(1)*(1+economy.n(year))^(year-1)*(1+economy.ag.AL_growth_fss)^(year-1)*(1+economy.ag.A_growth_adj_fss)^(year-1);

        else % we simply update gov spending with population growth and productivity growth
             % Note the adjustment if we use Hicks Neutral, rather than
             % Harrod neutral (labor enhancing)
            
            gov.spend.amt(year) = gov.spend.amt(year-1)*(1+economy.n(year))*(1+economy.ag.AL_growth(year))*(1+economy.ag.A_growth(year))^(1/(1-prod.epsilon(year)));
            
        end
        
    end
    
    % Total Government Revenue = Taxes plus change in the debt. It is used to pay for 
    % Government spending plus interest payments
    gov.tax.revt(year) = gov.spend.amt(year) + ...
                         prices.r_gov(year)*gov.debt(year)*economy.ag.K(year);
    
    % Also add in the alternative formulation of debt, i.e., debt not as a
    % percentage of the capital stock
    
    %gov.tax.revt(year) = gov.tax.revt(year) + prices.r_gov(year)*gov.debt_alt(year);
                     
     
    % Now, calculate the total revenue that needs to be collected from
    % taxes
    gov.tax.tax_rev_req(year) = gov.tax.revt(year)*(1-gov.deficit(year));
    
    % Now, add up tax revenue from exogenous taxes
    
    % In order to do so, we need to calculate average tax rates
    
    % Initailize Exogenous tax variable
    gov.tax.exog_rev(year) = 0;
    
    % Loop through types
    for type = 1:economy.num_types 

        % Loop through ages
        for age=1:ps(type).demog.lifespan
            
            % Income Tax
            
                % Only calculate it if the tax rate is exogenous
                if (gov.tax.ype(year) ~= 1 && gov.tax.yge(year)~= 1)

                    ps(type).tax.ya(year,age) = run_schedule.tdamp * (gov.tax.yp(year) + ...
                                                gov.tax.yg(year)*...
                                                (ps(type).opt.inc(year,age) + ps(type).opt.cap_inc(year,age))/2) ... % Im Changing this now; used to be ps(type).opt.a(year,age)*prices.r(year)
                                                + (1-run_schedule.tdamp)*ps(type).tax.ya(year,age);

                    ps(type).tax.ym(year,age) = run_schedule.tdamp * (gov.tax.yp(year) + ...
                                                gov.tax.yg(year)*...
                                                (ps(type).opt.inc(year,age) + ps(type).opt.cap_inc(year,age) )) ... % used to be ps(type).opt.a(year,age)*prices.r(year)
                                                + (1-run_schedule.tdamp)*ps(type).tax.ym(year,age);

                    % Add to exogenous revenue from the income tax
                    gov.tax.exog_rev(year) = gov.tax.exog_rev(year) + ps(type).tax.ya(year,age)*(ps(type).opt.inc(year,age) + ps(type).opt.cap_inc(year,age) )*ps(type).demog.pop(year,age);

                end
                
            % Wage Tax
            
                % Only calculate it if the tax rate is exogenous
                if (gov.tax.wpe(year) ~= 1 && gov.tax.wge(year)~= 1)

                    ps(type).tax.wa(year,age) = run_schedule.tdamp * (gov.tax.wp(year) + ...
                                                gov.tax.wg(year)*...
                                                (ps(type).opt.inc(year,age))/2) ...
                                                + (1-run_schedule.tdamp)*ps(type).tax.wa(year,age);

                    % Year Born -- we need this because this is the index which ss_age
                    % is stored in

                    year_born = year - age + 1;
                    if year_born < 1
                        year_born = 1;
                    end
                    
                    % Question -- why do we take into account social
                    % security here? They seem to do it in AK as well.
                    % Mystery which will remain for now.
                    ps(type).tax.wm(year,age) = run_schedule.tdamp * (gov.tax.wp(year) + ...
                                                gov.tax.wg(year)*ps(type).opt.inc(year,age) ...
                                                + ps(type).tax.sst(year,age)*(age<gov.ss_age(year_born) )) ...
                                                + (1-run_schedule.tdamp)*ps(type).tax.wm(year,age);


                    % Add to exogenous revenue from the income tax
                    gov.tax.exog_rev(year) = gov.tax.exog_rev(year) + ps(type).tax.wa(year,age)*(ps(type).opt.inc(year,age))*ps(type).demog.pop(year,age);

                end
                
            % Capital Tax
            
                % Only calculate it if the tax rate is exogenous
                if (gov.tax.kpe(year) ~= 1 && gov.tax.kg(year)~= 1)

                    ps(type).tax.ka(year,age) = run_schedule.tdamp * (gov.tax.kp(year) + ...
                                                gov.tax.kg(year)*...
                                                (ps(type).opt.cap_inc(year,age))/2) ... % Changing it here
                                                + (1-run_schedule.tdamp)*ps(type).tax.ka(year,age);
                                        
                                            
                    ps(type).tax.km(year,age) = run_schedule.tdamp * (gov.tax.kp(year) + ...
                                                gov.tax.kg(year)*...
                                                (ps(type).opt.cap_inc(year,age))) ... % changing it here
                                                + (1-run_schedule.tdamp)*ps(type).tax.km(year,age);

                    % Add to exogenous revenue from the income tax
                    gov.tax.exog_rev(year) = gov.tax.exog_rev(year) + ps(type).tax.ka(year,age)*(ps(type).opt.cap_inc(year,age))*ps(type).demog.pop(year,age);

                    
                end
                
                
                
        end % Age
        
    end % Type
    
    % Calculate Residual Revenue Requirement
    gov.tax.resid_req(year) = gov.tax.tax_rev_req(year) - gov.tax.exog_rev(year);
    
    
    % Now, calculate endogenous tax rates! 
    
        % Wage
        
        if (gov.tax.wpe(year) == 1)
            
            gov.tax.wp(year) = (gov.tax.resid_req(year)*gov.tax.wrevpro(year) - gov.tax.wg(year)*economy.ag.wbasesq(year)*.5)/economy.ag.wbase(year);

            % Now, calculate average wage taxes
            
            % Loop through types
            for type = 1:economy.num_types 

                % Loop through ages
                for age=1:ps(type).demog.lifespan
            
                    ps(type).tax.wa(year,age) = run_schedule.tdamp * (gov.tax.wp(year) + ...
                                                gov.tax.wg(year)*...
                                                (ps(type).opt.inc(year,age))/2) ...
                                                + (1-run_schedule.tdamp)*ps(type).tax.wa(year,age);

                    % Year Born -- we need this because this is the index which ss_age
                    % is stored in

                    year_born = year - age + 1;
                    if year_born < 1
                        year_born = 1;
                    end

                    ps(type).tax.wm(year,age) = run_schedule.tdamp * (gov.tax.wp(year) + ...
                                                gov.tax.wg(year)*ps(type).opt.inc(year,age) ...
                                                + ps(type).tax.sst(year,age)*(age<gov.ss_age(year_born) )) ...
                                                + (1-run_schedule.tdamp)*ps(type).tax.wm(year,age);


                end
                
            end
            
        end
            
        % Income
        
        if (gov.tax.ype(year) == 1)
            
            gov.tax.yp(year) = (gov.tax.resid_req(year)*gov.tax.yrevpro(year) - gov.tax.yg(year)*economy.ag.ybasesq(year)*.5)/economy.ag.ybase(year);

            % Now, calculate average wage taxes
            
            % Loop through types
            for type = 1:economy.num_types 

                % Loop through ages
                for age=1:ps(type).demog.lifespan
            
                    ps(type).tax.ya(year,age) = run_schedule.tdamp * (gov.tax.yp(year) + ...
                                                gov.tax.yg(year)*...
                                                (ps(type).opt.inc(year,age) + ps(type).opt.cap_inc(year,age))/2) ... % Changing it
                                                + (1-run_schedule.tdamp)*ps(type).tax.ya(year,age);

                    ps(type).tax.ym(year,age) = run_schedule.tdamp * (gov.tax.yp(year) + ...
                                                gov.tax.yg(year)*...
                                                (ps(type).opt.inc(year,age) + ps(type).opt.cap_inc(year,age))) ... % Changing it
                                                + (1-run_schedule.tdamp)*ps(type).tax.ym(year,age);

                                            
                end
                
            end
            
        end
        
        % Capital Income
        
        if (gov.tax.kpe(year) == 1)
            
            gov.tax.kp(year) = (gov.tax.resid_req(year)*gov.tax.rrevpro(year) - gov.tax.kg(year)*economy.ag.rbasesq(year)*.5)/economy.ag.rbase(year);

            % Catch errors
            
            if gov.tax.kp(year) >= 1
                
                gov.tax.kp(year) = .5;
                
            end
            % Now, calculate average wage taxes
            
            % Loop through types
            for type = 1:economy.num_types 

                % Loop through ages
                for age=1:ps(type).demog.lifespan
            
                    ps(type).tax.ka(year,age) = run_schedule.tdamp * (gov.tax.kp(year) + ...
                                                gov.tax.kg(year)*...
                                                (ps(type).opt.cap_inc(year,age))/2) ... % Changing it
                                                + (1-run_schedule.tdamp)*ps(type).tax.ka(year,age);
                                        
                    ps(type).tax.km(year,age) = run_schedule.tdamp * (gov.tax.kp(year) + ...
                                                gov.tax.kg(year)*...
                                                (ps(type).opt.cap_inc(year,age))) ... % changing it
                                                + (1-run_schedule.tdamp)*ps(type).tax.km(year,age);

                end
                
            end
            
        end
        
        
        % Now, calculate debt!
        
        % We need to fill in debt for year 2
        
        if year==2
            
            gov.debt(2) = (gov.debt(1)*economy.ag.K(1)*(1+prices.r_gov(1)) + gov.spend.amt(1) - gov.tax.tax_rev_req(1))/economy.ag.K(2);
            
        end
        
            yp1 = year + 1;

            if year==129 || year==130
                
                
                eagle = 55;
                
            end
            
            if year < (policy.nt+2) % no need to calculate this beyond the final ss
                
                gov.debt(yp1) = (gov.debt(year)*economy.ag.K(year)*(1+prices.r_gov(year)) + gov.spend.amt(year) - gov.tax.tax_rev_req(year))/economy.ag.K(yp1);
           
            end
            
    
end

end

