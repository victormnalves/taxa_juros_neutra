function [economy,keep_running] = repeatfunc(ps,prod,gov,policy,prices,economy,run_schedule,tm)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Updates aggregate quantities and determines whether or not we will do %
% another iteration                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

keep_running = 0; % flag for whether we will repeat or not

K_ag = zeros(policy.nt +2,1); % Temporary storage for agg K
L_ag = zeros(policy.nt +2,1); % Temporary storage for agg Labor

% For partial equilibrium, get capital supplied
K_supply = zeros(policy.nt+2,1);

for year=tm.year_min: tm.year_max
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Step 1: Capital and Labor Convergence   %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Also, calculate aggregate consumption
        
        economy.ag.C(year) = 0;
        
        % And Aggregate Population
        
        economy.ag.pop(year) = 0;
        
        % Aggregate Bequests
        
        economy.ag.brt(year) = 0;
        
        economy.ag.brt_i(year) = 0;
        economy.ag.brt_u(year) = 0;
        
        economy.ag.bgt(year) = 0;
        
        economy.ag.bgt_i(year) = 0;
        economy.ag.bgt_u(year) = 0;
        
        % Annuity Firm Revenue
        
        economy.ag.annuity_rev_ind(year) = 0;
        economy.ag.annuity_rev_firm(year) = 0;
        economy.ag.annuity_rev_gov(year) = 0;
        
        economy.ag.annuity_cost(year) = 0;
        
        % Personal debt
        
        economy.ag.personal_debt(year) = 0;
        
        % Savings
        
        economy.ag.savings(year) = 0;
        
        % Loop through types
        for type=1:economy.num_types

            % Loop through ages
            for age=1:ps(type).demog.lifespan

                % Add up assets, multiply by population
                % One more thing. When individuals die, they still have
                % capital -- it is simply held by the annuity firms. 
                K_ag(year) = K_ag(year) + ps(type).opt.a(year,age) * ps(type).demog.pop(year,age)*(1/ps(type).demog.mhs(year,age));

                % Add up labor; multiply hours by human capital and hte
                % population
                L_ag(year) = L_ag(year) + (1-ps(type).opt.l(year,age))*ps(type).demog.hc(year,age) * ps(type).demog.pop(year,age);

                % Multiply consumption by population
                economy.ag.C(year) =  economy.ag.C(year) + ps(type).opt.C(year,age) * ps(type).demog.pop(year,age);
                
                % Get Population
                
                economy.ag.pop(year) = economy.ag.pop(year) + ps(type).demog.pop(year,age);
                
                % Get Savings
                
                economy.ag.savings(year) = economy.ag.savings(year) + ps(type).demog.pop(year,age)*ps(type).opt.savings(year,age);
                
                % Add up debt in the economy
                
                economy.ag.personal_debt(year) = economy.ag.personal_debt(year) + (ps(type).opt.a(year,age) < 0)*ps(type).opt.a(year,age) * ps(type).demog.pop(year,age)*(1/ps(type).demog.mhs(year,age));
                                
                % Get Annuity Revenue / costs
                
                % The annuity firm has revenue from lending to individuals
                % at the annuity rate. It only gets paid if the individual
                % survives, of course. 
                
                if ps(type).opt.a(year,age) < 0

                    % WARNING -- I may have to change the index on annuity
                    % pariticpation if I ever do a transition with this
                    % variable
                    
                    economy.ag.annuity_rev_ind(year) = economy.ag.annuity_rev_ind(year) + -1*ps(type).util.annuity_participation(year,age)*ps(type).opt.a(year,age) * ( ps(type).demog.pop(year,age) / ps(type).demog.mhs(year,age))*ps(type).demog.mhs(year,age)*(1+ps(type).prices.r_annuity(year,age));
                
                    economy.ag.annuity_cost(year) = economy.ag.annuity_cost(year) + -1*ps(type).util.annuity_participation(year,age)*ps(type).opt.a(year,age) * ( ps(type).demog.pop(year,age) / ps(type).demog.mhs(year,age))*(1+prices.r(year));
                else 
                    
                    % The annuity firm takes the assets from the
                    % individuals and lends to the firm, getting revenue of
                    % the risk free rate
                    economy.ag.annuity_rev_firm(year) = economy.ag.annuity_rev_firm(year) + ps(type).util.annuity_participation(year,age)*ps(type).opt.a(year,age) * ( ps(type).demog.pop(year,age) / ps(type).demog.mhs(year,age))*(1+prices.r(year));
                    
                    % With positive assets, the annuity firm is borrowing
                    % from the individual. Thus the cost is the interest
                    % payments, less the probability the individual dies
                
                    economy.ag.annuity_cost(year) = economy.ag.annuity_cost(year) + ps(type).util.annuity_participation(year,age)*ps(type).opt.a(year,age) * ( ps(type).demog.pop(year,age) / ps(type).demog.mhs(year,age))*(ps(type).demog.mhs(year,age))*(1+ps(type).prices.r_annuity(year,age));
                    
                end
                
            end
                
            K_supply(year) = K_ag(year);
            
            % Bequests
            
            economy.ag.brt(year) = economy.ag.brt(year) + ps(type).opt.brt(year);
            
            economy.ag.brt_i(year) = economy.ag.brt_i(year) + ps(type).opt.brt_i(year);
            economy.ag.brt_u(year) = economy.ag.brt_u(year) + ps(type).opt.brt_u(year);
            
            economy.ag.bgt(year) = economy.ag.bgt(year) + ps(type).opt.bgt(year);
            
            economy.ag.bgt_i(year) = economy.ag.bgt_i(year) + ps(type).opt.bgt_i(year);
            economy.ag.bgt_u(year) = economy.ag.bgt_u(year) + ps(type).opt.bgt_u(year);

        end
        
        % Adjust Capital for Government Debt alt
        
       %K_ag(year) = K_ag(year) - gov.debt_alt(year);
       
       % Adjust Capital for foreign demand for capital
        
        % Add up adjustments from transfer authority
        
        
        % The debt that the LSRA issues must be subtracted from the capital
        if year==2
            
            economy.ag.vdebt(year) = 0;
            
        % At the end of period 2, the authority
        % issues XG + V(2) of debt to individuals to pay for the transfers
        % of XG + V(2); this is the debt at the beginning of period 3, and
        % must be subtracted from capital
        elseif year==3
            
            economy.ag.vdebt(year) = economy.ag.xg;
            
             % Loop through types
            for type=1:economy.num_types
                
                economy.ag.vdebt(year) = economy.ag.vdebt(year) + ps(type).opt.V(year-1)*ps(type).demog.pop(year-1,1);
            
            end
            
        elseif year > 3
            
            economy.ag.vdebt(year) = economy.ag.vdebt(year-1)*(1+prices.r_gov(year-1)) ;
            
             % Loop through types
            for type=1:economy.num_types
                
                economy.ag.vdebt(year) = economy.ag.vdebt(year) + ps(type).opt.V(year-1)*ps(type).demog.pop(year-1,1);
                
            end
            
        end
        
        % Alternate Debt
        
       % for t=2:(year-1
        
        %%%% 
        % Error Catch for Debt from Transfer Authority
        %%%%
        
        % if lsra debt is greater than 1/2 of capital stock, adjust it
        if abs( economy.ag.vdebt(year) / economy.ag.K(year) ) > 1
            
            economy.ag.vdebt(year) = economy.ag.K(year)*1*sign(economy.ag.vdebt(year));
            
            disp(['Warning: LSRA debt is too high in year '  num2str(year)])

        end
        
        % We don't want to double count transfers; thus we need to subtract
        % them from capital
        economy.ag.vc(year) = 0;
        
        for type=1:economy.num_types
            
            economy.ag.vc(year) = economy.ag.vc(year) + ps(type).opt.V(year) ...
                                  * ps(type).demog.pop(year,1) ...
                                  / (1 + prices.r_gov(year)*(1-ps(type).tax.ya(year,1) - ps(type).tax.ka(year,1) ) );
                              
        end
        
        
        % Adjust VC for XG
        
        if year==2
            
             economy.ag.vc(year) = economy.ag.vc(year) + economy.ag.xg / (1 + prices.r_gov(year)*(1-ps(type).tax.ya(year,1) - ps(type).tax.ka(year,1) ));
            
        end
        
        % Adjust Capital for Q and Debt and LSRA
        
        
            % Set NFA
            
            if run_schedule.nfa_gdp==1 % i.e. we specify net foreign assets are a % of gdp
                
                economy.nfa(year) = economy.nfa_gdp(year)*economy.ag.Y(year)/economy.ag.K(year);
                
            end
        
        % We get this from the following formula:
            % Total Assets Held by Residents = K + gov.debt*K - nfa*K
            % Thus K = Total Assets by Residents / (1 + gov.debt - nfa)
            
            % Now, how about the relative price of capital? How does this
            % change things? 
            
            
            % Remember, in this program we think of the individual only
            % holding bonds that pay the real interest rate, not capital.
            % The equilibrium condition is that the value of bonds the
            % individuals holds equals the value of the assets in the
            % economy
            
            % Thus 
            % Value of Assets Held by Residents = relp*K + gov.debt*K -
            % nfa.K 
            
            % Notice that we store government debt as gov.debt*K, but we do
            % not multiply by the relative price :)
            
            % Thus K = VAL/(relp + gov.debt - nfa)
            
        % Year Minus 1 variable;
        if tm.year_min==tm.year_max
            
            ym1 = year;
            
        else
            
            ym1 = year - 1;
            
        end
        
        % Question -- why is this year minus 1?
        % Answer -- because the capital choice at the beginning of period t
        % was made in period t-1 :)
        K_ag(year) = (K_ag(year) - economy.ag.vc(year) - economy.ag.vdebt(year) )/(economy.Q(year)*economy.relp(ym1) + gov.debt(year) - economy.nfa(year));

        % Test for convernge in full equilibrium
        
        if run_schedule.pe.run==0
            % Test for convergence. If stop = 1, we stop
            keep_running = max( (abs(economy.ag.K(year) - K_ag(year)) > run_schedule.tol),keep_running);
            keep_running = max( (abs(economy.ag.L(year) - L_ag(year)) > run_schedule.tol),keep_running);

        else % In Partial Equilibrium
            
            keep_running = max( (abs(economy.ag.K_supply(year) - K_supply(year)) > run_schedule.tol),keep_running);
            
        end
        
        % Update Capital and Labor

        % General Equilibrium
        if run_schedule.pe.run==0
            
            economy.ag.K(year) = run_schedule.kldamp * K_ag(year) + (1-run_schedule.kldamp) * economy.ag.K(year);
        
            economy.ag.L(year) = run_schedule.kldamp * L_ag(year) + (1-run_schedule.kldamp) * economy.ag.L(year);

        
        else % Partial
            
            economy.ag.K_supply(year) = K_supply(year);
        
            economy.ag.L(year) = L_ag(year);
        end
            
        % Now, add the revenue of the annuity company to lending to the firms
        % Eventually, may have to add in financial frictions
        %economy.ag.annuity_rev_firm(year) = economy.ag.annuity_rev_firm(year) + (1+prices.r(year))*K_ag(year);

        % Now, add the revenue of the annuity firm lending to the
        % government
        
        %economy.ag.annuity_rev_gov(year) = economy.ag.annuity_rev_gov(year) + (1+prices.r_gov(year))*gov.debt(year)*K_ag(year);
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Step 2: Check for Shadow Wage Convergence   %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        for type=1:economy.num_types

            for age=1:ps(type).demog.lifespan

                if (ps(type).opt.sw(year,age) > run_schedule.tol) && ( ps(type).opt.l(year,age) > 1 + run_schedule.tol)

                    keep_running = 1;

                end

            end

        end
        

end % year

end

