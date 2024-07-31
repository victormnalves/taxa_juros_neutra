function [ nipa_bank ] = nipa(nipa_bank,ps,prod,gov,policy,economy,prices,run_schedule,year,adj)
% Calculate NIPA Statistics which might be useful

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1: Actual NIPA Statistics %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Account 1                      %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Domestic Income
        
            a1 = zeros(36,1);
            % l1 Compensation of employees, paid
            a1(1) = economy.ag.wbase(year);

                % l2 wages and salaries
                a1(2) = economy.ag.wbase(year);

                    % l3 domestic
                    a1(3) = economy.ag.wbase(year);

                    % l4 rest of the world
                    a1(4) = 0;

                % Supplements to wages and salaries
                a1(5) = 0;

            % l6 Taxes on Production and Imports
            a1(6) = 0;

            % l7 Less subsidies
            a1(7) = 0;

            % l8 Net Operating Surplus
            % We only want to include interest income -- the rbase also
            % includes some of the stuff from dying, which isn't really
            % interest income; 
            

                % l9 private enterprises
                a1(9) = economy.ag.K(year)*prices.r(year) + economy.ag.profit(year);

                % l10 Current surplus of govenrment enterprises    
                a1(10) = 0; % will have to adjust this
                
            a1(8) = a1(9) + a1(10) ;

            % l11 Consumption of fixed capital
            a1(11) = prod.deprec(year)*economy.ag.K(year);

            % l12 Gross Domestic Income
            a1(12) = a1(1) + a1(6) + a1(7) + a1(8) + a1(11);

            % l13 Statistical Discrepancy
            a1(13) = 0;

            % l14 GDP
            a1(14) = a1(12) + a1(13);

        % Product Account

            % l15 Consumption
            a1(15) = economy.ag.C(year);
            
            % l16 Goods
            a1(16) = economy.ag.C(year);
            
            % l17 Durable goods
            a1(17) = 0;
            
            % l18 nondurable goods
            a1(18) = economy.ag.C(year);
            
            % l19 services
            a1(19) = 0;

            % Investment
            if year==1

                % l20 Investment 
                a1(20) = economy.ag.K(year)*(prod.deprec(year) + economy.n(year) + economy.ag.AL_growth_iss + economy.n(year)*economy.ag.AL_growth_iss + economy.ag.A_growth_adj_iss + economy.n(year)*economy.ag.A_growth_adj_iss  );

            elseif year==(policy.nt+2)
                
                a1(20) = economy.ag.K(year)*(prod.deprec(year) + economy.n(year) + economy.ag.AL_growth_fss + economy.n(year)*economy.ag.AL_growth_fss + economy.ag.A_growth_adj_fss + economy.n(year)*economy.ag.A_growth_adj_fss  );

            else

                a1(20) = economy.ag.K(year+1) - economy.ag.K(year) + prod.deprec(year)*economy.ag.K(year);
            end
            
            % l21 fixed investment
            a1(21) = a1(20);
            
            % l22 Nonresidential
            a1(22) = a1(20);
            
            % l23 Structures
            a1(23) = 0;
            
            % l24 Equipment
            a1(24) = a1(20);
            
            % l25 Intelectual property
            a1(25) = 0;
            
            % l26 residential
            a1(26) = 0;
            
            % l27 change in private inventories
            a1(27) = 0;
            
            % l28 net exports
            
            % Essentially, we assume that foreign governments use their
            % capital income in order to satisfy economy.nfa, and ship the
            % rest of income home to be consumed. 
            
            foreign_capital = economy.nfa(year)*economy.ag.K(year);
            foreign_income = prices.rentk(year)*foreign_capital;
            
            if year==1
                
                % If year is 1, net foreign assets grow at (1+n)(1+g)
                                
                foreign_investment = foreign_capital*(prod.deprec(year) + economy.n(year) + economy.ag.AL_growth_iss + economy.n(year)*economy.ag.AL_growth_iss + economy.ag.A_growth_adj_iss + economy.n(year)*economy.ag.A_growth_adj_iss);
                
                
                
            elseif year==(policy.nt+2)
                
                foreign_investment = foreign_capital*(prod.deprec(year) + economy.n(year) + economy.ag.AL_growth_fss + economy.n(year)*economy.ag.AL_growth_fss + economy.ag.A_growth_adj_fss + economy.n(year)*economy.ag.A_growth_adj_fss  );
                
                
            else
                
                foreign_investment =  economy.nfa(year+1)*economy.ag.K(year+1) -  economy.nfa(year)*economy.ag.K(year) + prod.deprec(year)* economy.nfa(year)*economy.ag.K(year);
            end
            
            foreign_consumption = foreign_income - foreign_investment;
            
            
            % l29 exports
            a1(29) = foreign_consumption*(foreign_consumption > 0);
            
            % l30 imports
            a1(30) = -1*foreign_consumption*(foreign_consumption < 0);
            
            % Now net exports
            a1(28) = a1(29) - a1(30) ;
            
            % l31 Government Spending
            a1(31) = gov.spend.amt(year);
            
            % l32 Federal
            a1(32) = a1(31);
            
            % l33 National defence
            a1(33) = 0;
            
            % l34 Nondefence
            a1(34) = 0;

            % GDP
            a1(36) = a1(15) + a1(20) + a1(28) + a1(31);

            nipa_bank.a1(year,:) = a1';
            
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Account 3                      %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        % Personal Income and Outlay Account
        a3 = zeros(26,1);
        
        % l10 Compensation of employees, received
        a3(10) = a1(1);
        
        % l17 Proprietor's income
        a3(17) = economy.ag.profit(year);
        
        % l20 personal interest income
        % This DOES NOT INCLUDE interest income from the annuity market
        a3(20) = economy.ag.personal_pure_interest_income(year);
        
        % l22 Transfer inclues ss payments, bequests
        % And especially, annuity payments from people who have died
            % Annuity payments (non interest) are calculated as total personal interest
            % income - pure interest income
            
            economy.ag.pure_annuity_payments(year) = ( (economy.ag.personal_interest_income(year) - economy.ag.personal_interest_payments(year) ) - a1(8) );
            economy.ag.pure_annuity_receipts(year) = economy.ag.pure_annuity_payments(year);
            
        a3(22) = gov.total_ss_pay(year) + economy.ag.personal_transfer_receipts(year) + economy.ag.pure_annuity_receipts(year) ; 
        
        % l25 Less contributions for government social insurance. We have here
        % social security taxes
        a3(25) = gov.total_ss_tax(year);
        
        % l26 Personal Income
        a3(26) = a3(10) + a3(17) + a3(20) + a3(22) - a3(25);
        
        % l1 Personal current taxes
        a3(1) = gov.tax.tax_rev_req(year);
        
        % l3 personal consumption
        a3(3) = a1(15);
        
        % l4 personal interest payments (pure interest)
        a3(4) = -economy.ag.personal_pure_interest_payments(year) ;
        
        % l5 personal transfer payments
        a3(5) = economy.ag.personal_transfer_payments(year) + economy.ag.pure_annuity_payments(year) ;
        
        % l2 personal outlays
        a3(2) = a3(3) + a3(4) + a3(5);
        
        % l8 personal savings
        a3(8) = a3(26) - (a3(1) + a3(2));
        
        % l9 personal taxes, outpays and savings
        a3(9) = a3(1) + a3(2) + a3(8);
        
        nipa_bank.a3(year,:) = a3';
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Other Statistics               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Savings Rates, Etc             %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %personal_savings_rate = 1 - (economy.ag.C(year)) / (economy.ag.wbase(year)*(1-ps.tax.wa(year)) + economy.ag.rbase(year) );
        nipa_bank.personal_savings_rate(year) = a3(8) / a3(9);
        
        % This is, in fact, nominal K/Y, to align ourselves with the data
        nipa_bank.KY(year) = economy.relp(year)*economy.ag.K(year) / economy.ag.Y(year);
        
        % Consumption per person
        nipa_bank.con_pp(year) = economy.ag.C(year)/economy.ag.pop(year);
        
        % Consumption per worker
        nipa_bank.con_pw(year) = economy.ag.C(year)/economy.ag.L(year);
        
        % Output per person
        nipa_bank.Y_pw(year) = economy.ag.Y(year)/economy.ag.L(year);
        
        % Consumption per person
        nipa_bank.Y_pp(year) = economy.ag.Y(year)/economy.ag.pop(year);
        
        % Output per worker
        
        
        % Capital to Labor Ratio
        nipa_bank.KL(year) = economy.ag.K(year)/economy.ag.L(year);
        
        % Productivity Adjusted Capital to Labor Ratio
        
        nipa_bank.KL_adj(year) = (economy.ag.AK(year)*economy.ag.K(year))/(economy.ag.AL(year)*economy.ag.L(year));
        
        % Output per worker
        
        nipa_bank.y_pw(year) = economy.ag.Y(year)/economy.ag.L(year);
        
        % Golden Rule 
        
        nipa_bank.golden_rule(year) = (prod.epsilon(year)/(economy.n(year) + prod.deprec(year)))^(1/(1-prod.epsilon(year))) ;
        
        % Support Ratio
        
        % nipa_bank.support_ratio(year) = economy.ag.L(year)/economy.ag.pop(year);
        % We will do the non-productivity adjusted support ratio
        
            num_workers = 0;

            for type=1:economy.num_types
                
                for age=1:40 % currently people work until they are 40

                    num_workers = num_workers + ps(type).demog.pop(year,age);

                end
                
            end
            
            nipa_bank.support_ratio(year,1) = num_workers / economy.ag.pop(year);
            
        % Average age of population
            % Percent of people age 1-18, 19-36, 37-56
            
        
            sum_age = 0;
            sum_120 = 0;
            sum_2140 = 0;
            sum_4156 = 0;

            for type=1:economy.num_types
                
                for age=1:ps(type).demog.lifespan

                    sum_age = sum_age + ps(type).demog.pop(year,age)*age;
                    
                    if age>=1 && age<=20
                        
                        sum_120 = sum_120 + ps(type).demog.pop(year,age);
                        
                    elseif age>=21 && age<=40
                        
                        sum_2140 = sum_2140 + ps(type).demog.pop(year,age);
                        
                    elseif age>=41 && age<=56
                        
                        sum_4156 = sum_4156 + ps(type).demog.pop(year,age);
                        
                    end

                end
                
            end
            
            nipa_bank.average_age(year,1) = sum_age / economy.ag.pop(year);
            
            nipa_bank.percent_120(year,1) = sum_120 / economy.ag.pop(year);
            
            nipa_bank.percent_2140(year,1) = sum_2140 / economy.ag.pop(year);
            
            nipa_bank.percent_4156(year,1) = sum_4156 / economy.ag.pop(year);
            
        % Population Growth Rate
        
        if year==1
            
            nipa_bank.pop_growth(year,1) = economy.n(1);
            
        elseif year==(policy.nt+2)
            
            nipa_bank.pop_growth(year,1) = economy.n(policy.nt+2);
            
        else
            
            nipa_bank.pop_growth(year,1) = (economy.ag.pop(year+1) / economy.ag.pop(year))-1;
            
        end
                
        % Bequests as % of Cap Stock
        
        nipa_bank.beq_cap(year) = economy.ag.brt(year)/economy.ag.K(year);
        
        % Bequests as % of Income
        
        nipa_bank.beq_inc(year) = economy.ag.brt(year)/economy.ag.Y(year);

        % Social Security as % of Income
        
        nipa_bank.ss_inc(year) = gov.total_ss_pay(year)/economy.ag.Y(year);
        
        % Personal Debt as % of Income
        
        nipa_bank.debt_inc(year) = -economy.ag.personal_debt(year)/economy.ag.Y(year);
        
        % Personal Debt as % of Capital Stock
        
        nipa_bank.debt_cap(year) = -economy.ag.personal_debt(year)/economy.ag.K(year);

        % Investment to output ratio
        
        nipa_bank.IY(year,1) = a1(20)/economy.ag.Y(year);
        
        % Nominal investment to output
        
        % Year Minus 1 variable;
        if year==1 || year==(policy.nt+2)
            
            ym1 = year;
            
        else
            
            ym1 = year - 1;
            
        end
        
        nipa_bank.IY_nom(year,1) = economy.relp(ym1)*a1(20)/economy.ag.Y(year);
        
    % Creating Gini Index of Wage Income

        total_wage_vector = [];
        total_savings_vector = [];

        for type = 1:economy.num_types

            for age = 1:ps(type).demog.lifespan

                % Need 105 values in order to get a vector that is 100
                share = round(1000*ps(type).demog.pop(year,age) / economy.ag.pop(year));

                total_wage_vector = [total_wage_vector ; (ps(type).opt.inc(year,age) + ps(type).opt.cap_inc(year,age))*ones(share,1)] ; 

                total_savings_vector = [total_savings_vector ; (ps(type).opt.a(year,age))*ones(share,1)] ;

            end

        end


        % Income Gini
        total_wage_vector = sort(total_wage_vector);

        lorenz_income = cumsum(total_wage_vector)./sum(total_wage_vector);
        nipa_bank.income_gini(year) = (.5 - sum(lorenz_income/(length(total_wage_vector))) ) / .5;

        % Savings Gini
        
        total_savings_vector = sort(total_savings_vector);

        lorenz_savings = cumsum(total_savings_vector)./sum(total_savings_vector);
        nipa_bank.wealth_gini(year) = (.5 - sum(lorenz_savings/(length(total_savings_vector))) ) / .5;

        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Income Shares                  %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
        % Labor Share
        nipa_bank.inc_share(year,1) = 100*a1(1) / a1(14);
        
        % Capital Income
        nipa_bank.inc_share(year,2) = 100*(economy.ag.K(year)*prices.r(year) + a1(11)) / a1(14);
        
        % Profit Share
        nipa_bank.inc_share(year,3) = 100*economy.ag.profit(year) / a1(14);
        
        % Also, calculate total average return on capital
        
        % Should we multiply depreciation by the relative price of capital?
        % JR -- I THinks so! Because how much are you getting as return?
        % You are getting K*MPK or revenue; and your cost is the
        % replacement cost, which is relp*delta*K
        nipa_bank.average_return(year,1) = 100*(economy.ag.profit(year) + economy.ag.K(year)*(prices.mpk(year) - prod.deprec(year)*economy.relp(year)))/(economy.relp(year)*economy.ag.K(year));
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Bequest Account                %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ym1_index = year + 1;
    prod_adj = 1;
    pop_adj = 1;

    if year==1

        ym1_index = 1;
        prod_adj = (1+economy.ag.AL_growth_iss)*(1+economy.ag.A_growth_adj_iss);
        pop_adj = (1+economy.n(1));
        
    elseif year==(policy.nt+2)

        ym1_index = policy.nt+2;
        prod_adj = (1+economy.ag.AL_growth_fss)*(1+economy.ag.A_growth_adj_fss);
        pop_adj = (1+economy.n(policy.nt+2));
    end
    
    % Intentional Bequests Received
    nipa_bank.beq_account(year,1) = economy.ag.brt_i(year);
    
    % Unintentional Bequests Received
    nipa_bank.beq_account(year,2) = economy.ag.brt_u(year);
    
    % Total Bequests Received
    nipa_bank.beq_account(year,3) = economy.ag.brt(year);
    
    % Intentional Bequests Given Last Period
    nipa_bank.beq_account(year,4) = economy.ag.bgt_i(ym1_index)*(1/prod_adj)*(1/pop_adj);
    
    % Unintentional Bequests Given Last Period
    nipa_bank.beq_account(year,5) = economy.ag.bgt_u(ym1_index)*(1/prod_adj)*(1/pop_adj);
    
    % Total Bequests Given Last Period
    nipa_bank.beq_account(year,6) = economy.ag.bgt(ym1_index)*(1/prod_adj)*(1/pop_adj);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Annuity Account                %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Annuity Revenue from individuals
    nipa_bank.ann_account(year,1) = economy.ag.annuity_rev_ind(year);
    
    % Annuity Revenue from firms
    nipa_bank.ann_account(year,2) = economy.ag.annuity_rev_firm(year);
    
    % Total Annuity Revenue
    nipa_bank.ann_account(year,3) = nipa_bank.ann_account(year,1) + nipa_bank.ann_account(year,2);
    
    % Total Annuity Cost
    nipa_bank.ann_account(year,4) = economy.ag.annuity_cost(year);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Government Account             %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Spending as % of GDP
    
        nipa_bank.gov_account(year,1) = gov.spend.amt(year) / economy.ag.Y(year);
        
    % Interest on Debt as % of GDP
    
        nipa_bank.gov_account(year,2) = (gov.debt(year)*economy.ag.K(year)*prices.r_gov(year)) / economy.ag.Y(year);
    
    % Total Spending
    
        nipa_bank.gov_account(year,3) = nipa_bank.gov_account(year,1) +  nipa_bank.gov_account(year,2);
        
    % Tax Revenue as % of GDP
    
        nipa_bank.gov_account(year,4) = gov.tax.tax_rev_req(year) / economy.ag.Y(year);
    
    % Deficit (% of GDp)
    
        nipa_bank.gov_account(year,5) = nipa_bank.gov_account(year,3)  - nipa_bank.gov_account(year,4);
         
    % Total
    
        nipa_bank.gov_account(year,6) = nipa_bank.gov_account(year,4) + nipa_bank.gov_account(year,5);
        
    % Debt as % of GDP
    
        nipa_bank.gov_account(year,7) = (gov.debt(year)*economy.ag.K(year)) / economy.ag.Y(year);
     
end

