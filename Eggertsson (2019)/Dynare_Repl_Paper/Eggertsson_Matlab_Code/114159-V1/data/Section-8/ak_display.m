function [ output] = ak_display(nipa_bank,ps,prod,gov,policy,economy,prices,run_schedule,year,adj)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ak_display
%
% Displays output as in AK 1987 Models

year = year + 1;
adjustment = 1;

if adj==1
    
    adjustment = 1/(1+economy.n(year))^(year-1);
    
end
% Line 1
disp(['Year  ' 'Capital  ' 'Labor  ' 'Consumption  ' 'Utility'])
disp([num2str(year-1) '  ' num2str(economy.ag.K(year)*adjustment) '  ' num2str(economy.ag.L(year)*adjustment) '  ' num2str(economy.ag.C(year)*adjustment) '  ' num2str(ps(1).opt.usave(year))])

disp (['Wage  ' 'Interest  ' 'Revenue  ' 'Income  ' 'Save Rate  '])
disp([num2str(prices.wages(year)) '  ' num2str(prices.r(year)) '  ' num2str(gov.tax.tax_rev_req(year)*adjustment) '  ' num2str(economy.ag.Y(year)*adjustment) '  ' ])

disp (['Debt  ' 'Vdebt  ' 'Rep Rate  ' 'SS Tax  ' 'Transfer'])
disp([num2str(gov.debt(year)) '      ' num2str(economy.ag.vdebt(year)*adjustment)  '      ' num2str(gov.rep(year)) '      ' num2str(gov.tax.sst(year)) '      ' num2str(ps(1).opt.V(year))])


disp(['Ret Age  ' 'Con Tax  ' 'Wage Tax  ' 'Inc Tax  ' 'Cap Tax'])
disp([num2str(1) '         ' num2str(0) '      ' num2str(gov.tax.wp(year)) '     ' num2str(gov.tax.yp(year)) '        ' num2str(gov.tax.kp(year)) '  '])

disp(['Bequest Level    '  'Bequests Pct Cap Stock      ' 'Bequests Pct Income'])
disp([num2str(ps(1).opt.br(year,ps(1).util.a_br(year))*ps(1).demog.num_kids(year)) '   ' num2str(nipa_bank.beq_cap(year)) '    ' num2str(nipa_bank.beq_inc(year))])

disp(['SS Pay Pct GDP    ' 'Debt Pct Income    ' 'Debt Pct Cap'  ])
disp([num2str(nipa_bank.ss_inc(year)) '  ' num2str(nipa_bank.debt_inc(year)) '  ' num2str(nipa_bank.debt_cap(year))])

disp(['Income Gini    '  'Weath Gini'])
disp([num2str(nipa_bank.income_gini(year)) '   ' num2str(nipa_bank.wealth_gini(year))])

disp(['Nom Investment Output    '  'Labor Share'])
disp([num2str(nipa_bank.IY_nom(year,1)) '   ' num2str(nipa_bank.inc_share(year,1))])

disp(['Rental K    '  'MPK        ' 'Pop Growth'])
disp([num2str(prices.rentk(year,1)) '   ' num2str(prices.mpk(year,1)) '   ' num2str(economy.n(year,1))])



disp(['   '])
disp(['   '])
disp(['Consumption per Person  '  'Consumption per Worker  ' 'Capital Labor  ' 'Output to Labor  ' ])
disp([num2str(nipa_bank.con_pp(year)) '   ' num2str(nipa_bank.con_pw(year)) '   ' num2str(nipa_bank.KL(year)) '   ' num2str(nipa_bank.y_pw(year))])

disp(['   '])
disp(['   '])
disp(['Golden Rule KL    ' 'Support Ratio   ' 'Utility from Con' ])
disp([num2str(nipa_bank.golden_rule(year)) '    ' num2str(nipa_bank.support_ratio(year)) '    ' num2str(ps(1).opt.u_new(1)) ])



% Steady State Consumption GRaph
if year==1
    if any(strcmp('graph',fieldnames(run_schedule)))==1
    % Plot Consumption against Gourinchas & Parker Consumption
    
        c_gp = [20346.44 ;20608.55 ;20553.39 ;20568.67 ;20832.85 ;20960.67 ;21566.72 ;21459.77 ;21282.63 ;22127.91 ;21929.88 ;21460.42 ;21619.97 ;22359.83 ;22107.36 ;23016.14 ;22424.62 ;22871.31 ;23250.59 ;23839.26 ;22803.64 ;22548.49 ;23354.26 ;22359.36 ;21651.77 ;21383.05 ;21787.46 ;21454.16 ;20358.78 ;19842.86 ;20311.13 ;20353.93 ;19331.62 ;19082.18 ;17613.98 ;19077.07 ;18321.34 ;18501.98 ;17788.04 ;18201.76];
    
        
        % Adjustment Factor for consumption
        af_1 = mean(c_gp); 
        af_2 = mean(ps(1).opt.C_iss(1:40)); % 1 to 40 b/c these are non-retirement years
        new_af = af_1/af_2;
        
        close all
        
        hold on
        plot(26:65,new_af*ps(1).opt.C_iss(1:40),'LineWidth',2)
        plot(26:65,c_gp,'o','LineWidth',1)
        legend('Model Consumption','CEX Consumption','Location','northoutside')
        xlabel('Age')
        ylabel('Consumption, 1987 Dollars')
        ylim([0 30000])
        print('figures/iss_con_profile','-dpng') 
        hold off
        
        close all
        
    end
    
end

% Demographic Pydamid Graph
if any(strcmp('graph',fieldnames(run_schedule)))==1 % If the field exists
    
    if year==1

        % Import US Population

            us_pop_data =  readtable('data/us_census_pop.xlsx');
            us_pop_1970 = table2array(us_pop_data(:,'us_pop_1970'));   

            us_pop_2015 = table2array(us_pop_data(:,'us_pop_2015'));   

            % Clean Pop Data

            us_pop_1970 = us_pop_1970 / sum(us_pop_1970);

            us_pop_2015 = us_pop_2015 / sum(us_pop_2015);

            % Make the plot for 1970

                subplot(1,2,1)

                barh(26:81,ps.demog.pop(1,:)/economy.ag.pop(1),'r')
                legend('Model Population','Location','northoutside')
                xlabel('% of Total Population')
                ylabel('Age')

                subplot(1,2,2)
                barh(26:81,us_pop_1970,'b')
                legend('US 1970 Population','Location','northoutside')
                xlabel('% of Total Population')
                ylabel('Age')

                print('figures/iss_dem_pyramid','-dpng') 

                hold off

                close all

            % Make the plot for 2015

                subplot(1,2,1)

                barh(26:81,ps.demog.pop(1,:)/economy.ag.pop(1),'r')
                legend('Model Population','Location','northoutside')
                xlabel('% of Total Population')
                ylabel('Age')

                subplot(1,2,2)
                barh(26:81,us_pop_2015,'b')
                legend('US 2015 Population','Location','northoutside')
                xlabel('% of Total Population')
                ylabel('Age')

                print('figures/iss_dem_pyramid_2014','-dpng') 

                hold off

                close all

end
end
output = 1;

end

