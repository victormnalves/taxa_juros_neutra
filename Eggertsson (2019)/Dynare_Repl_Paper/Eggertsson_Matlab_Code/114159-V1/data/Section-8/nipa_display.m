function [ output ] = nipa_display(nipa_bank,ps,prod,gov,policy,economy,prices,run_schedule,year,adj)
% Calculate NIPA Statistics which might be useful

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1: Display NIPA Tables    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Account 1                      %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % we want to adjust this for GDP , thus we will devide by line 36. 
    nipa_bank.a1(year,:) = 100*nipa_bank.a1(year,:)./nipa_bank.a1(year,36);
    
    disp(['Account 1. Domestic Income and Product Account'])
    disp([' '])
    disp(['For the year ' num2str(year)])
    disp([' '])
    disp(['Values expressed as a % of GDP'])
    disp([' '])
    disp(['Line'])
    disp(['1.Compensation of employees,paid...' num2str(nipa_bank.a1(year,1),'%05.1f') ' '  '15.Personal consumption expenditure...' num2str(nipa_bank.a1(year,15),'%05.1f')])
    disp(['2.  Wages and salaries.............' num2str(nipa_bank.a1(year,2),'%05.1f') ' '  '16.  Goods............................' num2str(nipa_bank.a1(year,16),'%05.1f')])
    disp(['3.    Domestic.....................' num2str(nipa_bank.a1(year,3),'%05.1f') ' '  '17.    Durable Goods..................' num2str(nipa_bank.a1(year,17),'%05.1f')])
    disp(['4.    Rest of the world............' num2str(nipa_bank.a1(year,4),'%05.1f') ' '  '18.    Nondurable Goods...............' num2str(nipa_bank.a1(year,18),'%05.1f')])
    disp(['5.  Supplements to wages...........' num2str(nipa_bank.a1(year,5),'%05.1f') ' '  '19.  Services.........................' num2str(nipa_bank.a1(year,19),'%05.1f')])
    disp(['6.Taxes on production and imports..' num2str(nipa_bank.a1(year,6),'%05.1f') ' '  '20.Gross Domestic Investment..........' num2str(nipa_bank.a1(year,20),'%05.1f')])
    disp(['7.Less subsidies...................' num2str(nipa_bank.a1(year,7),'%05.1f') ' '  '21.  Fixed Investment.................' num2str(nipa_bank.a1(year,21),'%05.1f')])
    disp(['8.Net operating surplus............' num2str(nipa_bank.a1(year,8),'%05.1f') ' '  '22.    Nonresidential.................' num2str(nipa_bank.a1(year,22),'%05.1f')])
    disp(['9.  Private enterprise.............' num2str(nipa_bank.a1(year,9),'%05.1f') ' '  '23.      Structures...................' num2str(nipa_bank.a1(year,23),'%05.1f')])
    disp(['10. Government surplus.............' num2str(nipa_bank.a1(year,10),'%05.1f') ' ' '24.      Equipment....................' num2str(nipa_bank.a1(year,24),'%05.1f')])
    disp(['11.Consumption of capital..........' num2str(nipa_bank.a1(year,11),'%05.1f') ' ' '25.      Intellectual property........' num2str(nipa_bank.a1(year,25),'%05.1f')])
    disp(['                                        '                                    ' ' '26.    Residential....................' num2str(nipa_bank.a1(year,26),'%05.1f')])
    disp(['12.Gross domestic income...........' num2str(nipa_bank.a1(year,12),'%05.1f') ' ' '27.  Change in inventories............' num2str(nipa_bank.a1(year,27),'%05.1f')])
    disp(['                                        '                                    ' ' '28.Net exports........................' num2str(nipa_bank.a1(year,28),'%05.1f')])
    disp(['13.Statistical discrepancy.........' num2str(nipa_bank.a1(year,13),'%05.1f') ' ' '29.  Exports..........................' num2str(nipa_bank.a1(year,29),'%05.1f')])
    disp(['                                        '                                    ' ' '30.  Imports..........................' num2str(nipa_bank.a1(year,30),'%05.1f')])
    disp(['                                        '                                    ' ' '31.Government con and inv.............' num2str(nipa_bank.a1(year,31),'%05.1f')])
    disp(['                                        '                                    ' ' '32.  Federal..........................' num2str(nipa_bank.a1(year,32),'%05.1f')])
    disp(['                                        '                                    ' ' '33.    National defence...............' num2str(nipa_bank.a1(year,33),'%05.1f')])
    disp(['                                        '                                    ' ' '34.    Nondefence.....................' num2str(nipa_bank.a1(year,34),'%05.1f')])
    disp(['                                        '                                    ' ' '35.  State and local..................' num2str(nipa_bank.a1(year,35),'%05.1f')])
    disp(['14.Gross domestic product..........' num2str(nipa_bank.a1(year,14),'%05.1f') ' ' '36.Gross domestic product.............' num2str(nipa_bank.a1(year,36),'%05.1f')])
    
    disp([' '])
    disp([' '])
    disp([' '])
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Account 3                      %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % we want to adjust this as % of personal income, thus devide by line
    % 26.
    nipa_bank.a3(year,:) = 100*nipa_bank.a3(year,:)./nipa_bank.a3(year,26);
    
    disp(['Account 3. Personal Income and Outlay Account'])
    disp([' '])
    disp(['For the year ' num2str(year)])
    disp([' '])
    disp(['Values expressed as a % of Personal Income'])
    disp([' '])
    disp(['Line'])
    disp(['1.Personal Current Taxes...........' num2str(nipa_bank.a3(year,1),'%05.1f') ' '  '10.Compensation of employees..........' num2str(nipa_bank.a3(year,10),'%05.1f')])
    disp(['2.Personal Outlays.................' num2str(nipa_bank.a3(year,2),'%05.1f') ' '  '11.  Wages and salaries...............' num2str(nipa_bank.a3(year,11),'%05.1f')])
    disp(['3.  Personal consumption...........' num2str(nipa_bank.a3(year,3),'%05.1f') ' '  '12.    Domestic.......................' num2str(nipa_bank.a3(year,12),'%05.1f')])
    disp(['4.  Personal interest payments.....' num2str(nipa_bank.a3(year,4),'%05.1f') ' '  '13.    Rest of world..................' num2str(nipa_bank.a3(year,13),'%05.1f')])
    disp(['5.  Personal transfer payments.....' num2str(nipa_bank.a3(year,5),'%05.1f') ' '  '14.  Supplements to wages.............' num2str(nipa_bank.a3(year,14),'%05.1f')])
    disp(['6.    To the government............' num2str(nipa_bank.a3(year,6),'%05.1f') ' '  '15.    Employer pension contribution..' num2str(nipa_bank.a3(year,15),'%05.1f')])
    disp(['7.    To the rest of the world.....' num2str(nipa_bank.a3(year,7),'%05.1f') ' '  '16.    Employer gov social insurance..' num2str(nipa_bank.a3(year,16),'%05.1f')])
    disp(['8.Personal Savings.................' num2str(nipa_bank.a3(year,8),'%05.1f') ' '  '17.Proprietors income.................' num2str(nipa_bank.a3(year,17),'%05.1f')])
    disp(['                                        '                                    ' ' '18.Rental income......................' num2str(nipa_bank.a3(year,18),'%05.1f')])
    disp(['                                        '                                    ' ' '19.Personal asset income..............' num2str(nipa_bank.a3(year,19),'%05.1f')])
    disp(['                                        '                                    ' ' '10.  Personal interest income.........' num2str(nipa_bank.a3(year,20),'%05.1f')])
    disp(['                                        '                                    ' ' '21.  Personal dividend income.........' num2str(nipa_bank.a3(year,21),'%05.1f')])
    disp(['                                        '                                    ' ' '22.Personal transfer receipts.........' num2str(nipa_bank.a3(year,22),'%05.1f')])
    disp(['                                        '                                    ' ' '23.  Government benefits..............' num2str(nipa_bank.a3(year,23),'%05.1f')])
    disp(['                                        '                                    ' ' '24.  From business....................' num2str(nipa_bank.a3(year,24),'%05.1f')])
    disp(['                                        '                                    ' ' '25.Less contr for gov soc ins.........' num2str(nipa_bank.a3(year,25),'%05.1f')])
    disp(['9.Personal taxes,outlays,savings...' num2str(nipa_bank.a3(year,9),'%05.1f')  ' ' '26.Personal Income....................' num2str(nipa_bank.a3(year,26),'%05.1f')])
    
    disp([' '])
    disp([' '])
    disp([' '])
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Random Shit                    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    disp(['Account YYY. Random statistics'])
    disp([' '])
    disp(['For the year ' num2str(year)])
    disp([' '])
    disp([' '])
    disp(['Line'])
    disp(['1.Personal Savings Rate............' num2str(nipa_bank.personal_savings_rate(year),'%05.2f') ' '  '9.Bequests / Capital................' num2str(nipa_bank.beq_cap(year),'%05.2f')])
    disp(['2.Capital Output...................' num2str(nipa_bank.KY(year),'%05.2f') ' '                     '10.Bequests / Income................' num2str(nipa_bank.beq_inc(year),'%05.2f')])
    disp(['3.Capital Labor....................' num2str(nipa_bank.KL(year),'%05.2f') ' '                     '11.Personal Debt / Capital...........' num2str(nipa_bank.debt_cap(year),'%05.2f')])
    disp(['4.Capital Labor Adj................' num2str(nipa_bank.KL_adj(year),'%05.2f') ' '                 '12.Personal Debt / Inc..............' num2str(nipa_bank.debt_inc(year),'%05.2f')])
    disp(['5.Income Gini......................' num2str(nipa_bank.income_gini(year),'%05.2f') ' '            '13.Support Ratio....................' num2str(nipa_bank.support_ratio(year),'%05.2f')])
    disp(['6.Wealth Gini......................' num2str(nipa_bank.wealth_gini(year),'%05.2f') ' '            '14.Labor Share......................' num2str(nipa_bank.inc_share(year,1),'%05.2f')])
    disp(['7.Capital Share....................' num2str(nipa_bank.inc_share(year,2),'%05.2f') ' '            '15.Profit Share.....................' num2str(nipa_bank.inc_share(year,3),'%05.2f')])
    disp(['8.Average Return...................' num2str(nipa_bank.average_return(year,1),'%05.2f') ' '            '16.Nope.............................' num2str(0,'%05.2f')])
    
    
    disp([' '])
    disp([' '])
    disp([' '])
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Bequest Account                %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    nipa_bank.beq_account(year,:) = 100 * nipa_bank.beq_account(year,:)  / nipa_bank.beq_account(year,3);
    disp(['Account BEQ. Bequest Account Received / Given LAST Period'])
    disp([' '])
    disp(['For the year ' num2str(year)])
    disp([' '])
    disp(['Values expressed as a % of Total Bequests'])
    disp([' '])
    disp(['Line'])
    disp(['1.Intentional Received.............' num2str(nipa_bank.beq_account(year,1),'%05.1f') ' '                     '4.Intentional Given.................' num2str(nipa_bank.beq_account(year,4),'%05.1f')])
    disp(['2.Unintentional Received...........' num2str(nipa_bank.beq_account(year,2),'%05.1f') ' '                     '5.Unintentional Given...............' num2str(nipa_bank.beq_account(year,5),'%05.1f')])
    disp(['3.Total Received...................' num2str(nipa_bank.beq_account(year,3),'%05.1f') ' '                     '6.Total Given.......................' num2str(nipa_bank.beq_account(year,6),'%05.1f')])

    disp([' '])
    disp([' '])
    disp([' '])
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Annuity                        %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    nipa_bank.ann_account(year,:) = 100 * nipa_bank.ann_account(year,:)  / nipa_bank.ann_account(year,3);
    disp(['Account ANN. Annuity Firm Revenue and Cost Statement'])
    disp([' '])
    disp(['For the year ' num2str(year)])
    disp([' '])
    disp(['Values expressed as a % of Total Annuity Firm Revenues'])
    disp([' '])
    disp(['Line'])
    disp(['1.Lending to individuals rev.......' num2str(nipa_bank.ann_account(year,1),'%05.1f') ' '                     '4.Payments to Individuals...........' num2str(nipa_bank.ann_account(year,4),'%05.1f')])
    disp(['2.Lending to firms rev.............' num2str(nipa_bank.ann_account(year,2),'%05.1f') ' '                     '5.Payments to firms.................' num2str(0,'%05.1f')])
    disp(['3.Total Rev........................' num2str(nipa_bank.ann_account(year,3),'%05.1f') ' '                     '6.Total payments....................' num2str(nipa_bank.ann_account(year,4),'%05.1f')])

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Government Account             %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    nipa_bank.gov_account(year,:) = 100 * nipa_bank.gov_account(year,:) ;% / economy.ag.Y(year);
    disp(['Account GOV: Government Spending and Revenue Account'])
    disp([' '])
    disp(['For the year ' num2str(year)])
    disp([' '])
    disp(['Values expressed as a % of GDP'])
    disp([' '])
    disp(['Line'])
    disp(['1.Gov Spending.....................' num2str(nipa_bank.gov_account(year,1),'%05.1f') ' '                     '4.Tax Revenue.......................' num2str(nipa_bank.gov_account(year,4),'%05.1f')])
    disp(['2.Interest on Debt.................' num2str(nipa_bank.gov_account(year,2),'%05.1f') ' '                     '5.Deficit...........................' num2str(nipa_bank.gov_account(year,5),'%05.1f')])
    disp(['3.Total............................' num2str(nipa_bank.gov_account(year,3),'%05.1f') ' '                     '6.Total.............................' num2str(nipa_bank.gov_account(year,6),'%05.1f')])
    disp(['Gov Debt...........................' num2str(nipa_bank.gov_account(year,7),'%05.1f') ' '                     '6.Blank.............................' num2str(0,'%05.1f')])
    
    output = 1;
end

