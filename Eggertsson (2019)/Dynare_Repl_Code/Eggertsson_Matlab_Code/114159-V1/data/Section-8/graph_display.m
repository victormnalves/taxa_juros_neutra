function [ output ] = graph_display(nipa_bank,ps,prod,gov,policy,economy,prices,run_schedule,aux,graph_name,comparison)

        % Aux is a set of auxiliary data which may be useful
        % graph_name is going to be the name of the graphs :)

        output = 1;
        
        csv_output = [];
                
        year_vec = [1:152]' + 1969;
        
        % Interest Rate Graph

        hold off
        
        figure
        plot(year_vec,100*prices.r,'LineWidth',2)
        legend('Interest Rates (%)','Location','northoutside')
        xlabel('Period')
        ylabel('Rates')
        title('Interest Rates')
        print(strcat('figures/',graph_name,'_transition','.png'),'-dpng') 
        
        table_output = table([1:1:policy.nt+2]','VariableNames',{'Year'});
        
        table_output = [table_output,table(prices.r,'VariableNames',{'interest_rates'})];
        
        csv_output = [csv_output,prices.r];

        % Plot Average Return
        
        figure
        plot(year_vec,nipa_bank.average_return,'LineWidth',2)
        legend('Average Return (%)','Location','northoutside')
        xlabel('Period')
        ylabel('Rates')
        title('Average Return')
        print(strcat('figures/',graph_name,'_averageret','.png'),'-dpng') 
        
        table_output = [table_output,table(nipa_bank.average_return,'VariableNames',{'average_return'})];
        
        % Add MPK and RETK to the output of the model

        table_output = [table_output,table(prices.mpk,'VariableNames',{'mpk'})];
        
        table_output = [table_output,table(prices.rentk,'VariableNames',{'rental_capital'})];
        
        % Set 1 of Figures
            % Output per Person
            % Output per Worker
            % Productivity Growth
            
        csv_output = [csv_output,nipa_bank.Y_pp./comparison.Y_pp,nipa_bank.Y_pw./comparison.Y_pw,economy.ag.AL_growth];

        table_output = [table_output,table(nipa_bank.Y_pp./comparison.Y_pp,'VariableNames',{'output_per_person'})];
        table_output = [table_output,table(nipa_bank.Y_pw./comparison.Y_pw,'VariableNames',{'output_per_worker'})];
        table_output = [table_output,table(economy.ag.AL_growth,'VariableNames',{'productivity_growth'})];
        
        figure
        
        subplot(3,1,1) 
        plot(year_vec,nipa_bank.Y_pp./comparison.Y_pp,'LineWidth',2)
        xlabel('Period')
        ylabel('GDP per Person')
        title('Output per Person Transition')
                
        subplot(3,1,2) 
        plot(year_vec,nipa_bank.Y_pw./comparison.Y_pw,'LineWidth',2)
        xlabel('Period')
        ylabel('GDP per Worker')
        title('Output per Worker Transition')
        
        subplot(3,1,3) 
        plot(year_vec,economy.ag.AL_growth,'LineWidth',2)
        xlabel('Period')
        ylabel('Prod Growth')
        title('Labor Productivity Growth Rate')
                
        print(strcat('figures/',graph_name,'_set1','.png'),'-dpng') 
        
        
        % Set 2 of figures
            % KY Ratio
            % IY Ratio
            % Nominal I/Y Ratio
            % Personal Savings Rate
        
        csv_output = [csv_output,nipa_bank.KY,nipa_bank.IY,nipa_bank.IY.*economy.relp,nipa_bank.personal_savings_rate ];

        table_output = [table_output,table(nipa_bank.KY,'VariableNames',{'capital_output'})];
        table_output = [table_output,table(nipa_bank.IY,'VariableNames',{'real_investment_output'})];
        table_output = [table_output,table(nipa_bank.IY.*economy.relp,'VariableNames',{'nom_investment_output'})];
        table_output = [table_output,table(nipa_bank.personal_savings_rate,'VariableNames',{'personal_savings'})];
        
        figure
        
        subplot(4,1,1) 
        plot(year_vec,nipa_bank.KY,'LineWidth',2)
        axis([min(year_vec) max(year_vec) 1.5 5])
        axis 'auto x'
        xlabel('Period')
        ylabel('K / Y')
        title('Capital Output Ratio')
        
        subplot(4,1,2) 
        plot(year_vec,nipa_bank.IY,'LineWidth',2)
        axis([min(year_vec) max(year_vec) .1 .5])
        axis 'auto x'
        xlabel('Period')
        ylabel('I / Y')
        title('Investment to Output Ratio')
        
        subplot(4,1,3) 
        plot(year_vec,nipa_bank.IY.*economy.relp,'LineWidth',2)
        axis([min(year_vec) max(year_vec) .1 .5])
        axis 'auto x'
        xlabel('Period')
        ylabel('Nominal I / Y')
        title('Nominal Investment to Output Ratio')
        
        subplot(4,1,4) 
        plot(year_vec,nipa_bank.personal_savings_rate,'LineWidth',2)
        axis([min(year_vec) max(year_vec) -.05 .3])
        axis 'auto x'
        xlabel('Period')
        ylabel('Personal Savings')
        title('Personal Savings Rate')
        
        print(strcat('figures/',graph_name,'_set2','.png'),'-dpng') 
        
        % Set 3 of figures
            % Non Interest government spending / GDP
            % Total Gov spending / GDP
            % Tax revenue to GDP
            % Deficit to GDP
            % Debt to GDP
        
        csv_output = [csv_output,nipa_bank.gov_account(:,1),nipa_bank.gov_account(:,3),nipa_bank.gov_account(:,4),nipa_bank.gov_account(:,5),nipa_bank.gov_account(:,7)];

        table_output = [table_output,table(nipa_bank.gov_account(:,1),'VariableNames',{'non_interest_gov_gdp'})];
        table_output = [table_output,table(nipa_bank.gov_account(:,3),'VariableNames',{'total_gov_gdp'})];
        table_output = [table_output,table(nipa_bank.gov_account(:,4),'VariableNames',{'tax_rev_gdp'})];
        table_output = [table_output,table(nipa_bank.gov_account(:,5),'VariableNames',{'deficit_gdp'})];
        table_output = [table_output,table(nipa_bank.gov_account(:,7),'VariableNames',{'debt_gdp'})];
        
        
        figure
        
        subplot(3,1,1) 
        plot(year_vec,nipa_bank.gov_account(:,1),'LineWidth',2)
        axis([min(year_vec) max(year_vec) .1 .3])
        axis 'auto x'
        xlabel('Period')
        ylabel('Gov Spending / GDP')
        title('Non Interest Government Spending to Output Ratio')
        
        subplot(3,1,2) 
        plot(year_vec,nipa_bank.gov_account(:,3),'LineWidth',2)
        axis([min(year_vec) max(year_vec) .1 .3])
        axis 'auto x'
        xlabel('Period')
        ylabel('Gov Spending / GDP')
        title('Total Government Spending to Output Ratio')
        
        subplot(3,1,3) 
        plot(year_vec,economy.nfa_gdp,'LineWidth',2)
        axis([min(year_vec) max(year_vec) -.05 .5])
        axis 'auto x'
        xlabel('Period')
        ylabel('Net Foreign Asset / GDP')
        title('Net Foreign Asset Ownership to Output Ratio')
        
        print(strcat('figures/',graph_name,'_set3','.png'),'-dpng') 
        
        % Set 4 of figures
        
        
        figure
        
        subplot(3,1,1) 
        plot(year_vec,nipa_bank.gov_account(:,4),'LineWidth',2)
        axis([min(year_vec) max(year_vec) .1 .3])
        axis 'auto x'
        xlabel('Period')
        ylabel('Tax / GDP')
        title('Tax Revenue to Output')
        
        subplot(3,1,2) 
        plot(year_vec,nipa_bank.gov_account(:,5),'LineWidth',2)
        axis([min(year_vec) max(year_vec) 0 .1])
        axis 'auto x'
        xlabel('Period')
        ylabel('Deficit')
        title('Deficit to Output')
        
        subplot(3,1,3) 
        plot(year_vec,nipa_bank.gov_account(:,7),'LineWidth',2)
        axis([min(year_vec) max(year_vec) .1 1.5])
        axis 'auto x'
        xlabel('Period')
        ylabel('Debt / GDP')
        title('Debt to Output Ratio')
                
        print(strcat('figures/',graph_name,'_set4','.png'),'-dpng') 
        
        
        % Set 5 of figures
        
            % Fertility
            % Population Growth
            % Support Ratio
        csv_output = [csv_output,ps(1).demog.num_kids,nipa_bank.pop_growth,nipa_bank.support_ratio];

        table_output = [table_output,table(ps(1).demog.num_kids,'VariableNames',{'total_fertility'})];
        table_output = [table_output,table(nipa_bank.pop_growth,'VariableNames',{'pop_growth'})];
        table_output = [table_output,table(nipa_bank.support_ratio,'VariableNames',{'support_ratio'})];

        figure
        
        subplot(3,1,1) 
        plot(year_vec,ps(1).demog.num_kids,'LineWidth',2)
        axis([min(year_vec) max(year_vec) .7 2])
        axis 'auto x'
        xlabel('Period')
        ylabel('TFR')
        title('Total Fertility Rate')
        
        subplot(3,1,2) 
        plot(year_vec,nipa_bank.pop_growth,'LineWidth',2)
        axis([min(year_vec) max(year_vec) -.02 .05])
        axis 'auto x'
        xlabel('Period')
        ylabel('Popgrowth')
        title('Population Growth Rate')
        
        subplot(3,1,3) 
        plot(year_vec,nipa_bank.support_ratio,'LineWidth',2)
        axis([min(year_vec) max(year_vec) .5 1])
        axis 'auto x'
        xlabel('Period')
        ylabel('Support Ratio')
        title('Support Ratio')
                
        print(strcat('figures/',graph_name,'_set5','.png'),'-dpng') 
       
        % Set 6 of figures
        	% Average age
            % Percent 1-20
            % Percent 21-40
            % Percent 41-56
        csv_output = [csv_output,nipa_bank.average_age,nipa_bank.percent_120,nipa_bank.percent_2140,nipa_bank.percent_4156];
        
        table_output = [table_output,table(nipa_bank.average_age,'VariableNames',{'average_age'})];
        table_output = [table_output,table(nipa_bank.percent_120,'VariableNames',{'pct_120'})];
        table_output = [table_output,table(nipa_bank.percent_2140,'VariableNames',{'pct_2140'})];
        table_output = [table_output,table(nipa_bank.percent_4156,'VariableNames',{'pct_4156'})];
        
        figure
        
        subplot(4,1,1) 
        plot(year_vec,nipa_bank.average_age,'LineWidth',2)
        axis([min(year_vec) max(year_vec) 20 30])
        axis 'auto x'
        xlabel('Period')
        ylabel('Average Age')
        title('Average Age')
        
        subplot(4,1,2) 
        plot(year_vec,nipa_bank.percent_120,'LineWidth',2)
        axis([min(year_vec) max(year_vec) .3 .6])
        axis 'auto x'
        xlabel('Period')
        ylabel('% Age 1-20')
        title('Percent Ages 1-20')
        
        subplot(4,1,3) 
        plot(year_vec,nipa_bank.percent_2140,'LineWidth',2)
        axis([min(year_vec) max(year_vec) .3 .6])
        axis 'auto x'
        xlabel('Period')
        ylabel('% Age 21-40')
        title('Percent Ages 21-40')
        
        subplot(4,1,4) 
        plot(year_vec,nipa_bank.percent_4156,'LineWidth',2)
        axis([min(year_vec) max(year_vec) .1 .3])
        axis 'auto x'
        xlabel('Period')
        ylabel('% Age 41-56')
        title('Percent Ages 41-56')
                
        print(strcat('figures/',graph_name,'_set6','.png'),'-dpng') 
       
        % Now Create Plots of the Demographic Pyramid
        
        figure
        
        subplot(1,4,1)
        barh(26:81,ps.demog.pop(1,:)/economy.ag.pop(1),'r')
        xlabel('% of Total Population')
        ylabel('Age')
        title('Dem Pyramid Year 1970')
        
        subplot(1,4,2)
        barh(26:81,ps.demog.pop(25,:)/economy.ag.pop(25),'b')
        xlabel('% of Total Population')
        ylabel('Age')
        title('Dem Pyramid Year 1995')
        
        subplot(1,4,3)
        barh(26:81,ps.demog.pop(50,:)/economy.ag.pop(50),'g')
        xlabel('% of Total Population')
        ylabel('Age')
        title('Dem Pyramid Year 2020')
        
        subplot(1,4,4)
        barh(26:81,ps.demog.pop(75,:)/economy.ag.pop(75),'y')
        xlabel('% of Total Population')
        ylabel('Age')
        title('Dem Pyramid Year 245')
        
        print(strcat('figures/',graph_name,'_set7','.png'),'-dpng') 
        
       % Now Create Plots of the Demographic Pyramid
        
        figure
        
        subplot(1,3,1)
        barh(26:81,ps.demog.pop(100,:)/economy.ag.pop(100),'r')
        xlabel('% of Total Population')
        ylabel('Age')
        title('Dem Pyramid Year 100')
        
        subplot(1,3,2)
        barh(26:81,ps.demog.pop(125,:)/economy.ag.pop(125),'b')
        xlabel('% of Total Population')
        ylabel('Age')
        title('Dem Pyramid Year 125')
        
        subplot(1,3,3)
        barh(26:81,ps.demog.pop(152,:)/economy.ag.pop(152),'g')
        xlabel('% of Total Population')
        ylabel('Age')
        title('Dem Pyramid Year 150')
                
        print(strcat('figures/',graph_name,'_set8','.png'),'-dpng') 
        
        % Set 9 of figures
        
            % Labor Share
            % Capital Share
            % Profit Share

        table_output = [table_output,table(nipa_bank.inc_share(:,1),'VariableNames',{'labor_share'})];
        table_output = [table_output,table(nipa_bank.inc_share(:,2),'VariableNames',{'capital_share'})];
        table_output = [table_output,table(nipa_bank.inc_share(:,3),'VariableNames',{'profit_share'})];

        figure
        
        subplot(3,1,1) 
        plot(year_vec,nipa_bank.inc_share(:,1),'LineWidth',2)
        axis([min(year_vec) max(year_vec) 0 100])
        axis 'auto x'
        xlabel('Period')
        ylabel('Lab Sh')
        title('Labor Share')
        
        subplot(3,1,2) 
        plot(year_vec,nipa_bank.inc_share(:,2),'LineWidth',2)
        axis([min(year_vec) max(year_vec) 0 100])
        axis 'auto x'
        xlabel('Period')
        ylabel('Cap Sh')
        title('Capital Share')
        
        subplot(3,1,3) 
        plot(year_vec,nipa_bank.inc_share(:,3),'LineWidth',2)
        axis([min(year_vec) max(year_vec) 0 100])
        axis 'auto x'
        xlabel('Period')
        ylabel('Profit Sh')
        title('Profit Share')
                
        print(strcat('figures/',graph_name,'_set9','.png'),'-dpng') 
        
        
        % Create new population pyramid, including the population sizes of
        % indidiuals age 1 to 25
        
            pop_matrix = zeros(size(ps.demog.pop,1),size(ps.demog.pop,2)+25);
        
            % Fill in population 26 and above
            pop_matrix(:,26:81) = ps.demog.pop;
            
            % Fill in Population Age 1 to 25
            for year=1:(policy.nt+2)
                                
                for age=1:25
                                        
                    % Indivdiuals who are age `age' in year `year' will age 26 in
                    % year 'year + 26 - 'age'' 
                    
                    year_index = year + 26 - age;
                    
                    if year_index <=policy.nt+2
                        
                        pop_matrix(year,age) = pop_matrix(year_index,26);
                    
                    end
                    
                end

            end
            
            pop_matrix_ag = sum(pop_matrix,2);
            
            % Now Create New Demographic Pyramid

            
            figure
        
            subplot(1,4,1)
            barh(1:81,pop_matrix(1,:)/pop_matrix_ag(1),'r')
            hold on
            plot([0;.04],[26;26],'c','LineWidth',1.5)
            xlabel('% of Total Population')
            ylabel('Age')
            title('Dem Pyramid Year 1970')
            hold off
            
            subplot(1,4,2)
            barh(1:81,pop_matrix(25,:)/pop_matrix_ag(25),'b')
            hold on
            plot([0;.04],[26;26],'c','LineWidth',1.5)
            hold off
            xlabel('% of Total Population')
            ylabel('Age')
            title('Dem Pyramid Year 1995')

            subplot(1,4,3)
            barh(1:81,pop_matrix(50,:)/pop_matrix_ag(50),'g')
            hold on
            plot([0;.04],[26;26],'c','LineWidth',1.5)
            hold off
            xlabel('% of Total Population')
            ylabel('Age')
            title('Dem Pyramid Year 2020')

            subplot(1,4,4)
            barh(1:81,pop_matrix(75,:)/pop_matrix_ag(75),'y')
            hold on
            plot([0;.04],[26;26],'c','LineWidth',1.5)
            hold off
            xlabel('% of Total Population')
            ylabel('Age')
            title('Dem Pyramid Year 2045')

            print(strcat('figures/',graph_name,'_new_pop','.png'),'-dpng') 

        
        % Output to CSV File
        
        csvwrite(strcat('figures/',graph_name,'_figures','.csv'),csv_output)
        
        writetable(table_output,strcat('figures/',graph_name,'_tables','.csv'))
        
        close all
         

end