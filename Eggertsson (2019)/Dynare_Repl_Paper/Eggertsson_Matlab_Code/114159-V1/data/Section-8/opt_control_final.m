function [ ps ] = opt_control_final(ps,gov,prod,policy,prices,economy,run_schedule,tm)

% ISS optimization works by calculating the optimal quantities of
% individuals who are born at period 1. All information about the economy
% in period 1 is stored in row 1 of the vectors. How do we get information
% about what the wage will be after period 1, which will be necessary to
% continue the optimization? Since we are proposing an initial steady
% state, we can simply update the wages by the rate of productivity growth.
% This is done with the prod_adjustment program. 

% Now we have the optimal consumption / assets for an individual born in
% period 1. But in order to calculate the equilibrium for the ISS, we need
% information on all of the generations alive in period 1. Similarly, we
% need to fill in our matrices which store the optimal decisions of all
% individuals alive at age 1. Once again, since we are in an ISS, we can
% simply take the optimal consumption decisions of individuals born in
% period 1 and update. 

% This principle is the most complicated with dealing with bequests. 


%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop through Types     %
%%%%%%%%%%%%%%%%%%%%%%%%%%

for type = 1:economy.num_types
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Loop through Years     %
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % This indexes people born during this year,
    % or during the transition, people whose optimizaiton begins during
    % this year
    
    for initial_year=tm.year_min:tm.year_max
        
       
        % For the first year of the transition, we loop through all ages
        % allive
        
        if initial_year==2
            
            for initial_age=1:ps(type).demog.lifespan
                        
                tm.initial_age = initial_age;
                tm.initial_year = initial_year;
                
                
                % Recursive Variables
                rec.final_age = ps(type).demog.lifespan;
                rec.final_debt = 0;
                
                % Finds the optimal consumption
                [C,a] = recursive_solve(ps(type),gov,prod,policy,prices,economy,run_schedule,tm,rec);   
                
                % Stores and compiles the information
                ps(type) = opt_lb_master(ps(type),gov,prod,policy,prices,economy,run_schedule,tm,rec,C,a);
                          
            end
        
        else % Else, we only do for people born at the beginning of the year
                   
                tm.initial_age = 1;
                tm.initial_year = initial_year;
                
                % Recursive Variables
                rec.final_age = ps(type).demog.lifespan;
                rec.final_debt = 0;
                
                % Finds the optimal consumption
                
                if tm.initial_year==130
                    
                    eagle = 55;
                end
                [C,a] = recursive_solve(ps(type),gov,prod,policy,prices,economy,run_schedule,tm,rec);   
                
                % Stores and compiles the information
                ps(type) = opt_lb_master(ps(type),gov,prod,policy,prices,economy,run_schedule,tm,rec,C,a);
                          
            %[C,a] = recursive_solve(ps(type),gov,prod,policy,prices,economy,run_schedule,initial_year,1,year_min,year_max,ps(type).demog.lifespan,0);

            %ps(type) = opt_lb_master(ps(type),gov,prod,policy,prices,economy,run_schedule,initial_year,1,year_min,year_max,C,a);
                       
        end
        
    end % Year
        
end % Type




%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now Fill in Bequests   %
%%%%%%%%%%%%%%%%%%%%%%%%%%

% First, we add up the bequests that people meant to give
% Then, we add in "accidental" bequests from assets

% During the first and last years, we may have to adjust for productivity
% and for population

[ ps ] = bequest_control(ps,gov,prod,policy,prices,economy,run_schedule,tm);

% for type = 1:economy.num_types
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%
%     % Loop through Years     %
%     %%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%     
%     for year=tm.year_min:tm.year_max
%         
%         % Total Up Bequests That are meant to be given
%         
%         ps(type).opt.bgt(year) = ps(type).opt.bg_adj(year,ps(type).util.a_bg(1));
%                 
%         % Now add up the people who died
%         % The bequests given this period will be the initial assets next
%         % period, thus need to get yp1_index
%         % In addition, in period 152 will have to adjust for productivity
%         
%         yp1_index = year + 1;
%         prod_adj = 1;
% 
%         if year==1
% 
%             yp1_index = 1;
%             prod_adj = (1+economy.ag.AL_growth_iss)*(1+economy.ag.A_growth_adj_iss);
% 
%         elseif year==(policy.nt+2)
% 
%             yp1_index = policy.nt+2;
%             prod_adj = (1+economy.ag.AL_growth_fss)*(1+economy.ag.A_growth_adj_fss);
% 
%         end
% 
%         for age=2:ps(type).demog.lifespan
%             
%             % I think the age indexes are correct -- be careful about
%             % changing them. Essentially, the loop is for people's ages
%             % next year, thus the pop is their age-1, but this year. 
%             ps(type).opt.bgt(year) = ps(type).opt.bgt(year) + (1-ps(type).util.annuity_participation(year,age))*ps(type).opt.a(yp1_index,age)*ps(type).demog.pop(year,age-1)*(1-ps(type).demog.ms(year,age-1))*prod_adj;
%         
%         end
%         
%         
%         % Now bequests received total
%         ps(type).opt.brt(year) = ps(type).opt.br_adj(year,ps(type).util.a_br(1));
%                 
%         ym1_index = year - 1;
%         pop_adj = 1;
%         
%         if year==1
% 
%             pop_adj = (1+economy.n(1));
%             ym1_index = 1;
% 
%         % More Warning!
%         elseif year==(policy.nt+2)
% 
%             pop_adj = (1+economy.n(policy.nt+2));
%             ym1_index = policy.nt + 2;
% 
%         end
%         
%         for age=2:ps(type).demog.lifespan
%             
%             ps(type).opt.brt(year) = ps(type).opt.brt(year) + (1-ps(type).util.annuity_participation(ym1_index,age-1))*ps(type).opt.a(year,age)*ps(type).demog.pop(ym1_index,age-1)*(1-ps(type).demog.ms(ym1_index,age-1))*(1/pop_adj);
%         
%         end
%         
%         % Bequests received per person, we simply devide by the population
%         % that is alive
%         
%         ps(type).opt.br(year,ps(type).util.a_br(1)) = ps(type).opt.brt(year) / ps(type).demog.pop(year,ps(type).util.a_br(1));
%         
%     end
%     
% end

end

