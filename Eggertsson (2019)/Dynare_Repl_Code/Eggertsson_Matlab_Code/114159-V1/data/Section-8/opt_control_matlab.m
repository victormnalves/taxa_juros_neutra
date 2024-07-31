function [ ps ] = opt_control_matlab(ps,gov,prod,policy,prices,economy,run_schedule,tm)

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop through Types     %
%%%%%%%%%%%%%%%%%%%%%%%%%%


% [C_opt, A_opt] = recursive_solve(ps,gov,prod,policy,prices,economy,run_schedule,1,1,1,1,1,ps.demog.lifespan,0);

for type = 1:economy.num_types
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % Loop through Years     %
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % This indexes people born during this year,
    % or during the transition, people whose optimizaiton begins during
    % this year
    
    for initial_year=tm.year_min:tm.year_max
        
       
        % Else, we only do for people born at the beginning of the year

        tm.initial_age = 1;
        tm.initial_year = initial_year;
        
        % Get Initial Guess :)
        
            x0 = zeros(ps(type).demog.lifespan + 1,1);

            for age= tm.initial_age:ps(type).demog.lifespan

                tm.age = age;
                [year_index,year_born] = create_index(tm,policy);

                x0(age) = prod_adjustment(tm,age,economy,policy)*ps(type).opt.C(year_index,age);

                % Guess for Bequests Given
                if age==ps(type).util.a_bg(year_born)
                    
                    % Also need to adjust bequests for productivity
                    
                    x0(ps(type).demog.lifespan + 1) = prod_adjustment(tm,ps(type).util.a_bg(year_born),economy,policy)*ps(type).opt.bgo(year_index,ps(type).util.a_bg(year_born));
                    
                end
                
            end

            % If bequests are zero, set it to a small level so we can run
            % things.  
            if x0(ps(type).demog.lifespan + 1)==0
                
                x0(ps(type).demog.lifespan + 1) = .0001;
                
            end
            
            % On the first iteration, we guess ones. 
            if x0(1) ==0
                
                x0 = ones(ps(type).demog.lifespan + 1,1);
                
            end
            
            % x0 = ones(ps(type).demog.lifespan + 1,1);
            
        % Now Optimize
        
        options = optimset('MaxFunEvals',10000);
        
        %crrautil(x0,ps(type))
        %test_crrautil(x0(1:56),1/(1+ps.util.delta(1)),ps.util.gamma(1))
        %[cn,ceqn] = cb_nonl_cons(x0,ps(type),prices,economy,policy,tm)
        %[c,ceq] = testing_nl(x0(1:56),ps.demog.hc(1,:),prices.r(1),prices.r(1),ps.opt.dlim(1))
        
        C = fmincon(@(CB) crrautil(CB,ps(type),tm,policy),x0,-eye(ps(type).demog.lifespan+1),zeros(ps(type).demog.lifespan + 1,1),[],[],[],[],@(CB) cb_nonl_cons(CB,ps(type),prices,economy,policy,tm),options); 
        
        [~,~,a] = cb_nonl_cons(C,ps(type),prices,economy,policy,tm);
        % C_alt = fmincon(@(CB) test_crrautil(CB,1/(1+ps.util.delta(1)),ps.util.gamma(1)),x0(1:56),-eye(ps(type).demog.lifespan),zeros(ps(type).demog.lifespan,1),[],[],[],[],@(CB) testing_nl(CB,ps.demog.hc(1,:),prices.r(1),prices.r(1),ps.opt.dlim(1)),options); 
        
        % Recursive Variables
        rec.final_age = ps(type).demog.lifespan;
        rec.final_debt = 0;
        
        ps(type) = opt_lb_master(ps(type),gov,prod,policy,prices,economy,run_schedule,tm,rec,C(1:end-1),a);

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


% 
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
%             ps(type).opt.bgt(year) = ps(type).opt.bgt(year) + ps(type).opt.a(yp1_index,age)*ps(type).demog.pop(year,age-1)*(1-ps(type).demog.ms(year,age-1))*prod_adj;
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
%         elseif year==(policy.nt+2)
% 
%             pop_adj = (1+economy.n(policy.nt+2));
%             ym1_index = policy.nt + 2;
% 
%         end
%         
%         for age=2:ps(type).demog.lifespan
%             
%             ps(type).opt.brt(year) = ps(type).opt.brt(year) + ps(type).opt.a(year,age)*ps(type).demog.pop(ym1_index,age-1)*(1-ps(type).demog.ms(ym1_index,age-1))*(1/pop_adj);
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

