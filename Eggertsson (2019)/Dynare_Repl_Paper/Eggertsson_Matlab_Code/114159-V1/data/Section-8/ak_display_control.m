function [ ret] = ak_display_control(nipa_bank,ps,prod,gov,policy,economy,prices,run_schedule,adj)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ak_display_control
%
% Displays output as in AK 1987 Models

for year=1:21
    
    temp = ak_display(nipa_bank,ps,prod,gov,policy,economy,prices,run_schedule,year,adj);
    
end

for decade=3:15
    
    temp = ak_display(nipa_bank,ps,prod,gov,policy,economy,prices,run_schedule,decade*10,adj);
    
end

ret = 1;
end

