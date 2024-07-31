function [ util] = test_crrautil(C,beta,gamma)

% Basic Log Utility Function, with discount rate beta

    lifespan = length(C);

    util = 0;

    for t=1:lifespan

        util = util + beta^(t-1)*C(t)^(1-1/gamma);

    end
    
    util = util*1/(1-1/gamma);

    util = -util;

end

