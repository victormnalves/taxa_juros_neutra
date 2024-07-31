function [ util] = logutil(C,beta)

% Basic Log Utility Function, with discount rate beta

    lifespan = length(C);

    util = 0;

    for t=1:lifespan

        util = util + beta^(t-1)*log(C(t));

    end

    util = -util;

end

