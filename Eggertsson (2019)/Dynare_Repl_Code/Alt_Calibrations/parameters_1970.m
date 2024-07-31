function [param] = parameters_1970

global run_alt_number

if run_alt_number==1

    n=0.013549868016229;
    fert_26=1.4;
    e=1.3;
    b=0.422844989041041;
    D=0.236636866976686;
    theta=7.563866678744176; 
    AL_growth=0.020209202089143;
    AL=1;
    AK=1;

    param = [n, fert_26, e, b, D, theta, AL_growth, AL, AK];

elseif run_alt_number==2
    
    n=0.013549868016229;
    fert_26=1.4;
    e=1.3;
    b=0.422844989041041;
    D=0.154842731526628;
    theta=7.651468691003310; 
    AL_growth=0.020209202089143;
    AL=1;
    AK=1;

    param = [n, fert_26, e, b, D, theta, AL_growth, AL, AK];
    
elseif run_alt_number==3
    
    n=0.013549868016229;
    fert_26=1.4;
    e=1.3;
    b=0.422844989041041;
    D=0.151570145996862;
    theta=7.616852619513112; 
    AL_growth=0.020209202089143;
    AL=1;
    AK=1;

    param = [n, fert_26, e, b, D, theta, AL_growth, AL, AK];
        
elseif run_alt_number==4
    
    n=0.013549868016229;
    fert_26=1.4;
    e=1.3;
    b=0.422844989041041;
    D=0.145580300664455;
    theta=8.747100034823779; 
    AL_growth=0.020209202089143;
    AL=1;
    AK=1;

    param = [n, fert_26, e, b, D, theta, AL_growth, AL, AK];
            
end

