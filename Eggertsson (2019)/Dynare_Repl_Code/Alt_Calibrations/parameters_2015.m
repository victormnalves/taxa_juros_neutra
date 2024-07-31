function [param] = parameters_2015

global run_alt_number

if run_alt_number==1
    
    n=-0.002578211534468;
    fert_26=0.9375;
    e=1.009368277002316;
    b=1.172880806665463;
    D=0.345045166791504;
    theta=4.623205788121817;
    AL_growth=0.006460448230438;
    AL=2.99049071248668;
    AK=1;

    param = [n, fert_26, e, b, D, theta, AL_growth, AL, AK];

elseif run_alt_number==2
    
    n=-0.002578211534468;
    fert_26=0.9375;
    e=1.009368277002316;
    b=1.172880806665463;
    D=0.233418475865886;
    theta=4.898957630032689;
    AL_growth=0.006460448230438;
    AL=2.99049071248668;
    AK=1;

    param = [n, fert_26, e, b, D, theta, AL_growth, AL, AK];
    
elseif run_alt_number==3
    
    n=-0.002578211534468;
    fert_26=0.9375;
    e=1.009368277002316;
    b=1.172880806665463;
    D=0.278544836590899;
    theta=5.480271248212643;
    AL_growth=0.006460448230438;
    AL=2.99049071248668;
    AK=1;

    param = [n, fert_26, e, b, D, theta, AL_growth, AL, AK];
        
elseif run_alt_number==4
    
    n=-0.002578211534468;
    fert_26=0.9375;
    e=1.009368277002316;
    b=1.172880806665463;
    D=0.261173987577748;
    theta=4.917840559496626;
    AL_growth=0.006460448230438;
    AL=2.99049071248668;
    AK=1;

    param = [n, fert_26, e, b, D, theta, AL_growth, AL, AK];
            
end

