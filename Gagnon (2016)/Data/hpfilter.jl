# This file implements a Hodrick-Prescott filter. The inputs are:
# y  Vector of data to be filtered
# w  Lagrange multiplier


# Loading package
using LinearAlgebra

function hpfilter(y,w)
    # Making sure the unfiltered data are in vector form
    if size(y,1)<size(y,2)
       y=y';
    end

    t=size(y,1);
    a=6*w+1;
    b=-4*w;
    c=w;
    d=hcat(c,b,a);
    d=ones(t,1)*d;
    m=Diagonal(d[:,3]) + vcat(hcat(zeros(t-1),Diagonal(d[1:t-1,2])),zeros(t)')+ vcat(zeros(t)',hcat(Diagonal(d[1:t-1,2]),zeros(t-1)));
    m=m+vcat(hcat(zeros(t-2),zeros(t-2),Diagonal(d[1:t-2,1])),zeros(t)',zeros(t)') + vcat(zeros(t)',zeros(t)',hcat(Diagonal(d[1:t-2,1]),zeros(t-2),zeros(t-2)));


    m[1,1]=1+w;       m[1,2]=-2*w;
    m[2,1]=-2*w;      m[2,2]=5*w+1;
    m[t-1,t-1]=5*w+1; m[t-1,t]=-2*w;
    m[t,t-1]=-2*w;    m[t,t]=1+w;

    out = inv(Matrix(m))*y;

    return out
end
    
