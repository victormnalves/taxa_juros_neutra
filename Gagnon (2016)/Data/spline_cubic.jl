# This function computes n cubic splines given the data {x[i],y[i]}, with i=1,...,n+1.
# We assume that {x[i]} are ordered from the smallest to the largest; if not the case, then
# the program sorts the data such that x[i]<=x[i+1]. We also assume that there are no
# repeated value of x[i]; if not the case, then the program computes the average y of repeated 
# values of x and retain only this average as an observation. For n+1 non-repeated values of 
# x, the program thus computes n splines. 
# 
# For each points[j], the program evaluates the i-th spline such that x[i]<=points[j]<x[i+1].
# If points[j] is either smaller than x[1] or larger than x[n+1], the program returns a Taylor-series 
# expansion of order 'offorder' around x[1] or x[n+1], respectively. The default option is to return 
# the values at x[1] and x[n+1].
# 
# The computation of the splines requires two identifying restrictions at x[1] a and x[n+1]. We  
# offer two options through the variable 'method':
# method=0    (default option) Computes 'natural' splines with derivatives at x[1] and x[end] equal to zero. 
# method=1    Computes 'secant-Hermite' splines with first derivative at x[1] equal to slope between x[1] and x[2] 
#             and first derivative at x[end] equal to slope between x[end-1] and x[end]. 
# 
# The output data should always be checked when extrapolating outside of the support of the function 
# to ensure that the extrapolation is well behaved. For a discussion of cubic splines, see
# chapter 6 in Judd (1998), "Numerical Methods in Economics."
#
# Coded for Gagnon, Johannsen, and Lopez-Salido (2016).

# Loading packages
using Statistics

function spline_cubic(x, y, points, method=0, offorder=0) 
    
    # Ordering the data and dealing with repetitions
    x , y = ordernorep(x,y)
    
    # Creating the vector containing the spline coefficients and the matric containing the identifying restrictions
    n = length(x) - 1         # Number of polynomials between data points
    lhs = zeros(4 * n)        # Vector of intercepts of identifying restrictions
    rhs = zeros(4 * n, 4 * n) # Matrix of coefficients in identifying restrictions
    idx = 0
    
    # Right-spline evaluation at x[i] returns y[i]
    for i = 1:n
        idx = idx + 1
        lhs[idx] = y[i+1]
        rhs[idx, i] = 1.0
        rhs[idx, n+i] = x[i+1]
        rhs[idx, 2*n+i] = x[i+1]^2.0
        rhs[idx, 3*n+i] = x[i+1]^3.0
    end
    for i = 1:n  # Left evaluation returns y[i] 
        idx = idx + 1
        lhs[idx] = y[i]
        rhs[idx, i] = 1.0
        rhs[idx, n+i] = x[i]
        rhs[idx, 2*n+i] = x[i]^2.0
        rhs[idx, 3*n+i] = x[i]^3.0
    end
    for i = 1:(n-1) # same first derivative on left and right
        idx = idx + 1
#        lhs[idx] = 0
        rhs[idx, n+i] = 1.0
        rhs[idx, 2*n+i] = 2.0*x[i+1]
        rhs[idx, 3*n+i] = 3.0*x[i+1]^2.0
        rhs[idx, n+i+1] = -1.0
        rhs[idx, 2*n+i+1] = -2.0*x[i+1]
        rhs[idx, 3*n+i+1] = -3.0*x[i+1]^2.0
    end
    for i = 1:(n-1)  # Same second derivative on left and right
        idx = idx + 1
#        lhs[idx] = 0.0
        rhs[idx, 2*n+i] = 2.0
        rhs[idx, 3*n+i] = 6.0*x[i+1]
        rhs[idx, 2*n+i+1] = -2.0
        rhs[idx, 3*n+i+1] = -6.0*x[i+1]
    end
    
    # Extra two identifying restrictions
    idx = idx + 1
    if method==0 # natural splines
        lhs[idx] = 0.0
    elseif method==1  # secant-Hermite
        lhs[idx] = (y[2]-y[1])/(x[2]-x[1])
    end
    rhs[idx, n+1] = 1.0
    rhs[idx, 2*n+1] = 2.0*x[1]
    rhs[idx, 3*n+1] = 3.0*x[1]^2.0
    idx = idx + 1
    if method==0 # natural splines
        lhs[idx] = 0.0
    elseif method==1  # secant-Hermite
        lhs[idx] = (y[n+1]-y[n])/(x[n+1]-x[n])
    end
    rhs[idx, 2*n] = 1.0
    rhs[idx, 3*n] = 2.0*x[n+1]
    rhs[idx, 4*n] = 3.0*x[n+1]^2.0
    
    # Solving linear system of coefficients
    parameters = rhs \ lhs
    a = parameters[1:n]
    b = parameters[(n+1):(2*n)]
    c = parameters[(2*n+1):(3*n)]
    d = parameters[(3*n+1):(4*n)]
    
    # Extrapolating the solution at evaluation points
    vals = collect(points) * 0.0
    for j = 1:length(vals)
        ngt = sum(x .<= points[j])
        if ngt == 0   # points[j] < x[1]
            # Zero-order approximation
            vals[j] = a[1] + b[1]*x[1] + c[1]*x[1]^2.0 + d[1]*x[1]^3.0
            if offorder>=1 # First-order approximation
                vals[j] = vals[j] + (b[1] + 2.0*c[1]*x[1] + 3.0*d[1]*x[1]^2.0)*(points[j]-x[1])
                if offorder>=2 # Second-order approximation
                    vals[j] = vals[j] + (c[1] + 3.0*d[1]*x[1])*(points[j]-x[1])^2.0
                    if offorder>=3 # Third-order approximation
                        vals[j] = vals[j] + d[1]*x[1]*(points[j]-x[1])^3.0
                    end
                end
            end
        elseif ngt == (n+1) # points[j] > x[end]
            vals[j] = a[n] + b[n]*x[n+1] + c[n]*x[n+1]^2.0 + d[n]*x[n+1]^3.0
            if offorder>=1
                vals[j] = vals[j] + (b[n] + 2.0*c[n]*x[n+1] + 3.0*d[n]*x[n+1]^2.0)*(points[j]-x[n+1])
                if offorder>=2
                    vals[j]=vals[j] + (c[n] + 3.0*d[n]*x[n+1])*(points[j]-x[n+1])^2.0
                    if offorder>=3
                        vals[j]=vals[j] + d[n]*(points[j]-x[n+1])^3.0
                    end
                end
            end
        else # x[1] <= points[j] <= x[end]
            vals[j] = a[ngt] + b[ngt] * points[j] + c[ngt] * points[j]^2.0 + d[ngt] * points[j]^3.0
        end
    end
    
    return vals
end
