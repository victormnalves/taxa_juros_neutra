# Function to prepare input data matrices for expected demographics under backward-looking expectations
function matrix_prep(i::Int64, hh::Int64, mat::Array{Float64, 2})
    """
       i is the iteration of the for loop this function will be called in,
       it sets which column is the last one selected of the input data

       hh is the number of columns backwards which we want to select

       mat is the input matrix

       When i<=hh, we don't have sufficient historical data to compute expectations.
       In that case, we fill in demographic data before the initial period using
       data from the initial period.

    """
    out = ones(size(mat)) # always return a matrix of the same size as input
    if(size(out)[2] == 1) # Vector case
        if(i < hh)
            temp = ((hh-i)/hh).*mat[1] + (1/hh)*sum(mat[1:i])
        else
            temp = mean(mat[(i-hh+1):i])
        end
    else # two-dimensional case (matrix, with each column corresponding to a period)
        if (i < hh)
            # We have to duplicate the first column the requisite number of times
            # for example, if i = 2 but hh = 4 we will need to get the average of
            # three times the values in column 1 and one time the values in column 2
#            temp = ((hh - i)/hh).*mat[:, 1] + (1/hh).*mapslices(sum, mat[:, 1:i], 2)
            temp = ((hh - i)/hh).*mat[:, 1] + (1/hh).*mapslices(sum, mat[:, 1:i], dims=2)   # EG: 20200824
        else
            # In the case where i >= hh we just take a straight average
            # over the columns (i-hh+1:i)
#            temp = (1/hh).*mapslices(sum, mat[:,(i-hh+1):i], 2)
            temp = (1/hh).*mapslices(sum, mat[:,(i-hh+1):i], dims=2)   # EG: 20200824
        end
    end
    return out.*temp
end
