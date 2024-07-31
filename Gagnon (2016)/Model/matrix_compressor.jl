# Function to prepare input data matrices for backwards looking expectations
function matrix_compressor(i::Int64, hh::Int64, mat::data_t)
    """
       i is the iteration of the for loop this function will be called in,
       it sets which column is the last one selected of the input data

       hh is the number of columns backwards which we want to select

       mat is the input matrix
    """
    out = ones(size(mat)) # always return a matrix of the same size as input
    if (i < hh)
        # We have to duplicate the first column the requisite number of times
        # for example, if i = 2 but hh = 4 we will need to get the average of
        # 3 column 1s and 1 column 2
        temp = ((hh - i)/hh).*mat[:, 1] + 1/hh.*mapslices(sum, mat[:, 2:i], 2)
    else
        # In the case where i >= hh we just take a straight average
        # over the columns (i-hh+1:i)
        temp = (1/hh).*mapslices(sum, mat[:,(i-hh+1):i], 2)
    end
    return out.*temp
end
