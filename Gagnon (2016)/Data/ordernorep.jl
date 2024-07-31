# This function orders the data couples {x[i],y[i]}, i=1:n, such that x[i]<=x[i+1]. If there are
# repeated values of x[i], then the function keeps on one couple, with the value of y[i] set 
# to the average y[j] across all values of y[j] such that x[j]==x[i].
# 
# Note: If the average of repeated numbers is in a different type than the y inputed, then 
# one may get a "Load error inexct... while loading".  The solution is to ensure that y is 
# in a floating format (say by replacing "1" with "1.0" in the elements of y).


function ordernorep(x,y)

    # Ordering the data and dealing with repetitions
#    println("Ordering the data and dealing with repetitions")
    xorder = sortperm(x)
    x = x[xorder]
    y = y[xorder]
    
    # Checking that there are no repeated values of x
#    println("Checking that there are no repeated values of x")
    chk0=sum(x[1:(end-1)] .== x[2:end])
#    @printf "chk0= %f\n" chk0
    while (sum(x[1:(end-1)] .== x[2:end])>0)
#        println("Here 1")
        tmpn=1:(length(x)-1)
        tmpndiff=tmpn[x[1:(end-1)].==x[2:end]]
        chk1=length(tmpndiff)
#        @printf "chk1 = %f\n" chk1
#        @printf "first rep: %f\n" tmpndiff[1]
#        println("Here 2")
#        @printf "lx = %f; ly = %f\n" length(x) length(y) 
#        for aa=1:length(x)
#            @printf "x[%d]=%f y[%d]=%f\n" aa x[aa] aa y[aa]
#        end
        #tmpy=mean(y)
        #tmpy=(x[tmpndiff[2]])
        tmpy=mean(y[x.==x[tmpndiff[1]]])
 #       @printf "tmpy=%f    tmpndiff=%f  xx=%f \n" tmpy tmpndiff[1] x[tmpndiff[1]]
        y[x.==x[tmpndiff[1]]] .= tmpy
 #       println("Here 3")
 #       for aa=1:length(x)
 #           @printf "x[%d]=%f y[%d]=%f\n" aa x[aa] aa y[aa]
 #       end
        # Trimming x and y

        tmpn2=1:length(x)
#        chk2=(x.!=x[tmpndiff[1]])
#        chk2=(tmpn2.==tmpndiff[1])
#        chk2=(x.!=x[tmpndiff[1]]) | (tmpn2.==tmpndiff[1])
#        @printf "chk2 = %f\n" length(chk2)
#        for aa=1:length(x)
#            @printf "x[%d]=%f chk2[%d]=%f\n" aa x[aa] aa chk2[aa]
#        end
        
        y=y[(x.!=x[tmpndiff[1]]) .| (tmpn2.==tmpndiff[1])]
        x=x[(x.!=x[tmpndiff[1]]) .| (tmpn2.==tmpndiff[1])]
 #       @printf "lx = %f; ly = %f\n" length(x) length(y) 
 #       for aa=1:length(x)
 #           @printf "x[%d]=%f y[%d]=%f\n" aa x[aa] aa y[aa]
 #       end
    end

    return x , y

end
