function DINEOF_pmQC(X,XA,errmap,musquare)

# poor mans quality check. Here XA is the EOF filtered field. Normally for consistency the OI-EOF analysed field should be used 
# error maps is variance of error
# TODO: add cloud and value indicators as in Aida's paper. Exploit multidimensional approach of bubble generation to make the ND generalisation
#
validpoints=.!isnan.(X)

    @show size(X[validpoints])

    o=abs.(X[validpoints].-XA[validpoints])./sqrt.(max.(musquare.-errmap[validpoints],1E-10*musquare))

    om=median(o)

    del=1.4826*median(abs.(o.-om))
    
    OO=zeros(Float64,size(X))
    OO[validpoints]=abs.(o.-om)./del
    
    return OO

end


