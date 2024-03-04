using JLD2
using FermiCG
using Plots
using LinearAlgebra
using Printf


#@load("data_tpsci_013.jld2")
#alist = [i for i in 9:15]
#dlist  = [3,5,7,9,11,13,15]


function run()
    for geom in 5:30
        file1 = @sprintf("data_tpsci_%03i.jld2", geom)
        @eval @load $file1
        vproj = FermiCG.overlap(vguess, v0);
        weights = diag(vproj'*vproj)

        @printf("geom_%03i: ",geom)
        [@printf("%12.8f ", i^.5) for i in weights]
        println()
    end
end

run()
