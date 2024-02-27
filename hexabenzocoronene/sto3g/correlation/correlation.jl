using JLD2
using FermiCG
using Plots
using LinearAlgebra
using Printf


function run()
    @load "tucker_thresh_5_2e4.jld2"
    @load "hexabenzocoronene_sto3g.jld2"
    #cluster_ops = FermiCG.compute_cluster_ops(cluster_bases, ints);

    #FermiCG.add_cmf_operators!(cluster_ops, cluster_bases, ints, d1.a, d1.b);

    cmf_state = TPSCIstate(clusters, FockConfig(init_fspace))

    display(v0b)
    C_f = correlation_functions(v0b, cmf_state)
    @save "correlation.jld2" C_f
end

run()



