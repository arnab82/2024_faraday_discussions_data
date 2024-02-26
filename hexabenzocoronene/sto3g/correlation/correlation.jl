using JLD2
using FermiCG
using Plots
using LinearAlgebra
using Printf


function run()
    @load("/Users/arnab/arnab/workspace/2024_faraday_discussions_data/hexabenzocoronene/sto3g/correlation/tucker_thresh_5_2e4.jld2")
    @load("/Users/arnab/arnab/workspace/2024_faraday_discussions_data/hexabenzocoronene/sto3g/correlation/hexabenzocoronene_sto3g.jld2")
    cluster_ops = FermiCG.compute_cluster_ops(cluster_bases, ints);

    FermiCG.add_cmf_operators!(cluster_ops, cluster_bases, ints, d1.a, d1.b);

    display(v0b)
    C_f=correlation_functions(v0b, cluster_ops)
    @save "/Users/arnab/arnab/workspace/2024_faraday_discussions_data/hexabenzocoronene/sto3g/correlation/correlation.jld2" C_f
end

run()



