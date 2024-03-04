using JLD2
using FermiCG
using Plots
using LinearAlgebra
using Printf


function run()
    @load("/Users/arnab/arnab/workspace/2024_faraday_discussions_data/hexabenzocoronene/correlations_ccpvdz/fermicluster_integrals/correlation/est/tucker_thresh_5_1e4.jld2")
    @load("/Users/arnab/arnab/workspace/2024_faraday_discussions_data/hexabenzocoronene/correlations_ccpvdz/fermicluster_integrals/hexabenzocoronene_old.jld2")
    init_fspace = [(3,3),(3,3),(3,3),(3,3),(3,3),(3,3),(3,3)]
    display(v0b)
    ecore = ints.h0

    cluster_bases = FermiCG.compute_cluster_est_basis(ints, clusters,d1.a,d1.b; 
                thresh_schmidt=1e-6, thresh_orb=1e-8, thresh_ci=1e-6,
                do_embedding=true,verbose=1,init_fspace=FermiCG.FockConfig(init_fspace),
                est_nr=1, est_max_cycles=200, est_thresh=1e-4)
    clustered_ham = FermiCG.extract_ClusteredTerms(ints, clusters)
    cluster_ops = FermiCG.compute_cluster_ops(cluster_bases, ints);

    FermiCG.add_cmf_operators!(cluster_ops, cluster_bases, ints, d1.a, d1.b);

    nroots = 1
    ci_vector = FermiCG.TPSCIstate(clusters, FermiCG.FockConfig(init_fspace), R=nroots);
    ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([1,1,1,1,1,1,1])] = zeros(Float64,nroots)
    FermiCG.eye!(ci_vector)
    # @load "cmfstate.jld2"
    eci, v = FermiCG.tps_ci_direct(ci_vector, cluster_ops, clustered_ham);
    @save "/Users/arnab/arnab/workspace/2024_faraday_discussions_data/hexabenzocoronene/correlations_ccpvdz/fermicluster_integrals/correlation/est/tpsstate_initial.jld2" v
    # @load "tpsstate_initial.jld2" 
    
    # eci, v0b_ = FermiCG.tps_ci_direct(v0b, cluster_ops, clustered_ham);
    display(v)
    C_f=correlation_functions(v0b, v)
    # C_f_=correlation_functions(v0b, ci_vector)
    @save "/Users/arnab/arnab/workspace/2024_faraday_discussions_data/hexabenzocoronene/correlations_ccpvdz/fermicluster_integrals/correlation/est/correlation_QO.jld2" C_f
end

run()