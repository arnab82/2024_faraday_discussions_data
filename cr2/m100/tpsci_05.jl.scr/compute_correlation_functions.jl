using FermiCG, JLD2

@load "out.jld2"


cluster_ops = FermiCG.compute_cluster_ops(cluster_bases, ints);

cf = correlation_functions(v0a, cluster_ops);
@save "cf.jld2"

#function run(v::TPSCIstate{T,N,R}, ints, cluster_bases) where {T,N,R}
#
#    cluster_ops = FermiCG.compute_cluster_ops(cluster_bases, ints);
#
#    cf = correlation_functions(v, cluster_ops);
#    @save "cf.jld2"
#
#end
#
#run(v0a, ints, cluster_bases)
