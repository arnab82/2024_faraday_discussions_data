using QCBase
using Printf
using FermiCG
using NPZ
using InCoreIntegrals
using RDM
using JLD2


@load "data_cmf_pi_cg_4.jld2"

init_fspace = [(3,3),(3,3),(3,3),(3,3),(3,3),(3,3),(3,3)]

ecore = ints.h0

cluster_bases = FermiCG.compute_cluster_est_basis(ints, clusters,d1.a,d1.b; 
                thresh_schmidt=1e-6, thresh_orb=1e-8, thresh_ci=1e-6,
                do_embedding=true,verbose=1,init_fspace=FermiCG.FockConfig(init_fspace),delta_elec=[2,2,2,2,2,2,2],
                est_nr=1, est_max_cycles=200, est_thresh=1e-4)
clustered_ham = FermiCG.extract_ClusteredTerms(ints, clusters)
cluster_ops = FermiCG.compute_cluster_ops(cluster_bases, ints);

FermiCG.add_cmf_operators!(cluster_ops, cluster_bases, ints, d1.a, d1.b);

nroots = 1
ci_vector = FermiCG.TPSCIstate(clusters, FermiCG.FockConfig(init_fspace), R=nroots);
ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([1,1,1,1,1,1,1])] = zeros(Float64,nroots)
FermiCG.eye!(ci_vector)

eci, v = FermiCG.tps_ci_direct(ci_vector, cluster_ops, clustered_ham);

tucker = [8e-3,6e-3,4e-3,2e-3,1e-3,8e-4,6e-4,4e-4,2e-4,1e-4,8e-5]

e0b, v0b = FermiCG.tpsci_ci(v, cluster_ops, clustered_ham,
                                incremental  = true,
                                thresh_cipsi = 4e-4,
                                thresh_foi   = 1e-6,
                                thresh_asci  = -1,
				nbody        = 4,
                                max_mem_ci = 200.0);

@time e2 = FermiCG.compute_pt2_energy(v0b, cluster_ops, clustered_ham, thresh_foi=1e-8);
name = "tucker_thresh_4e4.jld2"
clustered_S2 = FermiCG.extract_S2(ci_vector.clusters)
@time s2 = FermiCG.compute_expectation_value_parallel(v0b, cluster_ops, clustered_S2)
@save "/home/arnab22/SPF-data/hexabenzocoronene/sto3g/tpsci/est_basis/tpsci_nbody4/deltaelec2/jld2/name" cluster_bases e0b v0b e2 s2

rotations = FermiCG.hosvd(v0b, cluster_ops)
for ci in clusters
    FermiCG.rotate!(cluster_ops[ci.idx], rotations[ci.idx])
    FermiCG.rotate!(cluster_bases[ci.idx], rotations[ci.idx])
    FermiCG.check_basis_orthogonality(cluster_bases[ci.idx])
end
e0b, v0b = FermiCG.tpsci_ci(v0b, cluster_ops, clustered_ham,
                                incremental  = true,
                                thresh_cipsi = 2.5e-4,
                                thresh_foi   = 1e-6,
                                thresh_asci  = -1,
				nbody        = 4,
                                max_mem_ci = 200.0);

@time e2 = FermiCG.compute_pt2_energy(v0b, cluster_ops, clustered_ham, thresh_foi=1e-8);
name = "tucker_thresh_5_2e4.jld2"
clustered_S2 = FermiCG.extract_S2(ci_vector.clusters)
@time s2 = FermiCG.compute_expectation_value_parallel(v0b, cluster_ops, clustered_S2)
@save "/home/arnab22/SPF-data/hexabenzocoronene/sto3g/tpsci/est_basis/tpsci_nbody4/deltaelec2/jld2/name" cluster_bases e0b v0b e2 s2

rotations = FermiCG.hosvd(v0b, cluster_ops)
for ci in clusters
    FermiCG.rotate!(cluster_ops[ci.idx], rotations[ci.idx])
    FermiCG.rotate!(cluster_bases[ci.idx], rotations[ci.idx])
    FermiCG.check_basis_orthogonality(cluster_bases[ci.idx])
end


e0b, v0b = FermiCG.tpsci_ci(v0b, cluster_ops, clustered_ham,
                                incremental  = true,
                                thresh_cipsi = 2e-4,
                                thresh_foi   = 1e-6,
                                thresh_asci  = -1,
				nbody        = 4,
                                max_mem_ci = 200.0);

@time e2 = FermiCG.compute_pt2_energy(v0b, cluster_ops, clustered_ham, thresh_foi=1e-8);
name = "tucker_thresh_2e4.jld2"
clustered_S2 = FermiCG.extract_S2(ci_vector.clusters)
@time s2 = FermiCG.compute_expectation_value_parallel(v0b, cluster_ops, clustered_S2)
@save "/home/arnab22/SPF-data/hexabenzocoronene/sto3g/tpsci/est_basis/tpsci_nbody4/deltaelec2/jld2/name" cluster_bases e0b v0b e2 s2


