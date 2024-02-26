using FermiCG
using JLD2
using LinearAlgebra
using Printf
using QCBase
using RDM
@load  "data_cmf_pi_cg_4.jld2"
ref_fspace = FockConfig(init_fspace)
ecore = ints.h0
@load  "tucker_thresh_5_2e4.jld2"
clustered_ham = FermiCG.extract_ClusteredTerms(ints, clusters)
cluster_ops = FermiCG.compute_cluster_ops(cluster_bases, ints);
FermiCG.add_cmf_operators!(cluster_ops, cluster_bases, ints, d1.a, d1.b);

e0, v1 = FermiCG.tps_ci_direct(v0b, cluster_ops, clustered_ham)
@time e2 = FermiCG.compute_pt2_energy(v1, cluster_ops, clustered_ham, thresh_foi=1e-8);
clustered_S2 = FermiCG.extract_S2(v1.clusters)
@time s2 = FermiCG.compute_expectation_value_parallel(v1, cluster_ops, clustered_S2)
@printf(" %5s %12.8f %12.8f %12.8f\n",1, e0[1]+ecore,e2[1]+ecore, abs(s2[1]))
@save "/home/arnab22/SPF-data/hexabenzocoronene/sto3g/tpsci/est_basis/tpsci_nbody2/deltaelec0/tpsci_tucker_nbody2_deltaelec0.jl.scr/clip/tucker_clip_0.00025.jld2" clusters cluster_bases e0 v1 e2 ecore s2

FermiCG.clip!(v0b, thresh=0.0004)
e0, v11 = FermiCG.tps_ci_direct(v0b, cluster_ops, clustered_ham)
@time e2 = FermiCG.compute_pt2_energy(v11, cluster_ops, clustered_ham, thresh_foi=1e-8);
clustered_S2 = FermiCG.extract_S2(v11.clusters)
@time s2 = FermiCG.compute_expectation_value_parallel(v11, cluster_ops, clustered_S2)
@printf(" %5s %12.8f %12.8f %12.8f\n",1, e0[1]+ecore,e2[1]+ecore, abs(s2[1]))
@save "/home/arnab22/SPF-data/hexabenzocoronene/sto3g/tpsci/est_basis/tpsci_nbody2/deltaelec0/tpsci_tucker_nbody2_deltaelec0.jl.scr/clip/tucker_clip_0.0004.jld2" clusters cluster_bases e0 v11 e2 ecore s2


FermiCG.clip!(v11, thresh=0.0006)
e0, v2 = FermiCG.tps_ci_direct(v11, cluster_ops, clustered_ham)
@time e2 = FermiCG.compute_pt2_energy(v2, cluster_ops, clustered_ham, thresh_foi=1e-8);
clustered_S2 = FermiCG.extract_S2(v2.clusters)
@time s2 = FermiCG.compute_expectation_value_parallel(v2, cluster_ops, clustered_S2)
@printf(" %5s %12.8f %12.8f %12.8f\n",1, e0[1]+ecore,e2[1]+ecore, abs(s2[1]))
@save "/home/arnab22/SPF-data/hexabenzocoronene/sto3g/tpsci/est_basis/tpsci_nbody2/deltaelec0/tpsci_tucker_nbody2_deltaelec0.jl.scr/clip/tucker_clip_0.0006.jld2" clusters cluster_bases e0 v2 e2 ecore s2


FermiCG.clip!(v2, thresh=0.0008)
e0, v3 = FermiCG.tps_ci_direct(v2, cluster_ops, clustered_ham)
@time e2 = FermiCG.compute_pt2_energy(v3, cluster_ops, clustered_ham, thresh_foi=1e-8);
clustered_S2 = FermiCG.extract_S2(v3.clusters)
@time s2 = FermiCG.compute_expectation_value_parallel(v3, cluster_ops, clustered_S2)
@printf(" %5s %12.8f %12.8f %12.8f\n",1, e0[1]+ecore,e2[1]+ecore, abs(s2[1]))
@save "/home/arnab22/SPF-data/hexabenzocoronene/sto3g/tpsci/est_basis/tpsci_nbody2/deltaelec0/tpsci_tucker_nbody2_deltaelec0.jl.scr/clip/tucker_clip_0.0008.jld2" clusters cluster_bases e0 v3 e2 ecore s2


FermiCG.clip!(v2, thresh=0.001)
e0, v4 = FermiCG.tps_ci_direct(v3, cluster_ops, clustered_ham)
@time e2 = FermiCG.compute_pt2_energy(v4, cluster_ops, clustered_ham, thresh_foi=1e-8);
clustered_S2 = FermiCG.extract_S2(v4.clusters)
@time s2 = FermiCG.compute_expectation_value_parallel(v4, cluster_ops, clustered_S2)
@printf(" %5s %12.8f %12.8f %12.8f\n",1, e0[1]+ecore,e2[1]+ecore, abs(s2[1]))
@save "/home/arnab22/SPF-data/hexabenzocoronene/sto3g/tpsci/est_basis/tpsci_nbody2/deltaelec0/tpsci_tucker_nbody2_deltaelec0.jl.scr/clip/tucker_clip_0.001.jld2" clusters cluster_bases e0 v4 e2 ecore s2

