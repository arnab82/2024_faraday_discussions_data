using FermiCG, NPZ, JLD2
using LinearAlgebra
using Printf
using QCBase
using RDM
using ClusterMeanField

@load  "/home/nbraunsc/FermiCG-data/faraday/tetracene_dimer_1_2/cmf_diis_herr.jld2"

ints = deepcopy(ints_cmf)
C = deepcopy(C_cmf);
ecore = ints.h0

M = 150

ref_fock = FockConfig(init_fspace)
#
# Build Cluster basis
cluster_bases = FermiCG.compute_cluster_eigenbasis_spin(ints, clusters, d1, [5,5], ref_fock, max_roots=M, verbose=1);
#
# Build ClusteredOperator
clustered_ham = FermiCG.extract_ClusteredTerms(ints, clusters);

#
# Build Cluster Operators
cluster_ops = FermiCG.compute_cluster_ops(cluster_bases, ints);

#
# Add cmf hamiltonians for doing MP-style PT2 
FermiCG.add_cmf_operators!(cluster_ops, cluster_bases, ints, d1.a, d1.b, verbose=0);

nroots = 10
ci_vector = FermiCG.TPSCIstate(clusters, FermiCG.FockConfig(init_fspace), R=nroots);
#ci_vector = FermiCG.add_spin_focksectors(ci_vector)

# Add the lowest energy single exciton to basis
ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([1,1])] = zeros(Float64,nroots)
ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([2,1])] = zeros(Float64,nroots)
ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([1,2])] = zeros(Float64,nroots)
ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([3,1])] = zeros(Float64,nroots)
ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([1,3])] = zeros(Float64,nroots)
ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([4,1])] = zeros(Float64,nroots)
ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([1,4])] = zeros(Float64,nroots)

# TT states ms=0
ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig([2,2])] = zeros(Float64,nroots)

# Spin-flip states
fspace_0 = FermiCG.FockConfig(init_fspace)

## ba
tmp_fspace = FermiCG.replace(fspace_0, (1,2), ([6,4],[4,6]))
FermiCG.add_fockconfig!(ci_vector, tmp_fspace)
ci_vector[tmp_fspace][FermiCG.ClusterConfig([1,1])] = zeros(Float64,nroots)

## ab
tmp_fspace = FermiCG.replace(fspace_0, (1,2), ([4,6],[6,4]))
FermiCG.add_fockconfig!(ci_vector, tmp_fspace)
ci_vector[tmp_fspace][FermiCG.ClusterConfig([1,1])] = zeros(Float64,nroots)

FermiCG.eye!(ci_vector)

e_guess, v_guess = FermiCG.tps_ci_direct(ci_vector, cluster_ops, clustered_ham);
@time e2_guess = FermiCG.compute_pt2_energy(v_guess, cluster_ops, clustered_ham, thresh_foi=1e-8);
@save "/home/nbraunsc/FermiCG-data/faraday/tetracene_dimer_1_2/H_guess.jld2" e_guess v_guess e2_guess ecore

e0, v0 = FermiCG.tpsci_ci(v_guess, cluster_ops, clustered_ham,
                          thresh_asci =-1,     # Threshold of P-space configs to search from
                          thresh_foi  =1e-5,    # Threshold for keeping terms when defining FOIS
                          thresh_cipsi=0.0004, # Threshold for adding to P-space
                          max_iter=10,
                          max_mem_ci=260.0);

@time e2 = FermiCG.compute_pt2_energy(v0, cluster_ops, clustered_ham, thresh_foi=1e-8);
clustered_S2 = FermiCG.extract_S2(ci_vector.clusters)
@time s2 = FermiCG.compute_expectation_value_parallel(v0, cluster_ops, clustered_S2)
@save "/home/nbraunsc/FermiCG-data/faraday/tetracene_dimer_1_2/thresh_0004.jld2" clusters cluster_bases e0 v0 e2 s2 ecore

