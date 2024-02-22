using FermiCG, NPZ, JLD2
using LinearAlgebra
using Printf
using QCBase
using RDM

@load  "cmf_diis_kekule.jld2"

ints = deepcopy(ints_cmf)
C = deepcopy(C_cmf);
ecore = ints.h0

M = 150

ref_fock = FockConfig(init_fspace)
#
# Build Cluster basis
tmp = zeros(Int, 21)
fill!(tmp, 2)

cluster_bases = FermiCG.compute_cluster_eigenbasis_spin(ints, clusters, d1, tmp, ref_fock, max_roots=M, verbose=1);
#
# Build ClusteredOperator
clustered_ham = FermiCG.extract_ClusteredTerms(ints, clusters);

#
# Build Cluster Operators
cluster_ops = FermiCG.compute_cluster_ops(cluster_bases, ints);

#
# Add cmf hamiltonians for doing MP-style PT2 
FermiCG.add_cmf_operators!(cluster_ops, cluster_bases, ints, d1.a, d1.b, verbose=0);

nroots = 1
ci_vector = FermiCG.TPSCIstate(clusters, FermiCG.FockConfig(init_fspace), R=nroots);
#ci_vector = FermiCG.add_spin_focksectors(ci_vector)

# Add the lowest energy single exciton to basis
ci_vector[FermiCG.FockConfig(init_fspace)][FermiCG.ClusterConfig(fill!(tmp,1))] = [1.0]

#FermiCG.eye!(ci_vector)

e_guess, v_guess = FermiCG.tps_ci_direct(ci_vector, cluster_ops, clustered_ham);
@time e2_guess = FermiCG.compute_pt2_energy(v_guess, cluster_ops, clustered_ham, thresh_foi=1e-8);
@save "H_guess.jld2" e_guess v_guess e2_guess ecore

e0a, v0a = FermiCG.tpsci_ci(v_guess, cluster_ops, clustered_ham,
                            incremental  = true,
                            thresh_cipsi = 1e-3,
                            thresh_foi   = 1e-5,
                            thresh_asci  = -1, 
                            max_mem_ci = 100.0);


rotations = FermiCG.hosvd(v0a, cluster_ops)
for ci in clusters
    FermiCG.rotate!(cluster_ops[ci.idx], rotations[ci.idx])
    FermiCG.rotate!(cluster_bases[ci.idx], rotations[ci.idx])
    FermiCG.check_basis_orthogonality(cluster_bases[ci.idx])
end

tucker = [.8e-3, .6e-3, .4e-3, .2e-3, .1e-3]

for i in tucker
    e0b, v0b = FermiCG.tpsci_ci(ci_vector, cluster_ops, clustered_ham,
                                incremental  = true,
                                thresh_cipsi = i,
                                thresh_foi   = 1e-5,
                                thresh_asci  = -1,
                                max_mem_ci = 100.0);

    @time e2 = FermiCG.compute_pt2_energy(v0b, cluster_ops, clustered_ham, thresh_foi=1e-8);
    name = "tucker_thresh_"*string(i)*".jld2"

    println()
    println("	*======TPSCI results======*")
    @printf("TCI Thresh: %8.6f  Dim:%8d\n",i,size(v0b)[1])
    println()
    @printf("TCI %5s %12s %12s\n", "Root", "E(0)", "E(2)") 
    for r in 1:nroots
        @printf("TCI %5s %12.8f %12.8f\n",r, e0b[r] + ecore, e0b[r] + e2[r] + ecore)
    end

    clustered_S2 = FermiCG.extract_S2(ci_vector.clusters)

    println()
    println("	*======TPSCI S2 results======*")
    @printf(" %-50s", "Compute FINAL S2 expectation values: ")
    @time s2 = FermiCG.compute_expectation_value_parallel(v0b, cluster_ops, clustered_S2)

    @printf(" %5s %12s %12s\n", "Root", "Energy", "S2") 
    for r in 1:nroots
        @printf(" %5s %12.8f %12.8f\n",r, e0b[r]+ecore, abs(s2[r]))
    end
    @save ""*string(name) cluster_bases e0b v0b e2 s2 ecore

    rotations = FermiCG.hosvd(v0b, cluster_ops)
    for ci in clusters
        FermiCG.rotate!(cluster_ops[ci.idx], rotations[ci.idx])
        FermiCG.rotate!(cluster_bases[ci.idx], rotations[ci.idx])
        FermiCG.check_basis_orthogonality(cluster_bases[ci.idx])
    end
end


