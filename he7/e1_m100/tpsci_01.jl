using QCBase
using Printf
using FermiCG
using NPZ
using InCoreIntegrals
using RDM
using JLD2


function run_tpsci_scan()
    for geom in reverse(5:30)
        f_out = @sprintf("%03i.out",geom)
        f_err = @sprintf("%03i.err",geom)

        redirect_stdio(stdout=f_out, stderr=f_err) do

            file1 = @sprintf("../data_cmf_%03i.jld2", geom)
            println(" CMF File: ", file1)
            @eval @load $file1
            println(" Enuc: ", ints.h0)

            M = 100 


            ref_fspace = FockConfig(init_fspace)
            ecore = ints.h0

            cluster_bases = FermiCG.compute_cluster_eigenbasis_spin(ints, clusters, d1, [1,1,1,1,1,1,1], ref_fspace, max_roots=M, verbose=1);

            clustered_ham = FermiCG.extract_ClusteredTerms(ints, clusters)
            cluster_ops = FermiCG.compute_cluster_ops(cluster_bases, ints);

            FermiCG.add_cmf_operators!(cluster_ops, cluster_bases, ints, d1.a, d1.b);

            nroots = 15 

            ci_vector = FermiCG.single_excitonic_basis(clusters, ref_fspace, R=nroots, Nk=3)

            eguess, vguess, H = FermiCG.tps_ci_direct(ci_vector, cluster_ops, clustered_ham);
            @time e2guess = FermiCG.compute_pt2_energy(vguess, cluster_ops, clustered_ham, thresh_foi=1e-8);

            e0, v0 = FermiCG.tpsci_ci(vguess, cluster_ops, clustered_ham,
                                        incremental  = true,
                                        thresh_cipsi = 1e-3,
                                        thresh_foi   = 1e-5,
                                        thresh_spin  = 1e-5, 
                                        max_mem_ci = 50.0);
            
            e0, v0, H = FermiCG.tps_ci_direct(v0, cluster_ops, clustered_ham);


            @time e2 = FermiCG.compute_pt2_energy(v0, cluster_ops, clustered_ham, thresh_foi=1e-8);
    
            e0 .+= ints.h0
            eguess .+= ints.h0
            H .+= Matrix(I*ints.h0, size(H)) 

            print("*EGuess: ")
            [@printf("%12.8f ", i) for i in eguess]
            println()

            print("*E2Guess: ")
            [@printf("%12.8f ", i) for i in eguess.+e2guess]
            println()

            print("*E0: ")
            [@printf("%12.8f ", i) for i in e0]
            println()

            print("*E2: ")
            [@printf("%12.8f ", i) for i in e2.+e0]
            println()

            # Save H for plotting
            # 
            @save @sprintf("data_tpsci_%03i.jld2",geom) vguess v0  eguess e0 e2 e2guess H
        end
    end
end

run_tpsci_scan()

