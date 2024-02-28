using FermiCG, JLD2

@load "out.jld2"



function run(v::TPSCIstate{T,N,R}, ints, cluster_basis, opstring::String) where {T,N,R}

    cluster_ops = FermiCG.compute_cluster_ops(cluster_bases, ints);

    cf1 = zeros(T,N)
    cf2 = zeros(T,N,N)

    for (fock,configs) in v.data
        prob = 0
        fock_trans = fock - fock
        for (config, coeff) in configs
            prob += coeff[root]*coeff[root]

        end

        for ci in v.clusters
            n1[root][ci.idx] += prob * (fock[ci.idx][1] + fock[ci.idx][2])
            sz1[root][ci.idx] += prob * (fock[ci.idx][1] - fock[ci.idx][2]) / 2
            for cj in v.clusters
                ci.idx <= cj.idx || continue
                n2[root][ci.idx, cj.idx] += prob * (fock[ci.idx][1] + fock[ci.idx][2]) * (fock[cj.idx][1] + fock[cj.idx][2])
                sz2[root][ci.idx, cj.idx] += prob * (fock[ci.idx][1] - fock[ci.idx][2]) * (fock[cj.idx][1] - fock[cj.idx][2]) / 4
                n2[root][cj.idx, ci.idx] = n2[root][ci.idx, cj.idx]
                sz2[root][cj.idx, ci.idx] = sz2[root][ci.idx, cj.idx]
            end
        end
    end
end

run(v0a, 1)
