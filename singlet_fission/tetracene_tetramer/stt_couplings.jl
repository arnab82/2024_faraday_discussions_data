using JLD2
using FermiCG
using Plots
using LinearAlgebra
using Printf

function make_H(diabats::TPSCIstate, adiabats::TPSCIstate, ein)
    vproj = FermiCG.overlap(diabats, adiabats)
    println(" Vproj:")
    display(vproj)

    e = ein

    #e = e .- (sum(e)/length(e))

    S = vproj'*vproj
    println(" S:")
    display(S)
    vproj_orth = vproj * inv(sqrt(S))

    Heff = vproj_orth * diagm(e) * vproj_orth'

    println(" Heff")
    display(Heff)

    # shift diagonal to center at zero
    return Heff
end

function run()
    #hartree to meV
    conversion = 27.2114079527*1000

    @load "thresh_spin_0.0004.jld2"
    @load "H_guess_spin.jld2"
    @load "v_model_tetramer.jld2"

    display(FermiCG.get_vector(v_guess))
    v_barett = v_guess
    e_barett = e_guess
    e_exacttt = (e0)
    #e_exact = (e2 + e0)
    v_exacttt = deepcopy(v0)
    v_modeltt = deepcopy(v_guess)
    v_diabattt = deepcopy(ci_vector)
    FermiCG.eye!(v_diabattt)

    #plot tetramer data now in the same formatting
    # First make bare 
    println()
    println(" Make Model Space Bare Hamiltonian")
    Heff_bare_tt = make_H(v_diabattt, v_barett, e_barett)
    Heff_bare_tt .*= conversion

    # Now exact
    println()
    println(" Make Model Space Effective Hamiltonian")
    Heff_exact_tt = make_H(v_diabattt, v_exacttt, e_exacttt)
    Heff_exact_tt .*= conversion

    #make S-TT Heff
    #(1,2), (1,4)(2,3), (1,3)(2,4), (3,4)
    stt = [1,10,11,12,13,14,20,26,16,22,28,17,23,29,15,21,27,18,24,30,19,25,31]
    Heff_exact_tt = Heff_exact_tt[stt, stt]
    Heff_bare_tt = Heff_bare_tt[stt, stt]

    #Block Diagonalize S1S2 and TT
    U = zeros(length(stt), length(stt))
    U[1,1]=1.0
    U[2:5,2:5].=Matrix(1.0I, 4, 4)
    e, v = eigen(Heff_bare_tt[6:8, 6:8])
    U[6:8,6:8] .= v
    e, v = eigen(Heff_bare_tt[9:11, 9:11])
    U[9:11,9:11] .= v
    e, v = eigen(Heff_bare_tt[12:14, 12:14])
    U[12:14,12:14] .= v
    e, v = eigen(Heff_bare_tt[15:17, 15:17])
    U[15:17,15:17] .= v
    e, v = eigen(Heff_bare_tt[18:21, 18:21])
    U[18:21,18:21] .= v
    e, v = eigen(Heff_bare_tt[21:23,21:23])
    U[21:23, 21:23] .= v

    Heff_bare_tt = U'*Heff_bare_tt*U
    display(Heff_bare_tt)

    U = zeros(length(stt), length(stt))
    U[1,1]=1.0
    U[2:5,2:5].=Matrix(1.0I, 4, 4)
    e, v = eigen(Heff_exact_tt[6:8, 6:8])
    U[6:8,6:8] .= v
    e, v = eigen(Heff_exact_tt[9:11, 9:11])
    U[9:11,9:11] .= v
    e, v = eigen(Heff_exact_tt[12:14, 12:14])
    U[12:14,12:14] .= v
    e, v = eigen(Heff_exact_tt[15:17, 15:17])
    U[15:17,15:17] .= v
    e, v = eigen(Heff_exact_tt[18:21, 18:21])
    U[18:21,18:21] .= v
    e, v = eigen(Heff_exact_tt[21:23,21:23])
    U[21:23, 21:23] .= v

    Heff_exact_tt = U'*Heff_exact_tt*U
    display(Heff_exact_tt)

    Heff_exact_tt = Heff_exact_tt[[1,2,3,4,5,6,9,12,15,18,21], [1,2,3,4,5,6,9,12,15,18,21]]
    Heff_bare_tt = Heff_bare_tt[[1,2,3,4,5,6,9,12,15,18,21], [1,2,3,4,5,6,9,12,15,18,21]]
    println("Heff bare")
    display(Heff_bare_tt)
    println("Heff exact")
    display(Heff_exact_tt)

    #Heff_bare_tt = abs.(Heff_bare_tt - diagm(diag(Heff_bare_tt)))
    #Heff_exact_tt = abs.(Heff_exact_tt - diagm(diag(Heff_exact_tt)))

    Heff_bare_tt = Heff_bare_tt - diagm(diag(Heff_bare_tt))
    Heff_exact_tt = Heff_exact_tt - diagm(diag(Heff_exact_tt))

    max_val = max(maximum(abs.(Heff_bare_tt)), maximum(abs.(Heff_exact_tt)))

    plotd = heatmap(Heff_bare_tt; color=palette(:RdGy_9, 100), aspect_ratio=1, dpi=300, size=(300,300), right_margin = 10Plots.mm,  
                    clims=(-max_val, max_val), ticks = false,xaxis=false,yaxis=false, 
                    xlims = (0.5,11.5), ylims = (0.5,11.5))
    m, n = size(Heff_bare_tt)
    vline!(0.5:(n+0.5), c=:grey, label=false)
    hline!(0.5:(m+0.5), c=:grey, label=false)


    savefig(plotd,"Heff_bare_stt.png")



    plotd = heatmap(Heff_exact_tt; color=palette(:RdGy_9, 100), aspect_ratio=1, dpi=300, size=(300,300), right_margin = 10Plots.mm,  
                    clims=(-max_val, max_val), ticks = false, xaxis=false,yaxis=false,
                    xlims = (0.5,11.5), ylims = (0.5,11.5))
    m, n = size(Heff_exact_tt)
    vline!(0.5:(n+0.5), c=:grey, label=false)
    hline!(0.5:(m+0.5), c=:grey, label=false)
    savefig(plotd,"Heff_exact_stt.png")

    tmpbare = deepcopy(Heff_bare_tt)
    tmpeff = deepcopy(Heff_exact_tt)

    #(1,2), (1,4)(2,3), (1,3)(2,4), (3,4)
    dimers = [[1,2,3,6],[1,2,5,7], [1,3,4,8], [1,2,4,9], [1,3,5,10], [1,4,5,11]] 
    labels = ("12", "14", "23", "13", "24", "34")

    max_val = 0
    for i in 1:6
        tte = tmpeff[dimers[i],dimers[i]]
        max_val = max(max_val, maximum(abs.(tte)))
    end

    for i in 1:6
        ttb = tmpbare[dimers[i], dimers[i]]
        tte = tmpeff[dimers[i],dimers[i]]

        plotd = heatmap(ttb; color=palette(:RdGy_9, 100), aspect_ratio=1, dpi=300, size=(300,300), right_margin = 10Plots.mm,  
                        clims=(-max_val, max_val), ticks = false,xaxis=false,yaxis=false, 
                        xlims = (0.5,4.5), ylims = (0.5,4.5))
        m, n = size(ttb)
        vline!(0.5:(n+0.5), c=:grey, label=false)
        hline!(0.5:(m+0.5), c=:grey, label=false)
        savefig(plotd,"plots/y_Heff_bare_"*labels[i]*".png")

        plotd = heatmap(tte; color=palette(:RdGy_9, 100), aspect_ratio=1, dpi=300, size=(300,300), right_margin = 10Plots.mm,  
                        clims=(-max_val, max_val), ticks = false,xaxis=false,yaxis=false, 
                        xlims = (0.5,4.5), ylims = (0.5,4.5))
        m, n = size(tte)
        vline!(0.5:(n+0.5), c=:grey, label=false)
        hline!(0.5:(m+0.5), c=:grey, label=false)
        savefig(plotd,"plots/y_Heff_exact_"*labels[i]*".png")
    end
end
run()

