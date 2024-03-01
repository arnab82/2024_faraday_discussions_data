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


#hartree to meV
conversion = 27.2114079527*1000

@load "thresh_spin_0.0004.jld2"
@load "H_guess_spin.jld2"
@load "v_model_tetramer.jld2"

display(FermiCG.get_vector(v_guess))
v_barett = v_guess
e_barett = e_guess
e_exacttt = (e0)
#e_exacttt = (e2 + e0)
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

#show only singlets and biexitons
stt = [1,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31]
Heff_bare_tt = Heff_bare_tt - diagm(diag(Heff_bare_tt))
Heff_exact_tt = Heff_exact_tt - diagm(diag(Heff_exact_tt))

Heff_exact_tt = Heff_exact_tt[stt, stt]
Heff_bare_tt = Heff_bare_tt[stt, stt]
max_val = max(maximum(abs.(Heff_bare_tt)), maximum(abs.(Heff_exact_tt)))

#plotd = heatmap(Heff_bare_tt; color=palette(:bone_1, rev=true), aspect_ratio=1, dpi=300, size=(300,300), right_margin = 10Plots.mm,  
plotd = heatmap(Heff_bare_tt; color=palette(:RdGy_9, 100), aspect_ratio=1, dpi=300, size=(300,300), right_margin = 10Plots.mm,  
                clims=(-max_val, max_val), ticks = false,xaxis=false,yaxis=false, 
                xlims = (0.5,23.5), ylims = (0.5,23.5), yflip=true)
m, n = size(Heff_bare_tt)
vline!(0.5:(n+0.5), c=:white, label=false)
hline!(0.5:(m+0.5), c=:white, label=false)
tmp = [0.5, 1.5, 5.5, 23.5]
vline!(tmp, c=:grey, label=false)
hline!(tmp, c=:grey, label=false)


savefig(plotd,"plots/Heff_bare_spin.png")

#plotd = heatmap(Heff_exact_tt; color=palette(:bone_1, rev=true), aspect_ratio=1, dpi=300, size=(300,300), right_margin = 10Plots.mm,  
plotd = heatmap(Heff_exact_tt; color=palette(:RdGy_9, 100), aspect_ratio=1, dpi=300, size=(300,300), right_margin = 10Plots.mm,  
                clims=(-max_val, max_val), ticks = false, xaxis=false,yaxis=false,
                xlims = (0.5,23.5), ylims = (0.5,23.5), yflip=true)
m, n = size(Heff_bare_tt)
vline!(0.5:(n+0.5), c=:white, label=false)
hline!(0.5:(m+0.5), c=:white, label=false)

tmp = [0.5, 1.5, 5.5, 23.5]
vline!(tmp, c=:grey, label=false)
hline!(tmp, c=:grey, label=false)
savefig(plotd,"plots/Heff_exact_spin.png")



