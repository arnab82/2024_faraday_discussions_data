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

@load "thresh_0004_no_clusterbases.jld2"
@load "H_guess_0004.jld2"
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

Heff_bare_tt = abs.(Heff_bare_tt - diagm(diag(Heff_bare_tt)))
Heff_exact_tt = abs.(Heff_exact_tt - diagm(diag(Heff_exact_tt)))
max_val = max(maximum(abs.(Heff_bare_tt)), maximum(abs.(Heff_exact_tt)))

#plotd = heatmap(Heff*conversion; color=palette([:teal, :white, :orange], 10), aspect_ratio=1, dpi=300, size=(300,300), right_margin = 10Plots.mm)
plotd = heatmap(Heff_bare_tt; color=palette(:bone_1, rev=true), aspect_ratio=1, dpi=300, size=(300,300), right_margin = 10Plots.mm,  
                clims=(0, max_val), ticks = false,xaxis=false,yaxis=false, 
                xlims = (0.5,31.5), ylims = (0.5,31.5), title="Heff_bare")
m, n = size(Heff_bare_tt)
vline!(0.5:(n+0.5), c=:white, label=false)
hline!(0.5:(m+0.5), c=:white, label=false)


savefig(plotd,"Heff_bare.png")



plotd = heatmap(Heff_exact_tt; color=palette(:bone_1, rev=true), aspect_ratio=1, dpi=300, size=(300,300), right_margin = 10Plots.mm,  
                clims=(0, max_val), ticks = false, xaxis=false,yaxis=false,
                xlims = (0.5,31.5), ylims = (0.5,31.5), title="Heff_exact")
#plotd = heatmap(Heff*27.2114079527*1000; color=palette(:Blues, 10), aspect_ratio=1, dpi=300, size=(300,300), right_margin = 10Plots.mm, xlims = (0,5), ylims = (0,5))
m, n = size(Heff_bare_tt)
vline!(0.5:(n+0.5), c=:white, label=false)
hline!(0.5:(m+0.5), c=:white, label=false)
savefig(plotd,"Heff_exact.png")



