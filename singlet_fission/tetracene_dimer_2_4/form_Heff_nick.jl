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

@load "thresh_0004_nocluster_bases.jld2"
@load "../tetracene_dimer_1_4/v_model_dimer.jld2"
@load "H_guess.jld2"

display(FermiCG.get_vector(v_guess))
v_bare = v_guess
e_bare = e_guess
e_exact = (e0)
e_exact = (e2 + e0)
v_exact = deepcopy(v0)
v_model = deepcopy(v_guess)
v_diabat = deepcopy(ci_vector_dimer)
FermiCG.eye!(v_diabat)

# First make bare 
println()
println(" Make Model Space Bare Hamiltonian")
Heff_bare = make_H(v_diabat, v_bare, e_bare)
for i in 1:length(e_bare)
    @printf(" Adiabatic: %12.8f Diabatic: %12.8f\n", e_bare[i], Heff_bare[i,i])
end
Heff_bare .*= conversion

# Now exact
println()
println(" Make Model Space Effective Hamiltonian")
Heff_exact = make_H(v_diabat, v_exact, e_exact)
for i in 1:length(e_exact)
    @printf(" Adiabatic: %12.8f Diabatic: %12.8f\n", e_exact[i], Heff_exact[i,i])
end
Heff_exact .*= conversion


Heff_bare = Heff_bare - diagm(diag(Heff_bare))
Heff_exact = Heff_exact - diagm(diag(Heff_exact))

Heff_bare = abs.(Heff_bare - diagm(diag(Heff_bare)))
Heff_exact = abs.(Heff_exact - diagm(diag(Heff_exact)))

@load "../tetracene_tetramer/thresh_0004.jld2"
@load "../tetracene_tetramer/H_guess_0004.jld2"
@load "../tetracene_tetramer/v_model_tetramer.jld2"

display(FermiCG.get_vector(v_guess))
v_barett = v_guess
e_barett = e_guess
e_exacttt = (e0)
e_exact = (e2 + e0)
v_exacttt = deepcopy(v0)
v_modeltt = deepcopy(v_guess)
v_diabattt = deepcopy(ci_vector)
FermiCG.eye!(v_diabattt)

# First make bare 
println()
println(" Make Model Space Bare Hamiltonian")
Heff_bare_tt = make_H(v_diabattt, v_barett, e_barett)
Heff_bare_tt .*= conversion
diab = [1,5,3,9,7,13,11,18,30,24]
Heff_bare_tt = Heff_bare_tt[diab, diab]

# Now exact
println()
println(" Make Model Space Effective Hamiltonian")
Heff_exact_tt = make_H(v_diabattt, v_exacttt, e_exacttt)
Heff_exact_tt .*= conversion
diab = [1,5,3,9,7,13,11,18,30,24]
Heff_exact_tt = Heff_exact_tt[diab, diab]


Heff_bare_tt = Heff_bare_tt - diagm(diag(Heff_bare_tt))
Heff_exact_tt = Heff_exact_tt - diagm(diag(Heff_exact_tt))

Heff_bare_tt = abs.(Heff_bare_tt - diagm(diag(Heff_bare_tt)))
Heff_exact_tt = abs.(Heff_exact_tt - diagm(diag(Heff_exact_tt)))
max_val = max(maximum(abs.(Heff_bare)), maximum(abs.(Heff_exact)), maximum(abs.(Heff_bare_tt)), maximum(abs.(Heff_exact_tt)))
max_val = maximum(abs.(Heff_bare))

#plotd = heatmap(Heff*conversion; color=palette([:teal, :white, :orange], 10), aspect_ratio=1, dpi=300, size=(300,300), right_margin = 10Plots.mm)
#plotd = heatmap(Heff_bare; color=palette(:RdGy_10, 100), aspect_ratio=1, dpi=300, size=(300,300), right_margin = 10Plots.mm,  
plotd = heatmap(Heff_bare; color=palette(:bone_1, rev=true), aspect_ratio=1, dpi=300, size=(300,300), right_margin = 10Plots.mm,  
                clims=(0, max_val), ticks = false,xaxis=false,yaxis=false, 
                xlims = (0.5,10.5), ylims = (0.5,10.5), title="Heff_bare")
m, n = size(Heff_bare)
vline!(0.5:(n+0.5), c=:white, label=false)
hline!(0.5:(m+0.5), c=:white, label=false)

savefig(plotd,"Heff_bare.png")


max_val = maximum(abs.(Heff_exact))

plotd = heatmap(Heff_exact; color=palette(:bone_1, rev=true), aspect_ratio=1, dpi=300, size=(300,300), right_margin = 10Plots.mm,  
                clims=(0, max_val), ticks = false, xaxis=false,yaxis=false,
                xlims = (0.5,10.5), ylims = (0.5,10.5), title="Heff_exact")
#plotd = heatmap(Heff*27.2114079527*1000; color=palette(:Blues, 10), aspect_ratio=1, dpi=300, size=(300,300), right_margin = 10Plots.mm, xlims = (0,5), ylims = (0,5))
m, n = size(Heff_bare)
vline!(0.5:(n+0.5), c=:white, label=false)
hline!(0.5:(m+0.5), c=:white, label=false)
savefig(plotd,"Heff_exact.png")

max_val = maximum(abs.(Heff_bare_tt))
#plotd = heatmap(Heff*conversion; color=palette([:teal, :white, :orange], 10), aspect_ratio=1, dpi=300, size=(300,300), right_margin = 10Plots.mm)
plotd = heatmap(Heff_bare_tt; color=palette(:bone_1, rev=true), aspect_ratio=1, dpi=300, size=(300,300), right_margin = 10Plots.mm,  
                clims=(0, max_val), ticks = false,xaxis=false,yaxis=false, 
                xlims = (0.5,10.5), ylims = (0.5,10.5), title="Heff_bare from TT")
m, n = size(Heff_bare)
vline!(0.5:(n+0.5), c=:white, label=false)
hline!(0.5:(m+0.5), c=:white, label=false)


savefig(plotd,"Heff_bare_tt.png")



max_val = maximum(abs.(Heff_exact_tt))
plotd = heatmap(Heff_exact_tt; color=palette(:bone_1, rev=true), aspect_ratio=1, dpi=300, size=(300,300), right_margin = 10Plots.mm,  
                clims=(0, max_val), ticks = false, xaxis=false,yaxis=false,
                xlims = (0.5,10.5), ylims = (0.5,10.5), title="Heff_exact from TT")
#plotd = heatmap(Heff*27.2114079527*1000; color=palette(:Blues, 10), aspect_ratio=1, dpi=300, size=(300,300), right_margin = 10Plots.mm, xlims = (0,5), ylims = (0,5))
m, n = size(Heff_bare)
vline!(0.5:(n+0.5), c=:white, label=false)
hline!(0.5:(m+0.5), c=:white, label=false)
savefig(plotd,"Heff_exact_tt.png")

