using JLD2
using FermiCG
using Plots
using LinearAlgebra
using Printf


function run()
    @load("/Users/arnab/arnab/workspace/2024_faraday_discussions_data/hexabenzocoronene/correlations_ccpvdz/orbital_partioning_integrals/correlation/spin/tucker_thresh_1e3.jld2")


    display(v0b)
    n1, n2, sz1, sz2 = correlation_functions(v0b)
    display(n2[1])
    max_val = max(0, maximum(abs.(n2[1])))

    plotd = heatmap(n2[1]; color=palette(:RdGy_9, 100), aspect_ratio=1, dpi=300, size=(300,300), right_margin = 10Plots.mm,  
                        clims=(-max_val, max_val), ticks = false,xaxis=false,yaxis=false, 
                        xlims = (0.5,7.5), ylims = (0.5,7.5),yflip=true)
    m, n = size(n2[1])
    print(m)    
    print(n)
        #vline!(0.5:(n+0.5), c=:white, label=false)
        #hline!(0.5:(m+0.5), c=:white, label=false)
    vline!(0.5:(n+0.5), c=:grey, label=false)
    hline!(0.5:(m+0.5), c=:grey, label=false)


    savefig(plotd,@sprintf("/Users/arnab/arnab/workspace/2024_faraday_discussions_data/hexabenzocoronene/correlations_ccpvdz/orbital_partioning_integrals/correlation/spin/n_correlation_hbc.png"))
    max_val1 = max(0, maximum(abs.(sz2[1])))
    plotd = heatmap(sz2[1]; color=palette(:RdGy_9, 100), aspect_ratio=1, dpi=300, size=(300,300), right_margin = 10Plots.mm,  
                        clims=(-max_val1, max_val1), ticks = false,xaxis=false,yaxis=false, 
                        xlims = (0.5,7.5), ylims = (0.5,7.5),yflip=true)
    m, n = size(n2[1])
    vline!(0.5:(n+0.5), c=:grey, label=false)
    hline!(0.5:(m+0.5), c=:grey, label=false)


    savefig(plotd,@sprintf("/Users/arnab/arnab/workspace/2024_faraday_discussions_data/hexabenzocoronene/correlations_ccpvdz/orbital_partioning_integrals/correlation/spin/sz_correlation_hbc.png"))

end

run()

