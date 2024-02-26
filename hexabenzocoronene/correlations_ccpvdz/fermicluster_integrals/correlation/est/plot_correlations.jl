using JLD2
using FermiCG
using Plots
using LinearAlgebra
using Printf
using NPZ

function run()
    @load("/Users/arnab/arnab/workspace/2024_faraday_discussions_data/hexabenzocoronene/correlations_ccpvdz/fermicluster_integrals/correlation/est/tucker_thresh_4e4.jld2")


    display(v0b)
    n2=npzread("/Users/arnab/arnab/workspace/2024_faraday_discussions_data/hexabenzocoronene/correlations_ccpvdz/fermicluster_integrals/correlation/est/N2_correlation.npy")
    sz2=npzread("/Users/arnab/arnab/workspace/2024_faraday_discussions_data/hexabenzocoronene/correlations_ccpvdz/fermicluster_integrals/correlation/est/sz2_correlation.npy")
    display(n2)
    max_val = max(0, maximum(abs.(n2[1])))

    plotd = heatmap(n2; color=palette(:RdGy_9, 100), aspect_ratio=1, dpi=300, size=(300,300), right_margin = 10Plots.mm,  
                        clims=(-max_val, max_val), ticks = false,xaxis=false,yaxis=false, 
                        xlims = (0.5,7.5), ylims = (0.5,7.5),yflip=true)
    m, n = size(n2)
    print(m)    
    print(n)
        #vline!(0.5:(n+0.5), c=:white, label=false)
        #hline!(0.5:(m+0.5), c=:white, label=false)
    vline!(0.5:(n+0.5), c=:grey, label=false)
    hline!(0.5:(m+0.5), c=:grey, label=false)


    savefig(plotd,@sprintf("/Users/arnab/arnab/workspace/2024_faraday_discussions_data/hexabenzocoronene/correlations_ccpvdz/fermicluster_integrals/correlation/est/n_correlation_hbc.png"))
    max_val1 = max(0, maximum(abs.(sz2[1])))
    plotd = heatmap(sz2; color=palette(:RdGy_9, 100), aspect_ratio=1, dpi=300, size=(300,300), right_margin = 10Plots.mm,  
                        clims=(-max_val1, max_val1), ticks = false,xaxis=false,yaxis=false, 
                        xlims = (0.5,7.5), ylims = (0.5,7.5),yflip=true)
    m, n = size(n2)
    vline!(0.5:(n+0.5), c=:grey, label=false)
    hline!(0.5:(m+0.5), c=:grey, label=false)


    savefig(plotd,@sprintf("/Users/arnab/arnab/workspace/2024_faraday_discussions_data/hexabenzocoronene/correlations_ccpvdz/fermicluster_integrals/correlation/est/sz_correlation_hbc.png"))

end

run()

