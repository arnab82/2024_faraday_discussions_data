using JLD2
using FermiCG
using Plots
using LinearAlgebra
using Printf


function run(cf)

    nroots = length(cf["N"][1])

    ordering = [1,2,3,4,5,6,7]

    tmp = cf["N"][2]
    max_val = 0    
    for r in 1:nroots
        max_val = max(max_val, maximum(abs.(tmp[r])))
    end
    for r in 1:nroots
        tmpr = cf["N"][2][r]
        tmpr = tmpr[ordering,:][:,ordering]
        plotd = heatmap(tmpr; color=palette(:RdGy_9, 100), aspect_ratio=1, dpi=300, size=(300,300), right_margin = 10Plots.mm,  
                        clims=(-max_val, max_val), ticks = false,xaxis=false,yaxis=false, 
                        xlims = (0.5,5.5), ylims = (0.5,5.5),
                        yflip=true)
        m, n = size(tmpr)
        #vline!(0.5:(n+0.5), c=:white, label=false)
        #hline!(0.5:(m+0.5), c=:white, label=false)
        vline!(0.5:(n+0.5), c=:grey, label=false)
        hline!(0.5:(m+0.5), c=:grey, label=false)


        savefig(plotd,@sprintf("n_correlation_%1i.png", r))
    end



    tmp = cf["Sz"][2]
    max_val = 0    
    for r in 1:nroots
        max_val = max(max_val, maximum(abs.(tmp[r])))
    end
    for r in 1:nroots
        tmpr = tmp[r]
        tmpr = tmpr[ordering,:][:,ordering]
        plotd = heatmap(tmpr; color=palette(:RdGy_9, 100), aspect_ratio=1, dpi=300, size=(300,300), right_margin = 10Plots.mm,  
                        clims=(-max_val, max_val), ticks = false,xaxis=false,yaxis=false, 
                        xlims = (0.5,5.5), ylims = (0.5,5.5),
                        yflip=true)
        m, n = size(tmpr)
        #vline!(0.5:(n+0.5), c=:white, label=false)
        #hline!(0.5:(m+0.5), c=:white, label=false)
        vline!(0.5:(n+0.5), c=:grey, label=false)
        hline!(0.5:(m+0.5), c=:grey, label=false)


        savefig(plotd,@sprintf("sz_correlation_%1i.png", r))
    end



    tmp = cf["S2"][2]
    max_val = 0    
    for r in 1:nroots
        max_val = max(max_val, maximum(abs.(tmp[r])))
    end
    for r in 1:nroots
        tmpr = tmp[r]
        tmpr = tmpr[ordering,:][:,ordering]
        plotd = heatmap(tmpr; color=palette(:RdGy_9, 100), aspect_ratio=1, dpi=300, size=(300,300), right_margin = 10Plots.mm,  
                        clims=(-max_val, max_val), ticks = false,xaxis=false,yaxis=false, 
                        xlims = (0.5,5.5), ylims = (0.5,5.5),
                        yflip=true)
        m, n = size(tmpr)
        #vline!(0.5:(n+0.5), c=:white, label=false)
        #hline!(0.5:(m+0.5), c=:white, label=false)
        vline!(0.5:(n+0.5), c=:grey, label=false)
        hline!(0.5:(m+0.5), c=:grey, label=false)


        savefig(plotd,@sprintf("s2_correlation_%1i.png", r))
    end
end

run(cf)


