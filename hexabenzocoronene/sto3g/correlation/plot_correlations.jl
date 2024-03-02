using JLD2
using FermiCG
using Plots
using LinearAlgebra
using Printf


function run()
    @load("tucker_thresh_5_2e4.jld2")
    @load("hexabenzocoronene_sto3g.jld2")
    cmf_state = TPSCIstate(clusters, FockConfig(init_fspace))


    display(v0b)
    cf = correlation_functions(v0b, cmf_state)
    @save "correlation.jld2" cf
    n1 = cf["N"][1]
    n2 = cf["N"][2]
    sz1 = cf["Sz"][1]
    sz2 = cf["Sz"][2]
    q1 = cf["Q"][1]
    q2 = cf["Q"][2]
    
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
    display(plotd)

    savefig(plotd,@sprintf("n_correlation_hbc.png"))
    

    max_val1 = max(0, maximum(abs.(sz2[1])))
    plotd = heatmap(sz2[1]; color=palette(:RdGy_9, 100), aspect_ratio=1, dpi=300, size=(300,300), right_margin = 10Plots.mm,  
                        clims=(-max_val1, max_val1), ticks = false,xaxis=false,yaxis=false, 
                        xlims = (0.5,7.5), ylims = (0.5,7.5),yflip=true)
    m, n = size(n2[1])
    vline!(0.5:(n+0.5), c=:grey, label=false)
    hline!(0.5:(m+0.5), c=:grey, label=false)


    savefig(plotd,@sprintf("sz_correlation_hbc.png"))
    

    max_val1 = max(0, maximum(abs.(q2[1])))
    plotd = heatmap(q2[1]; color=palette(:RdGy_9, 100), aspect_ratio=1, dpi=300, size=(300,300), right_margin = 10Plots.mm,  
                        clims=(-max_val1, max_val1), ticks = false,xaxis=false,yaxis=false, 
                        xlims = (0.5,7.5), ylims = (0.5,7.5),yflip=true)
    m, n = size(q2[1])
    vline!(0.5:(n+0.5), c=:grey, label=false)
    hline!(0.5:(m+0.5), c=:grey, label=false)


    savefig(plotd,@sprintf("q_correlation_hbc.png"))


end

run()