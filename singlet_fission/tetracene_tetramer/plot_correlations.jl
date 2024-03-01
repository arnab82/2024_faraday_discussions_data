using JLD2
using FermiCG
using Plots
using LinearAlgebra
using Printf


function run()

    @load("thresh_spin_0.0004.jld2")
    cmf = FermiCG.TPSCIstate(clusters, FermiCG.FockConfig(init_fspace), R=nroots);
    nroots = 31
    
    display(v0)
    n1, n2, sz1, sz2 = correlation_functions(v0)

    cf = correlation_functions(v0, cmf)

    nclusters = length(v0.clusters)

    r_wanted = [1,6,7,10,12,19,21,24,25,28,31]
    
    max_val = 0    
    #for r in 1:nroots
    for r in r_wanted
        max_val = max(max_val, maximum(abs.(n2[r])))
    end

    for r in r_wanted
    #for r in 1:nroots
        ordering = collect(1:4)
        n2r = n2[r][ordering,:][:,ordering]
        display(n2r)
        #plotd = heatmap(n2r; color=palette(:RdGy_9, 100), aspect_ratio=1, dpi=300, size=(250,200), right_margin = 15Plots.mm,  
        #                clims=(-max_val, max_val),ticks = false,xaxis=false,yaxis=false, 
        #                xlims = (0.5,nclusters+.5), ylims = (0.5,nclusters+.5),
        #                yflip=true, title="State "*string(r), legend=true)
        plotd = heatmap(n2r; color=palette(:RdGy_9, 100), aspect_ratio=1, dpi=300, size=(200,200), right_margin = Plots.mm,  
                        clims=(-max_val, max_val), ticks = false,xaxis=false,yaxis=false, 
                        xlims = (0.5,nclusters+.5), ylims = (0.5,nclusters+.5),
                        yflip=true, title="State "*string(r), legend=false)
        m, n = size(n2r)

        vline!(0.5:(n+0.5), c=:grey, label=false)
        hline!(0.5:(m+0.5), c=:grey, label=false)


        #savefig(plotd,@sprintf("plots/n_1_label_%1i.png", r))
        #break
        savefig(plotd,@sprintf("plots/n_correlation_%1i.png", r))
    end
    

    max_val = 0    
    #for r in 1:nroots
    for r in r_wanted
        max_val = max(max_val, maximum(abs.(sz2[r])))
    end

    #for r in 1:nroots
    for r in r_wanted
        ordering = collect(1:4)
        sz2r = sz2[r][ordering,:][:,ordering]
        display(sz2r)
        
        #plotd = heatmap(sz2r; color=palette(:RdGy_9, 100), aspect_ratio=1, dpi=300, size=(250,200), right_margin = 12Plots.mm,  
        #                clims=(-max_val, max_val), ticks = false,xaxis=false,yaxis=false, 
        #                xlims = (0.5,nclusters+.5), ylims = (0.5,nclusters+.5),
        #                yflip=true, title="State "*string(r), legend=true)
        plotd = heatmap(sz2r; color=palette(:RdGy_9, 100), aspect_ratio=1, dpi=300, size=(200,200), right_margin = Plots.mm,  
                        clims=(-max_val, max_val), ticks = false,xaxis=false,yaxis=false, 
                        xlims = (0.5,nclusters+.5), ylims = (0.5,nclusters+.5),
                        yflip=true,title="State "*string(r), legend=false)
        m, n = size(sz2r)
        vline!(0.5:(n+0.5), c=:grey, label=false)
        hline!(0.5:(m+0.5), c=:grey, label=false)


        #savefig(plotd,@sprintf("plots/sz_1_label_%1i.png", r))
        #break
        savefig(plotd,@sprintf("plots/sz_correlation_%1i.png", r))
    end
    
    max_val = 0    
    #for r in 1:nroots
    for r in r_wanted
        max_val = max(max_val, maximum(abs.(cf["Q"][2][r])))
    end

    #for r in 1:nroots
    for r in r_wanted
        ordering = collect(1:4)
        cfr = cf["Q"][2][r][ordering,:][:,ordering]
        display(cfr)
        #plotd = heatmap(cfr; color=palette(:RdGy_9, 100), aspect_ratio=1, dpi=300, size=(250,200), right_margin = 12Plots.mm,  
        #                clims=(-max_val, max_val), ticks = false,xaxis=false,yaxis=false, 
        #                xlims = (0.5,nclusters+.5), ylims = (0.5,nclusters+.5),
        #                yflip=true, title="State "*string(r), legend=true)
        plotd = heatmap(cfr; color=palette(:RdGy_9, 100), aspect_ratio=1, dpi=300, size=(200,200), right_margin = Plots.mm,  
                        clims=(-max_val, max_val), ticks = false,xaxis=false,yaxis=false, 
                        xlims = (0.5,nclusters+.5), ylims = (0.5,nclusters+.5),
                        yflip=true, title="State "*string(r), legend=false)
        m, n = size(cfr)
        vline!(0.5:(n+0.5), c=:grey, label=false)
        hline!(0.5:(m+0.5), c=:grey, label=false)


        #savefig(plotd,@sprintf("plots/q_1_label_%1i.png", r))
        #break
        savefig(plotd,@sprintf("plots/q_correlation_%1i.png", r))
    end
    
end

run()


