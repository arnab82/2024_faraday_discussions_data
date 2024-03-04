using JLD2
using FermiCG
using Plots
using LinearAlgebra
using Printf


function run()
    @load "tucker_thresh_5_2e4.jld2"
    @load "hexabenzocoronene_sto3g.jld2"
    # cluster_ops = FermiCG.compute_cluster_ops(cluster_bases, ints);

    # FermiCG.add_cmf_operators!(cluster_ops, cluster_bases, ints, d1.a, d1.b);
    # display(v0b)
    # C_f = correlation_functions(v0b, cluster_ops)
    # @save "correlation_S2_H.jld2" C_f
    S2_covar= [0.05121371   0.01513706  -0.00055531  -0.00062580  -0.00058424   0.01505366   0.01856175;
    0.01513706   0.05122004   0.01506346  -0.00058314  -0.00062630  -0.00055687   0.01855274;
   -0.00055531   0.01506346   0.05151813   0.01504831  -0.00055639  -0.00055967   0.01891469;
   -0.00062580  -0.00058314   0.01504831   0.05120134   0.01513773  -0.00055599   0.01855204;
   -0.00058424  -0.00062630  -0.00055639   0.01513773   0.05122470   0.01507440   0.01854542;
    0.01505366  -0.00055687  -0.00055967  -0.00055599   0.01507440   0.05153637   0.01891449;
    0.01856175   0.01855274   0.01891469   0.01855204   0.01854542   0.01891449   0.11646180]
    max_val = max(0, maximum(abs.(S2_covar)))

    plotd = heatmap(S2_covar; color=palette(:RdGy_9, 100), aspect_ratio=1, dpi=300, size=(300,300), right_margin = 10Plots.mm,  
                        clims=(-max_val, max_val), ticks = false,xaxis=false,yaxis=false, 
                        xlims = (0.5,7.5), ylims = (0.5,7.5),yflip=true)
    m, n = size(S2_covar)
    print(m)    
    print(n)
    vline!(0.5:(n+0.5), c=:grey, label=false)
    hline!(0.5:(m+0.5), c=:grey, label=false)
    display(plotd)

    savefig(plotd,@sprintf("S2_correlation_hbc.png"))
    
end

run()



