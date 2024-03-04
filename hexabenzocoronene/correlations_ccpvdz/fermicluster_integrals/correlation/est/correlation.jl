function run()
    @load "/Users/arnab/arnab/workspace/2024_faraday_discussions_data/hexabenzocoronene/correlations_ccpvdz/fermicluster_integrals/correlation/est/tucker_thresh_5_2e4.jld2"
    @load "/Users/arnab/arnab/workspace/2024_faraday_discussions_data/hexabenzocoronene/correlations_ccpvdz/fermicluster_integrals/hexabenzocoronene_old.jld2"
    cluster_ops = FermiCG.compute_cluster_ops(cluster_bases, ints);

    FermiCG.add_cmf_operators!(cluster_ops, cluster_bases, ints, d1.a, d1.b);
    display(v0b)
    C_f = correlation_functions(v0b, cluster_ops)
    @save "/Users/arnab/arnab/workspace/2024_faraday_discussions_data/hexabenzocoronene/correlations_ccpvdz/fermicluster_integrals/correlation/est/correlation_S2_H.jld2" C_f
    
    # max_val = max(0, maximum(abs.(S2_covar)))

    # plotd = heatmap(S2_covar; color=palette(:RdGy_9, 100), aspect_ratio=1, dpi=300, size=(300,300), right_margin = 10Plots.mm,  
    #                     clims=(-max_val, max_val), ticks = false,xaxis=false,yaxis=false, 
    #                     xlims = (0.5,7.5), ylims = (0.5,7.5),yflip=true)
    # m, n = size(S2_covar)
    # print(m)    
    # print(n)
    # vline!(0.5:(n+0.5), c=:grey, label=false)
    # hline!(0.5:(m+0.5), c=:grey, label=false)
    # display(plotd)

    # savefig(plotd,@sprintf("S2_correlation_hbc.png"))
    
end

run()
