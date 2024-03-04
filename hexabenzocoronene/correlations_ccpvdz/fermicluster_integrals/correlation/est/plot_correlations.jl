using JLD2
using FermiCG
using Plots
using LinearAlgebra
using Printf
#the clusters are in order of 1,2,6,3,5,7,4 , the matices are reordered to 1,2,3,4,5,6,7 to plot n2, sz2, q2

function run()
    @load("tucker_thresh_5_1e4.jld2")
    @load("../../hexabenzocoronene_old.jld2")
    cmf_state = TPSCIstate(clusters, FockConfig(init_fspace))


    display(v0b)
    cf = correlation_functions(v0b, cmf_state)
    @save "correlations_ccpvdzcorrelation.jld2" cf
    n1 = cf["N"][1]
    n2 = cf["N"][2]
    sz1 = cf["Sz"][1]
    sz2 = cf["Sz"][2]
    q1 = cf["Q"][1]
    q2 = cf["Q"][2]
    
    # display(n2[1])
    n2=[0.02550829	-0.00759149	-0.00017893	-0.0004038	-0.00018108	-0.00759755	-0.00955544;
    -0.00759149	0.02551065	-0.00759694	-0.00018174	-0.00040329	-0.0001787	-0.0095585;
    -0.00017893	-0.00759694	0.02548365	-0.00759733	-0.00017889	-0.00040443	-0.00952713;
    -0.0004038	-0.00018174	-0.00759733	0.02550731	-0.00759024	-0.00017861	-0.00955559;
    -0.00018108	-0.00040329	-0.00017889	-0.00759024	0.02551021	-0.0075972	-0.0095595;
    -0.00759755	-0.0001787	-0.00040443	-0.00017861	-0.0075972	0.02548714	-0.00953065;
    -0.00955544	-0.0095585	-0.00952713	-0.00955559	-0.0095595	-0.00953065	0.05728681]
    max_val = max(0, maximum(abs.(n2)))

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
    display(plotd)

    savefig(plotd,@sprintf("n_correlation_hbc.png"))
    
    sz2=[0.00925446	-0.0027587	-0.00004766	-0.00013415	-0.00004725	-0.00276328	-0.00350342;
    -0.0027587	0.00925935	-0.00276279	-0.00004727	-0.00013475	-0.00004721	-0.00350864;
    -0.00013415	-0.00004727	0.00925467	-0.00276338	-0.00004716	-0.00013419	-0.0034995;
    -0.00004766	-0.00276279	-0.00276338	0.00925471	-0.00275859	-0.00004761	-0.00350371;
    -0.00004725	-0.00013475	-0.00004716	-0.00275859	0.00926062	-0.002763	-0.00350988;
    -0.00276328	-0.00004721	-0.00013419	-0.00004761	-0.002763	0.00925669	-0.00350141;
    -0.00350342	-0.00350864	-0.0034995	-0.00350371	-0.00350988	-0.00350141	0.02102656]
    max_val1 = max(0, maximum(abs.(sz2)))
    plotd = heatmap(sz2; color=palette(:RdGy_9, 100), aspect_ratio=1, dpi=300, size=(300,300), right_margin = 10Plots.mm,  
                        clims=(-max_val1, max_val1), ticks = false,xaxis=false,yaxis=false, 
                        xlims = (0.5,7.5), ylims = (0.5,7.5),yflip=true)
    m, n = size(n2)
    vline!(0.5:(n+0.5), c=:grey, label=false)
    hline!(0.5:(m+0.5), c=:grey, label=false)


    savefig(plotd,@sprintf("sz_correlation_hbc.png"))
    
    q2=[0.03307346	0.01079417	0.0001041	-0.00005672	0.00011049	0.01080873	0.01342415;
    0.01079417	0.03307551	0.01080557	0.00011064	-0.00005735	0.00010328	0.01342571;
    0.0001041	0.01080557	0.03307124	0.01080811	0.00010314	-0.00006024	0.01342283;
    -0.00005672	0.00011064	0.01080811	0.03307314	0.01079245	0.00010397	0.01342431;
    0.00011049	-0.00005735	0.00010314	0.01079245	0.03307654	0.01080523	0.01342799;
    0.01080873	0.00010328	-0.00006024	0.00010397	0.01080523	0.03307683	0.01342283;
    0.01342415	0.01342571	0.01342283	0.01342431	0.01342799	0.01342283	0.06866612]
    max_val1 = max(0, maximum(abs.(q2)))
    plotd = heatmap(q2; color=palette(:RdGy_9, 100), aspect_ratio=1, dpi=300, size=(300,300), right_margin = 10Plots.mm,  
                        clims=(-max_val1, max_val1), ticks = false,xaxis=false,yaxis=false, 
                        xlims = (0.5,7.5), ylims = (0.5,7.5),yflip=true)
    m, n = size(q2)
    vline!(0.5:(n+0.5), c=:grey, label=false)
    hline!(0.5:(m+0.5), c=:grey, label=false)


    savefig(plotd,@sprintf("q_correlation_hbc.png"))
    
    s2=[0.03086426	0.00917567	-0.00020937	-0.00016823	-0.00021138	0.00913799	0.0116173;
    0.00917567	0.03081745	0.00913992	-0.000211	-0.00017273	-0.00021103	0.01157484;
    -0.00020937	0.00913992	0.03099207	0.00913734	-0.00021131	-0.00015027	0.01175149;
    -0.00016823	-0.000211	0.00913734	0.03086289	0.00917414	-0.00020901	0.01161784;
    -0.00021138	-0.00017273	-0.00021131	0.00917414	0.03082064	0.00914194	0.01157813;
    0.00913799	-0.00021103	-0.00015027	-0.00020901	0.00914194	0.03099828	0.01175363;
    0.0116173	0.01157484	0.01175149	0.01161784	0.01157813	0.01175363	0.07135078]
    max_val1 = max(0, maximum(abs.(s2)))
    plotd = heatmap(s2; color=palette(:RdGy_9, 100), aspect_ratio=1, dpi=300, size=(300,300), right_margin = 10Plots.mm,  
                        clims=(-max_val1, max_val1), ticks = false,xaxis=false,yaxis=false, 
                        xlims = (0.5,7.5), ylims = (0.5,7.5),yflip=true)
    m, n = size(s2)
    vline!(0.5:(n+0.5), c=:grey, label=false)
    hline!(0.5:(m+0.5), c=:grey, label=false)


    savefig(plotd,@sprintf("s2_correlation_hbc.png"))
    

end

run()