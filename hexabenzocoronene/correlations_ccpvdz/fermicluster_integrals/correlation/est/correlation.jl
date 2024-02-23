using FermiCG
using JLD2
using NPZ
@load "tucker_thresh_4e4.jld2"
@time n1,n2,sz1,sz2=FermiCG.correlation_functions(v0b)
@save "correlation_cc_pvdz_fermicluster_4e4.jld2" n1 n2 sz1 sz2
