# Tetracene Tetramer

structure: same as excited state paper
basis set: 6-31g* 
integrals: generated using the cis-no method with PM localization (not psuedo-can), 4T and 4S averaged into the rdm

# TPSCI
no HOSVD bootstrapping for these results
M = 150
delta elec = 5 for each chromophore
cluster basis = cluster eigenbasis spin
nroots = 31
thresh fois = 1e-5
thresh cipsi = 0.0004
pt2 thresh foi = 1e-8

# Bare Hamiltonian
vguess obtained from tps ci direct
e2guess thresh foi = 1e-8
Hbare = VEV'

# JLD2 Files
Hguess0004.jld2: e_guess, v_guess, e2_guess, ecore, Hbare
thresh0004.jld2: clusters, cluster_bases, e0, v0, e2, s2, ecore


