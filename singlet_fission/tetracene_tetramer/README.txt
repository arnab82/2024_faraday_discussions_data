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
thresh fois = 4e-5
thresh cipsi = 0.0004
thresh spin = 0.0008
pt2 thresh foi = 1e-8

# Bare Hamiltonian
vguess obtained from tps ci direct
e2guess thresh foi = 1e-8
Hbare = VEV'

# JLD2 Files
H_guess_spin.jld2: e_guess, v_guess, e2_guess, ecore, Hbare
thresh_spin_0.0004.jld2: clusters, e0, v0, e2, s2, ecore, init_fspace
correlations.jld2: correlation dictonary with S2


