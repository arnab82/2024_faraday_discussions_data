# Tetracene Dimer (2,4) in tetramer same as (1,3)
This dimer has the parallel chromophores

# Mapping from tetracene model space
diab = [1,5,3,9,7,13,11,18,30,24]

structure: same as excited state paper
basis set: 6-31g*
integrals: generated using the cis-no method with Boys localization (not psuedo-can), 2T and 2S averaged into the rdm
NOTE: this parallel dimer was optimized with Boys instead of PM bec PM was not localizing orbs onto chromophores (i.e. mos 1 and 2 were delocalized across the dimer)
active space: (10o,10e) per cluster so a total of (20o,20e)

# CMF
cmf was optimized with cmf diis
rdm was initalized with overlap and density matrices
diis_start=3
convg_tol = default (1e-6)

# TPSCI
no HOSVD bootstrapping for these results
M = 150
delta elec = 5 for each chromophore
cluster basis = cluster eigenbasis spin
nroots = 10
thresh fois = 1e-5
thresh cipsi = 0.0004
pt2 thresh foi = 1e-8

# Bare Hamiltonian
vguess obtained from tps ci direct
e2guess thresh foi = 1e-8
Hbare = PVEPV'

# JLD2 Files
Hguess0004.jld2: e_guess, v_guess, e2_guess, ecore, Hbare
thresh0004.jld2: clusters, cluster_bases, e0, v0, e2, s2, ecore

# Generate plots
using the "form Heff nick" julia script to generate the plots for the Hbare, Heff, and the associated Hbare and Heff extracted
for the tetramer data

# Nuclear Energy
3356.0296680809283


