
import numpy as np
import pyscf
from pyscf import cc



basis = "aug-cc-pvdz"

## Loop over geometries
#for geom in range(40,41):
for geom in range(1,41):

    filename = "geom_%03i.xyz"%(geom)
    mol = pyscf.gto.M(atom=filename)

    mol.basis = basis
    mol.build()
    
    mf = pyscf.scf.RHF(mol)
    mf.run()

    mycc = cc.CCSD(mf).run()
    print('CCSD total energy', mycc.e_tot)
    et = mycc.ccsd_t()
    print('CCSD(T) total energy', mycc.e_tot + et)


    e_ee, c_ee = mycc.eeccsd(nroots=15)
    
    print(" EOM-CCSD geom: %03i %12.8f "%(geom, mycc.e_tot), end='')
    for i in e_ee:
        print("%12.8f "%(i+mycc.e_tot), end='')
    print()

