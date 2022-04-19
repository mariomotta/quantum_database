from qiskit                                        import *
from qiskit.chemistry.drivers                      import UnitsType,HFMethodType
from qiskit.chemistry.core                         import Hamiltonian,TransformationType,QubitMappingType
from qiskit.chemistry.components.initial_states    import HartreeFock
from qiskit.chemistry.components.variational_forms import UCCSD
from qiskit.aqua.components.optimizers             import L_BFGS_B,COBYLA,CG
from qiskit.aqua.algorithms                        import VQE
from qiskit.aqua                                   import QuantumInstance,aqua_globals
from qiskit.circuit.library                        import EfficientSU2

import logging
from   qiskit.chemistry import set_qiskit_chemistry_logging
from   qiskit.aqua      import set_qiskit_aqua_logging
set_qiskit_chemistry_logging(logging.INFO)
set_qiskit_aqua_logging(logging.INFO)
from qiskit.aqua.algorithms                        import NumPyEigensolver

def write_to_file(X,dX,fname):
    outf = open(fname,'w')
    for c,P in X._paulis:
        P = P.to_label()
        if(P=='I'*X.num_qubits):
           c += dX
        outf.write('%s %.12f\n' % (P,c.real)) 
    outf.close()


import sys
import numpy as np
sys.path.append('/Users/mario/Documents/GitHub/QITE/qite_es/scf_calculations//subroutines')
sys.path.append('/Users/mario/Documents/GitHub/QITE/qite_es/scf_calculations//subroutines/pyscfd')
from pyscfdriver   import *
from utils         import *
from CustomVarForm import *

dist_list = [0.7,0.9,1.1,1.3,1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.1,3.3]

for dist in dist_list:
    outfile = open('pes.txt','w')
    rho = np.load('/Users/mario/Documents/GitHub/QITE/qite_es/scf_calculations/SCF_FCI/h2o/res_fourth_pass.npy',allow_pickle=True).item()
    rho = rho[str(dist)]['rho_scf']
    z0 = 0.1177
    z1 = 0.47116
    y0 = 0.75545
    r0 = np.sqrt(y0**2 + (z0+z1)**2)
    ratio = dist/r0
    driver    = PySCFDriver(atom="O 0.0000 0.0000 {0}; H 0.0000 {1} -{2}; H 0.0000 -{1} -{2}".format(z0*ratio, y0*ratio, z1*ratio),
                            unit=UnitsType.ANGSTROM,charge=0,spin=0,basis='sto-6g',hf_method=HFMethodType.RHF,symgroup='c2v',
                            outfile=outfile,irrep_nelec={'A1':6,'A2':0,'B1':2,'B2':2},rho=rho)
    molecule  = driver.run()
    outfile.write("\nHartree-Fock energy %f \n" % (molecule.hf_energy))
    
    core      = Hamiltonian(transformation=TransformationType.FULL,qubit_mapping=QubitMappingType.PARITY,two_qubit_reduction=True,freeze_core=False,orbital_reduction=[0,1])
    H_op,A_op = core.run(molecule)
    
    z2syms,sqlist           = None,None
    H_op,A_op,z2syms,sqlist = taper(molecule,core,H_op,A_op,outfile)

    offsets = [core._energy_shift + core._ph_energy_shift + core._nuclear_repulsion_energy,4,0,0]

    write_to_file(H_op,offsets[0],'h2o_%s_h.txt'%dist)
    write_to_file(A_op[0],offsets[1],'h2o_%s_ne.txt'%dist)
    write_to_file(A_op[1],offsets[2],'h2o_%s_s2.txt'%dist)
    write_to_file(A_op[2],offsets[3],'h2o_%s_sz.txt'%dist)

