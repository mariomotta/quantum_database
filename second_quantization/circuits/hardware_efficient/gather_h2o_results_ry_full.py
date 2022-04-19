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

import sys
import numpy as np
sys.path.append('/Users/mario/Documents/GitHub/QITE/qite_es/scf_calculations/subroutines')
sys.path.append('/Users/mario/Documents/GitHub/QITE/qite_es/scf_calculations/subroutines/pyscfd')
from pyscfdriver import *
from utils       import *

def print_circuit(fname,circuit):
    outf = open(fname,'w')
    for gate in circuit:
        gate_name = gate[0].name
        qubits = [q.index for q in gate[1]]
        outf.write(str(gate_name)+' | '+' '.join([str(x) for x in qubits])+' | '+' '.join([str(x._symbol_expr) for x in gate[0].params])+'\n')
    outf.close()

def print_results(fname,res_vqe,res_ee,dN=0):
    outf = open(fname,'w')
    outf.write('%.12f %.12f \n' % (res_vqe[0],np.real(res_ee[0])))
    outf.write('%.12f %.12f \n' % (res_vqe[1]+dN,np.real(res_ee[1]+dN)))
    outf.write('%.12f %.12f \n' % (res_vqe[2],np.real(res_ee[2])))
    outf.write('%.12f %.12f \n' % (res_vqe[3],np.real(res_ee[3])))
    outf.close()

dist_list = [0.7,0.9,1.1,1.3,1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.1,3.3]
depth_list = [2] #[1,2,3,4,5]

for depth in depth_list:
  for dist in dist_list:

    outfile = open('pes.txt','w')
    
    outfile.write("Bond distance: {}\n".format(dist))
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
    
    init_state = HartreeFock(num_orbitals=core._molecule_info['num_orbitals'],qubit_mapping=core._qubit_mapping,two_qubit_reduction=core._two_qubit_reduction,num_particles=core._molecule_info['num_particles'],sq_list=sqlist)
    init_state.construct_circuit()
    x_hf = init_state._bitstr
 
    outfile.write("\nHartree-Fock energy %f \n" % (molecule.hf_energy))
    outfile.write("\nHartree-Fock circuit\n")
    outfile.write(str(init_state.construct_circuit().draw())+"\n")

    par = np.loadtxt('/Users/mario/Documents/GitHub/QITE/qite_es/scf_calculations/Ry_full_DEPTH_%d/h2o/R_%s/output_parameters.txt' % (depth,str(dist)))
    var_form  = EfficientSU2(num_qubits=H_op.num_qubits,reps=depth,entanglement='full',su2_gates=['ry'],initial_state=init_state)
    circuit = var_form.bind_parameters(par)

    optimizer = COBYLA(maxiter=0)
    algo      = VQE(H_op,var_form,optimizer,aux_operators=A_op,include_custom=True,initial_point=par)
    backend          = Aer.get_backend('statevector_simulator')
    quantum_instance = QuantumInstance(backend=backend)
    algo_result      = algo.run(quantum_instance)
    
    # -----------------------------------------------------------------------------------------------------------------------
    
    p1 = algo._ret['opt_params']
    res_vqe,res_ee = get_results(H_op,A_op,molecule,core,algo_result,outfile)
    print_results('ry_full_%d/h2o_%s_results.txt' % (depth,str(dist)),res_vqe,res_ee,dN=4)

    circuit = var_form.bind_parameters(p1)
    print_circuit('ry_full_%d/h2o_%s_circuit.txt' % (depth,str(dist)),circuit)

