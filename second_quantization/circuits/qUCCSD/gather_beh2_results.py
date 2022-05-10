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
from pyscfdriver   import *
from utils         import *
from CustomVarForm import *

def print_circuit(fname,circuit):
    def to_str(x):
        if(type(x)==np.float or type(x)==np.float64):
           return str(x)
        elif(type(x)==np.int):
           return str(float(x))
        else:
           return str(x._symbol_expr)
    outf = open(fname,'w')
    for gate in circuit:
        gate_name = gate[0].name
        qubits = [q.index for q in gate[1]]
        outf.write(str(gate_name)+' | '+' '.join([str(x) for x in qubits])+' | '+' '.join([to_str(x) for x in gate[0].params])+'\n')
    outf.close()

def print_results(fname,res_vqe,res_ee,dN=0):
    outf = open(fname,'w')
    outf.write('%.12f %.12f \n' % (res_vqe[0],np.real(res_ee[0])))
    outf.write('%.12f %.12f \n' % (res_vqe[1]+dN,np.real(res_ee[1]+dN)))
    outf.write('%.12f %.12f \n' % (res_vqe[2],np.real(res_ee[2])))
    outf.write('%.12f %.12f \n' % (res_vqe[3],np.real(res_ee[3])))
    outf.close()

dist_list = [0.7,0.9,1.1,1.3,1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.1,3.3,3.5,3.7,3.9,4.1,4.3,4.5]
depth_list = [1,2,3,4,5,6,7,8]

for depth in [1,2]:
  for dist in dist_list:
    outfile = open('pes.txt','w')
    
    outfile.write("Bond distance: {}\n".format(dist))
    rho = np.load('/Users/mario/Documents/GitHub/QITE/qite_es/scf_calculations/SCF_FCI/beh2/res_fourth_pass.npy',allow_pickle=True).item()
    rho = rho[str(dist)]['rho_scf']
    driver    = PySCFDriver(atom="Be 0.0000 0.0000 0.0000; H 0.0000 0.0000 -{0}; H 0.0000 0.0000 {0}".format(dist),
                            unit=UnitsType.ANGSTROM,charge=0,spin=0,basis='sto-6g',hf_method=HFMethodType.RHF,symgroup='dooh',
                            outfile=outfile,irrep_nelec={'A1g':4,'A1u':2},rho=rho)
    molecule  = driver.run()
    outfile.write("\nHartree-Fock energy %f \n" % (molecule.hf_energy))
    
    core      = Hamiltonian(transformation=TransformationType.FULL,qubit_mapping=QubitMappingType.PARITY,two_qubit_reduction=True,freeze_core=False,orbital_reduction=[0])
    H_op,A_op = core.run(molecule)
    
    z2syms,sqlist           = None,None
    H_op,A_op,z2syms,sqlist = taper(molecule,core,H_op,A_op,outfile)
    
    init_state = HartreeFock(num_orbitals=core._molecule_info['num_orbitals'],qubit_mapping=core._qubit_mapping,two_qubit_reduction=core._two_qubit_reduction,num_particles=core._molecule_info['num_particles'],sq_list=sqlist)
    
    outfile.write("\nHartree-Fock energy %f \n" % (molecule.hf_energy))
    outfile.write("\nHartree-Fock circuit\n")
    outfile.write(str(init_state.construct_circuit().draw())+"\n")
    
    # -----------------------------------------------------------------------------------------------------------------------

    if(depth==1):
       var_form  = UCCSD(num_orbitals=core._molecule_info['num_orbitals'],num_particles=core._molecule_info['num_particles'],active_occupied=None,active_unoccupied=None,
                         initial_state=init_state,qubit_mapping=core._qubit_mapping,two_qubit_reduction=core._two_qubit_reduction,num_time_slices=1,z2_symmetries=z2syms)
       p0 = np.loadtxt('/Users/mario/Documents/GitHub/QITE/qite_es/scf_calculations/qUCCSD/beh2/R_%s/output_parameters.txt' % str(dist))
    else:
       var_form  = UCCSD(reps=2,num_orbitals=core._molecule_info['num_orbitals'],num_particles=core._molecule_info['num_particles'],active_occupied=None,active_unoccupied=None,
                         initial_state=init_state,qubit_mapping=core._qubit_mapping,two_qubit_reduction=core._two_qubit_reduction,num_time_slices=2,z2_symmetries=z2syms)
       p0 = np.loadtxt('/Users/mario/Documents/GitHub/QITE/qite_es/scf_calculations/qUCCSD/beh2_reps2/R_%s/output_parameters.txt' % str(dist))

    optimizer = COBYLA(maxiter=0)
    algo      = VQE(H_op,var_form,optimizer,aux_operators=A_op,include_custom=True,initial_point=p0)
   
    backend          = Aer.get_backend('statevector_simulator')
    quantum_instance = QuantumInstance(backend=backend)
    algo_result      = algo.run(quantum_instance)
   
    # -----------------------------------------------------------------------------------------------------------------------
    
    p1 = algo._ret['opt_params']
    res_vqe,res_ee = get_results(H_op,A_op,molecule,core,algo_result,outfile)
    print_results('unrestricted/singles_doubles/trotter/beh2_reps_%d_%s_results.txt' % (depth,str(dist)),res_vqe,res_ee,dN=2)

    circuit = algo.get_optimal_circuit().decompose()
    print_circuit('unrestricted/singles_doubles/trotter/beh2_reps_%d_%s_circuit.txt' % (depth,str(dist)),circuit)
    
