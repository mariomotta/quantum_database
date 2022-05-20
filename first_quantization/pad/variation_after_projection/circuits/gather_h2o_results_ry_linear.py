from pyscf import gto,scf
import numpy as np
import sys
sys.path.append('/Users/mario/Documents/GitHub/QITE/qite_es/first_quantization_pad/src/')
from subroutines import get_hamiltonian_matrix,get_spin_square_matrix,project_along_irrep,project_along_singlet,pad_matrix,trim,get_orbital_info,permute_orbitals
from quantum import matrix_to_qubit_operator,make_ry_ansatz,do_VQE
from VAP_vqe import VQE

import numpy as np
sys.path.append('/Users/mario/Documents/GitHub/QITE/qite_es/scf_calculations/subroutines/')
from utils         import *

from qiskit import Aer
from qiskit.chemistry.drivers                      import UnitsType,HFMethodType
from qiskit.chemistry.core                         import Hamiltonian,TransformationType,QubitMappingType
from qiskit.chemistry.components.initial_states    import HartreeFock
from qiskit.chemistry.components.variational_forms import UCCSD
from qiskit.aqua.components.optimizers             import L_BFGS_B,COBYLA,CG
from qiskit.aqua                                   import QuantumInstance,aqua_globals
from qiskit.circuit.library                        import EfficientSU2

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
    outf.write('%.12f %.12f \n' % (res_vqe[1],np.real(res_ee[1])))
    outf.write('%.12f %.12f \n' % (res_vqe[2],np.real(res_ee[2])))
    outf.close()

def geometry(d):
    z0 = 0.1177
    z1 = 0.47116
    y0 = 0.75545
    r0 = np.sqrt(y0**2+(z0+z1)**2)
    r  = d/r0
    return "O 0.0000 0.0000 {0}; H 0.0000 {1} -{2}; H 0.0000 -{1} -{2}".format(z0*r,y0*r,z1*r)

dist_list = [0.7,0.9,1.1,1.3,1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.1,3.3]
depth_list = [6,7] #[3,4,5]

for depth in depth_list:
    for dist in dist_list:
        rho  = np.load('/Users/mario/Documents/GitHub/QITE/qite_es/scf_calculations/SCF_FCI/h2o/res_fourth_pass.npy',allow_pickle=True).item()
        rho  = rho[str(dist)]['rho_scf']
        mol = gto.Mole()
        mol.build(atom  = geometry(dist),
                  basis = 'sto-6g',symmetry = 'c2v',spin = 0,charge = 0,verbose = 4)
        mf       = scf.RHF(mol)
        mf       = scf.newton(mf)
        mf.kernel(rho)

        orb_info = get_orbital_info(mf)
        perm = permute_orbitals(orb_info,nf=1,sort_occupied=['B2','B1','A1','A1'],sort_virtual=['A1','B1'])
        H_fci,mycas = get_hamiltonian_matrix(mol,mf,nf=1,permutation=perm)
        S_fci       = get_spin_square_matrix(mol,mf,nf=1,nci=H_fci.shape[0])
        print("original CI Hamiltonian matrix size ..........",H_fci.shape[0])
        
        H_fci,S_fci = project_along_irrep(mol,mf,mycas,nf=1,H_fci=H_fci,S_fci=S_fci)
        print("after projection on GS irrep, matrix size ....",H_fci.shape[0])
        
        H_fci,S_fci = project_along_singlet(H_fci,S_fci)
        print("after projection on singlet, matrix size .....",H_fci.shape[0])
        
        H_fci,S_fci,Proj,nqubit = pad_matrix(H_fci,S_fci)
        print("after padding, matrix size ...................",H_fci.shape[0])
        
        H_op = matrix_to_qubit_operator(H_fci)
        P = matrix_to_qubit_operator(Proj)
        J = matrix_to_qubit_operator(np.dot(Proj,np.dot(H_fci,Proj)))
        A_op = [P,J]

        var_form = make_ry_ansatz(nqubit,layers=depth)
        p0 = np.loadtxt('/Users/mario/Documents/GitHub/QITE/qite_es/first_quantization_pad_VAP_Ry/h2o/cascade_depth_%d/R_%s/output_parameters.txt'%(depth,str(dist)))
        optimizer = COBYLA(maxiter=0)
        algo      = VQE(H_op,var_form,optimizer,aux_operators=A_op,include_custom=True,initial_point=p0)
        backend          = Aer.get_backend('statevector_simulator')
        quantum_instance = QuantumInstance(backend=backend)
        algo_result      = algo.run(quantum_instance)
       
        p1 = algo._ret['opt_params']
        outfile = open('ry_linear_results.txt','w')
        res_vqe,res_ee = get_results_first_quantization(H_op,A_op,algo_result,outfile)
        print_results('ry_linear_%d/h2o_%s_results.txt' % (depth,str(dist)),res_vqe,res_ee)
        circuit = var_form.bind_parameters(p1)
        print_circuit('ry_linear_%d/h2o_%s_circuit.txt' % (depth,str(dist)),circuit)
 
