from qiskit                                        import *
from qiskit.chemistry.drivers                      import UnitsType,HFMethodType
from qiskit.chemistry.core                         import Hamiltonian,TransformationType,QubitMappingType
from qiskit.chemistry.components.initial_states    import HartreeFock
from qiskit.aqua.components.optimizers             import L_BFGS_B,COBYLA,CG
from qiskit.aqua.algorithms                        import VQE
from qiskit.aqua                                   import QuantumInstance,aqua_globals
from qiskit.circuit.library                        import EfficientSU2

import logging
from   qiskit.chemistry import set_qiskit_chemistry_logging
from   qiskit.aqua      import set_qiskit_aqua_logging
set_qiskit_chemistry_logging(logging.INFO)
set_qiskit_aqua_logging(logging.INFO)

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

import sys
import numpy as np
sys.path.append('/Users/mario/Documents/GitHub/QITE/qite_es/scf_calculations/subroutines')
sys.path.append('/Users/mario/Documents/GitHub/QITE/qite_es/scf_calculations/qUCCSD')
sys.path.append('/Users/mario/Documents/GitHub/QITE/qite_es/scf_calculations/subroutines/pyscfd')
from pyscfdriver   import *
if('unrestricted'=='restricted_closed_shell'):
    from restricted_closed_shell_uccsd import UCCSD
else:
    from unrestricted_uccsd import UCCSD
from utils         import *
from CustomVarForm import *

dist    = 3.3
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

init_state = HartreeFock(num_orbitals=core._molecule_info['num_orbitals'],qubit_mapping=core._qubit_mapping,two_qubit_reduction=core._two_qubit_reduction,
                         num_particles=core._molecule_info['num_particles'],sq_list=sqlist)

outfile.write("\nHartree-Fock energy %f \n" % (molecule.hf_energy))
outfile.write("\nHartree-Fock circuit\n")
outfile.write(str(init_state.construct_circuit().draw())+"\n")

# -----------------------------------------------------------------------------------------------------------------------

rev = ('singles_doubles'=='doubles_singles')

var_form  = UCCSD(num_orbitals=core._molecule_info['num_orbitals'],num_particles=core._molecule_info['num_particles'],
                  active_occupied=None,active_unoccupied=None,initial_state=init_state,qubit_mapping=core._qubit_mapping,
                  two_qubit_reduction=core._two_qubit_reduction,num_time_slices=1,z2_symmetries=z2syms,
                  reverse_excitations=rev,expansion_mode='suzuki',reps=1)
optimizer = COBYLA(maxiter=0)
depth = 1
if(depth==1): p0 = np.loadtxt('/Users/mario/Documents/GitHub/QITE/qite_es/scf_calculations/qUCCSD/quccsd_flavors_h2o/unrestricted/suzuki/ds_%s/R_%s/output_parameters.txt' % (str(rev),str(dist)))
else:         p0 = np.loadtxt('/Users/mario/Documents/GitHub/QITE/qite_es/scf_calculations/qUCCSD/quccsd_flavors_h2o_reps_2/unrestricted/suzuki/ds_%s/R_%s/output_parameters.txt' % (str(rev),str(dist)))
algo      = VQE(H_op,var_form,optimizer,aux_operators=A_op,include_custom=True,initial_point=p0)
backend          = Aer.get_backend('statevector_simulator')
quantum_instance = QuantumInstance(backend=backend)
algo_result      = algo.run(quantum_instance)

# -----------------------------------------------------------------------------------------------------------------------

p1 = algo._ret['opt_params']
np.savetxt("output_parameters.txt",p1)

res_vqe,res_ee = get_results(H_op,A_op,molecule,core,algo_result,outfile)
np.save('results.npy',{'hf_results':molecule.hf_energy,'num_qubits':H_op.num_qubits,'num_parameters':var_form.num_parameters,'results_vqe':res_vqe,'results_eigensolver':res_ee},allow_pickle=True)

outfile.write("\nCircuit\n")
outfile.write(str(algo.get_optimal_circuit().draw())+"\n")

Stot = A_op[1]
for Ei,Oi in zip(var_form.excitations,var_form._hopping_ops):
    comm = (Stot*Oi-Oi*Stot).chop(1e-10)
    print("Excitation ",Ei," commutes with total spin ",comm.is_empty())

p1 = algo._ret['opt_params']
res_vqe,res_ee = get_results(H_op,A_op,molecule,core,algo_result,outfile)
depth = 1
print_results('unrestricted/singles_doubles/suzuki/h2o_reps_%d_%s_results.txt' % (depth,str(dist)),res_vqe,res_ee,dN=4)
circuit = algo.get_optimal_circuit()
print_circuit('unrestricted/singles_doubles/suzuki/h2o_reps_%d_%s_circuit.txt' % (depth,str(dist)),circuit)
