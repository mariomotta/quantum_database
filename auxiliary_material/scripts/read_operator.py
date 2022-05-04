def read_operator(filename):
    from qiskit.opflow.primitive_ops import PauliSumOp
    f = [fi.split() for fi in open(filename,'r').readlines()]
    coeff = [float(fi[1]) for fi in f]
    pauli = [fi[0] for fi in f]
    pauli_op = [([pauli,weight]) for pauli,weight in zip(pauli,coeff)]
    H = PauliSumOp.from_list([op for op in pauli_op])
    return H

def get_expectation_value(operator,state_circuit):
    import numpy as np
    from qiskit.quantum_info import Statevector
    energy = np.real(Statevector(state_circuit).data@operator.to_matrix()@Statevector(state_circuit).data)
    return energy

