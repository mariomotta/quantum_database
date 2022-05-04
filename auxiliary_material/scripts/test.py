from read_operator import read_operator,get_expectation_value
from read_circuit import read_circuit

qc = read_circuit('../../second_quantization/circuits/hardware_efficient/cascade_3/bh_2.3_circuit.txt')
B  = read_operator('../../second_quantization/operators/bh_2.3_h.txt')
E  = get_expectation_value(B,qc)

print("E ",E)
