from read_operator import read_operator,get_expectation_value
from read_circuit import read_circuit

qc = read_circuit('../../second_quantization/circuits/hardware_efficient/cascade_3/bh_2.3_circuit.txt')
B  = read_operator('../../second_quantization/operators/bh_2.3_h.txt')
E  = get_expectation_value(B,qc)

print("E ",E)

def wg(qc):
    d  = qc.depth()
    o  = qc.count_ops()
    g  = sum([o[k] for k in o.keys()])
    g2 = o['cx']
    return d,(g-g2,g2)

for molecule in ['bh','hf','beh2','h2o']:
    for ansatz in ['ry_linear_1','ry_full_1','cascade_1']:
        qc = read_circuit('../../second_quantization/circuits/hardware_efficient/%s/%s_2.3_circuit.txt' % (ansatz,molecule))
        print(ansatz,molecule,wg(qc))
    qc = read_circuit('../../second_quantization/circuits/qUCCSD/unrestricted/singles_doubles/trotter/%s_reps_1_1.3_circuit.txt' % molecule)
    print('quccsd_reps_1',molecule,wg(qc))

print("*"*53)

for molecule in ['bh','hf','beh2','h2o']:
    qc = read_circuit('../../first_quantization/trim/circuits/cascade_5/%s_2.5_circuit.txt' % (molecule))
    print(molecule,wg(qc))

print("*"*53)

for molecule in ['bh','hf','beh2','h2o']:
    for ansatz in ['ry_linear_5','cascade_5']:
        qc = read_circuit('../../first_quantization/pad/variation_after_projection/circuits/%s/%s_2.5_circuit.txt' % (ansatz,molecule))
        print(ansatz,molecule,wg(qc))

