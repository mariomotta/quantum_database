def read_circuit(filename='',num_qubits=None):
    from qiskit import QuantumCircuit
    f = open(filename,'r')
    f = [fi.split('|') for fi in f.readlines()]
    names = [ fi[0]                           for fi in f]
    qbits = [ [int(x) for x in fi[1].split()] for fi in f]
    if(num_qubits is None):
       num_qubits = max([max(qi) for qi in qbits])+1
    qc = QuantumCircuit(num_qubits)
    for i,fi in enumerate(f):
        name = names[i]
        qbit = qbits[i]
        name = name[:len(name)-1]
        import difflib
        if(name=='cx'):   qc.cx(qbit[0],qbit[1])
        elif(name=='x'):  qc.x(qbit[0])
        elif(name=='ry'): qc.ry(float(fi[2]),qbit[0])
        elif(name=='p'):  qc.p(float(fi[2]),qbit[0])
        elif(name=='u3'):
             fi2 = fi[2].split()
             qc.u3(float(fi2[0]),float(fi2[1]),float(fi2[2]),qbit[0])
        elif(name=='u'): 
             fi2 = fi[2].split()
             qc.u(float(fi2[0]),float(fi2[1]),float(fi2[2]),qbit[0])
        elif(name=='h'):  qc.h(qbit[0])
        else:
          exit()
    return qc

