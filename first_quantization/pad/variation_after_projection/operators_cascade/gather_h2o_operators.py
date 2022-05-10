from pyscf import gto,scf
import numpy as np
import sys
sys.path.append('/Users/mario/Documents/GitHub/QITE/qite_es/first_quantization_pad/src/')
from subroutines import get_hamiltonian_matrix,get_spin_square_matrix,project_along_irrep,project_along_singlet,pad_matrix,trim
from quantum import matrix_to_qubit_operator,make_cascade_ansatz,do_VQE

import numpy as np
sys.path.append('/Users/mario/Documents/GitHub/QITE/qite_es/scf_calculations/subroutines/')
from scipy import linalg as LA
from utils import *

def geometry(d):
    z0 = 0.1177
    z1 = 0.47116
    y0 = 0.75545
    r0 = np.sqrt(y0**2+(z0+z1)**2)
    r  = d/r0
    return "O 0.0000 0.0000 {0}; H 0.0000 {1} -{2}; H 0.0000 -{1} -{2}".format(z0*r,y0*r,z1*r)

def write_to_file(X,dX,fname):
    outf = open(fname,'w')
    for c,P in X._paulis:
        P = P.to_label()
        if(P=='I'*X.num_qubits):
           c += dX
        outf.write('%s %.12f\n' % (P,c.real))
    outf.close()

dist_list = [0.7,0.9,1.1,1.3,1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.1,3.3]

info = []
for dist in dist_list:
    rho  = np.load('/Users/mario/Documents/GitHub/QITE/qite_es/scf_calculations/SCF_FCI/h2o/res_fourth_pass.npy',allow_pickle=True).item()
    rho  = rho[str(dist)]['rho_scf']
    
    mol = gto.Mole()
    mol.build(atom     = geometry(dist),
              basis    = 'sto-6g',
              symmetry = 'c2v',
              spin     = 0,
              charge   = 0,
              verbose  = 4)
    mf       = scf.RHF(mol)
    mf       = scf.newton(mf)
    mf.kernel(rho)
    
    H_fci,mycas = get_hamiltonian_matrix(mol,mf,nf=1)
    S_fci       = get_spin_square_matrix(mol,mf,nf=1,nci=H_fci.shape[0])
    d1 = H_fci.shape[0]
    e1 = LA.eigh(H_fci)[0][0]

    H_fci,S_fci = project_along_irrep(mol,mf,mycas,nf=1,H_fci=H_fci,S_fci=S_fci)
    d2 = H_fci.shape[0]
    
    H_fci,S_fci = project_along_singlet(H_fci,S_fci)
    d3 = H_fci.shape[0]

    H_fci,S_fci,Proj,nqubit = pad_matrix(H_fci,S_fci)
    print("after padding, matrix size ...................",H_fci.shape[0])
    d4 = H_fci.shape[0]
    e4 = LA.eigh(H_fci)[0][0]

    H_op = matrix_to_qubit_operator(H_fci)
    H_op = matrix_to_qubit_operator(H_fci)
    P = matrix_to_qubit_operator(Proj)
    J = matrix_to_qubit_operator(np.dot(Proj,np.dot(H_fci,Proj)))

    write_to_file(H_op,0,'h2o_%s_h.txt'%dist)
    write_to_file(P,0,'h2o_%s_p.txt'%dist)
    write_to_file(J,0,'h2o_%s_php.txt'%dist)

    info.append([d1,d2,d3,d4,e1,e4])

outf = open('h2o_info.txt','w')
for d1,d2,d3,d4,e1,e4 in info:
    outf.write('%d %d %d %d %.12f %.12f \n' % (d1,d2,d3,d4,e1,e4))
outf.close()
 
