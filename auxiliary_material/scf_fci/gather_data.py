import numpy as np
import matplotlib.pyplot as plt
import sys,os
sys.path.append('/Users/mario/Documents/GitHub/QITE/qite_es/scf_calculations/SCF_FCI')

import pyscf
from   pyscf            import gto,scf,mcscf,fci
from   scipy            import linalg as LA
from   orbital_analysis import *

L        = 5.0
fig,ax   = plt.subplots(4,2,figsize=(2*L,0.5*4*L))
fig.subplots_adjust(left=0.2)
c_list = {'purple'       : '#B284BE',
          'jacaranda'    : '#888FC7',
          'light_yellow' : '#EEEA62',
          'gold'         : '#FFD300',
          'earwax'       : '#DAA520',
          'brown'        : '#7B3F00',
          'light_blue'   : '#bcd9ea',
          'cerulean'     : '#5ba4cf',
          'cobalt'       : '#0079bf',
          'dark_blue'    : '#055a8c',
          'light_green'  : '#acdf87',
          'yellow_green' : '#B9D146',
          'mid_green'    : '#68bb59',
          'dark_green'   : '#1e5631',
          'orange'       : '#FF8856',
          'red'          : '#DC343B',
          'light-gray'   : '#C0C0C0',
          'teal'         : '#008080',
          'palatinate'   : '#72246C'}

style = {'A1'  :{'color':c_list['red'],         's':'--','w':2,'marker':'o','ms':16},
         'E1x' :{'color':c_list['cobalt'],      's':'--','w':2,'marker':'P','ms':16},
         'E1y' :{'color':c_list['dark_blue'],   's':'--','w':2,'marker':'X','ms':16},
         'A1g' :{'color':c_list['orange'],      's':'-.','w':2,'marker':'>','ms':13},
         'A1u' :{'color':c_list['earwax'],      's':'-.','w':2,'marker':'<','ms':13},
         'E1ux':{'color':c_list['yellow_green'],'s':'-.','w':2,'marker':'v','ms':11},
         'E1uy':{'color':c_list['mid_green'],   's':'-.','w':2,'marker':'^','ms':11},
         'B1'  :{'color':c_list['purple'],      's':'--','w':2,'marker':'s','ms':16},
         'B2'  :{'color':c_list['jacaranda'],   's':'--','w':2,'marker':'D','ms':16}}

def fill_panel(pan,xlabel,xlim,xticks,xticklabels,ylabel,ylim,yticks,yticklabels,p=20.0,q=10.0):
    x0,x1 = xlim
    xlim  = [x0-(x1-x0)/p,x1+(x1-x0)/p]
    pan.set_xlabel(xlabel)
    pan.set_xlim(xlim)
    pan.set_xticks(xticks)
    pan.set_xticklabels(xticklabels)
    pan.set_ylabel(ylabel)
    y0,y1 = ylim
    ylim  = [y0-(y1-y0)/q,y1+(y1-y0)/q]
    pan.set_ylim(ylim)
    pan.set_yticks(yticks)
    pan.set_yticklabels(yticklabels)
    pan.tick_params(direction='in',which='both')

# ------------------------------------------------------------- 

def perform_computations(folder,make_geometry,dist_list,group,irrep_nelec):
    res    = np.load('%s/res_fourth_pass.npy'%folder,allow_pickle=True).item()
    folder = folder.split('/')
    folder = folder[len(folder)-1]
    mo_eig = []
    no_occ = []
    mo_irr = []
    no_irr = []
    E_scf  = np.zeros((len(dist_list),2))
    for i_dist,dist in enumerate(dist_list):
        mol = gto.Mole()
        mol.build(atom     = make_geometry(dist),
                  basis    = 'sto-6g',
                  symmetry = group,
                  verbose  = 3)
        mf             = scf.RHF(mol)
        mf.irrep_nelec = irrep_nelec
        mf             = scf.newton(mf)
        mf.max_cycle   = 1
        if(folder=='h2o' and dist>2 and dist<2.69):
          Ehf            = mf.kernel(res['2.7']['rho_scf'])
        else:
          Ehf            = mf.kernel(res[str(dist)]['rho_scf'])
        a              = mf.stability()[0]
        Ehf            = mf.kernel(a,mf.mo_occ)
        Cmo            = mf.mo_coeff
        mo_eig.append(mf.mo_energy)
        mo_irr.append(get_irreps(mol,Cmo))
      
        cisolver  = fci.FCI(mol,Cmo)
        e, fcivec = cisolver.kernel()
        dm1       = cisolver.make_rdm1(fcivec,mol.nao_nr(),(mol.nelectron//2,mol.nelectron//2))
        # -------------------------------------------------------------------------------------
        IRR       = get_irreps(mol,Cmo)
        IRR_SET   = set(IRR)
        n         = Cmo.shape[1]
        s         = np.zeros(n)
        Cno       = np.zeros((n,n))
        for i in IRR_SET:
            idx        = [k for k in range(n) if IRR[k]==i]
            dm1_i      = dm1[idx,:].copy()
            dm1_i      = dm1_i[:,idx]
            si,Ui      = LA.eigh(dm1_i)
            s[idx]     = si
            Cno[:,idx] = np.dot(Cmo[:,idx],Ui)
        for i in range(len(s)):
            for j in range(i+1,len(s)):
                if(s[j]>s[i]):
                   s[j],s[i] = s[i],s[j]
                   cj,ci = Cno[:,j].copy(),Cno[:,i].copy()
                   Cno[:,j],Cno[:,i] = ci,cj
        print(" >>> ",s,get_irreps(mol,Cno))
        no_occ.append(s)
        no_irr.append(get_irreps(mol,Cno))
        R = str(dist)
        E_scf[i_dist,0] = dist
        E_scf[i_dist,1] = Ehf
    return mo_eig,mo_irr,no_occ,no_irr,E_scf

# ------------------------------------------------------------- top: bh

def geometry_bh(d):
    return "B 0.0000 0.0000 0.0000; H 0.0000 0.0000 {0}".format(d)
dist_list = [0.7,0.9,1.1,1.3,1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.1,3.3,3.5,3.7,3.9,4.1,4.3,4.5,4.7,4.9,5.1,5.3]
mo_eig,mo_irr,no_occ,no_irr,E_scf = perform_computations('/Users/mario/Documents/GitHub/QITE/qite_es/scf_calculations/SCF_FCI/bh',geometry_bh,dist_list,'coov',{'A1':6})

outf = open('bh_scf_data.txt','w')
for R,mo_eig_R,mo_irr_R,no_occ_R,no_irr_R in zip(dist_list,mo_eig,mo_irr,no_occ,no_irr):
    outf.write('%s\n'%R)
    for ei,irri,nj,irrj in zip(mo_eig_R,mo_irr_R,no_occ_R,no_irr_R):
        outf.write('%.12f %s %.12f %s\n'%(ei,irri,nj,irrj))
outf.close()
np.savetxt('bh_scf_energy.txt',E_scf)

def geometry_hf(d):
    return "F 0.0000 0.0000 0.0000; H 0.0000 0.0000 {0}".format(d)
dist_list = [0.7,0.9,1.1,1.3,1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.1,3.3,3.5,3.7,3.9,4.1,4.3,4.5]
mo_eig,mo_irr,no_occ,no_irr,E_scf = perform_computations('/Users/mario/Documents/GitHub/QITE/qite_es/scf_calculations/SCF_FCI/hf',geometry_hf,dist_list,'coov',{'A1':6,'E1x':2,'E1y':2})

outf = open('hf_scf_data.txt','w')
for R,mo_eig_R,mo_irr_R,no_occ_R,no_irr_R in zip(dist_list,mo_eig,mo_irr,no_occ,no_irr):
    outf.write('%s\n'%R)
    for ei,irri,nj,irrj in zip(mo_eig_R,mo_irr_R,no_occ_R,no_irr_R):
        outf.write('%.12f %s %.12f %s\n'%(ei,irri,nj,irrj))
outf.close()
np.savetxt('hf_scf_energy.txt',E_scf)

# ------------------------------------------------------------- third: beh2

def geometry_beh2(d):
    return "Be 0.0000 0.0000 0.0000; H 0.0000 0.0000 -{1}; H 0.0000 0.0000 {1}".format(d,d)
dist_list = [0.7,0.9,1.1,1.3,1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.1,3.3,3.5,3.7,3.9,4.1,4.3,4.5]
mo_eig,mo_irr,no_occ,no_irr,E_scf = perform_computations('/Users/mario/Documents/GitHub/QITE/qite_es/scf_calculations/SCF_FCI/beh2',geometry_beh2,dist_list,'dooh',{'A1g':4,'A1u':2})

outf = open('beh2_scf_data.txt','w')
for R,mo_eig_R,mo_irr_R,no_occ_R,no_irr_R in zip(dist_list,mo_eig,mo_irr,no_occ,no_irr):
    outf.write('%s\n'%R)
    for ei,irri,nj,irrj in zip(mo_eig_R,mo_irr_R,no_occ_R,no_irr_R):
        outf.write('%.12f %s %.12f %s\n'%(ei,irri,nj,irrj))
outf.close()
np.savetxt('beh2_scf_energy.txt',E_scf)

# ------------------------------------------------------------- fourth: h2o

def geometry_h2o(d):
    z0 = 0.1177
    z1 = 0.47116
    y0 = 0.75545
    r0 = np.sqrt(y0**2+(z0+z1)**2)
    r  = d/r0
    return "O 0.0000 0.0000 {0}; H 0.0000 {1} -{2}; H 0.0000 -{1} -{2}".format(z0*r,y0*r,z1*r)
dist_list = [0.7,0.9,1.1,1.3,1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.1,3.3]
mo_eig,mo_irr,no_occ,no_irr,E_scf = perform_computations('/Users/mario/Documents/GitHub/QITE/qite_es/scf_calculations/SCF_FCI/h2o',geometry_h2o,dist_list,'c2v',{'A1':6,'B1':2,'B2':2})

outf = open('h2o_scf_data.txt','w')
for R,mo_eig_R,mo_irr_R,no_occ_R,no_irr_R in zip(dist_list,mo_eig,mo_irr,no_occ,no_irr):
    outf.write('%s\n'%R)
    for ei,irri,nj,irrj in zip(mo_eig_R,mo_irr_R,no_occ_R,no_irr_R):
        outf.write('%.12f %s %.12f %s\n'%(ei,irri,nj,irrj))
outf.close()
np.savetxt('h2o_scf_energy.txt',E_scf)


