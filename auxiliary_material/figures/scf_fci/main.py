import numpy as np
import os
import matplotlib.pyplot as plt
import itertools
from matplotlib import rc

rc('font',**{'family':'serif','serif':['cmu serif'],'size':12})
rc('text', usetex=True)

L        = 5.0
fig,ax   = plt.subplots(4,2,figsize=(2*L,0.5*4*L))
fig.subplots_adjust(left=0.2)

c_list = {'yellow'      : '#DAA520',
          'orange'      : '#FF8856',
          'red'         : '#DC343B',
          'magenta'     : '#B284BE',
          'purple'      : '#888FC7',
          'blue'        : '#0079BF',
          'cyan'        : '#5BA4CF',
          'light-green' : '#B9D146',
          'dark-green'  : '#32826E'}

style = {'A1'  :{'color':c_list['red'],        's':'-', 'w':2,'marker':'o','ms':16,'label':'A$_1$'},
         'E1x' :{'color':c_list['blue'],       's':'--','w':2,'marker':'P','ms':16,'label':'E$_{1x}$'},
         'E1y' :{'color':c_list['cyan'],       's':'--','w':2,'marker':'X','ms':16,'label':'E$_{1y}$'},
         'A1g' :{'color':c_list['light-green'],'s':'-', 'w':2,'marker':'>','ms':13,'label':'A$_{1g}$'},
         'A1u' :{'color':c_list['dark-green'], 's':'--','w':2,'marker':'<','ms':13,'label':'A$_{1u}$'},
         'E1ux':{'color':c_list['magenta'],    's':'-.','w':2,'marker':'v','ms':11,'label':'E$_{1ux}$'},
         'E1uy':{'color':c_list['purple'],     's':'-.','w':2,'marker':'^','ms':11,'label':'E$_{1uy}$'},
         'B1'  :{'color':c_list['orange'],     's':'-.','w':2,'marker':'s','ms':16,'label':'B$_1$'},
         'B2'  :{'color':c_list['yellow'],     's':':', 'w':2,'marker':'D','ms':16,'label':'B$_1$'}}

def read_scf_fci_data(fname):
    inf = open(fname,'r')
    inf = [x.split() for x in inf.readlines()]
    data = {'geometries':[],
            'molecular_orbital_energies':[],
            'molecular_orbital_irreps':[],
            'natural_orbital_occupancies':[],
            'natural_orbital_irreps':[]}
    ind = []
    for j,r in enumerate(inf):
        if(len(r)==1):
           data['geometries'].append(float(r[0]))
           ind.append(j)
    data['geometries'] = np.array(data['geometries'])
    nR = len(data['geometries'])
    no = ind[1]-ind[0]-1
    for j in ind:
        data['molecular_orbital_energies'].append([None]*no)
        data['molecular_orbital_irreps'].append([None]*no)
        data['natural_orbital_occupancies'].append([None]*no)
        data['natural_orbital_irreps'].append([None]*no)
    for l,j in enumerate(ind):
        for k in range(no):
            moe,moi,noo,noi = inf[j+1+k]
            data['molecular_orbital_energies'][l][k] = float(moe)
            data['molecular_orbital_irreps'][l][k] = moi
            data['natural_orbital_occupancies'][l][k] = float(noo)
            data['natural_orbital_irreps'][l][k] = noi
    for l,j in enumerate(ind):
        data['molecular_orbital_energies'][l] = np.array(data['molecular_orbital_energies'][l])
        data['natural_orbital_occupancies'][l] = np.array(data['natural_orbital_occupancies'][l])
    return data

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

def organize_by_irrep(data):
    R_list  = data['geometries']
    mo_eig  = data['molecular_orbital_energies']
    mo_irr  = data['molecular_orbital_irreps']
    no_occ  = data['natural_orbital_occupancies']
    no_irr  = data['natural_orbital_irreps']
    nR      = len(R_list)
    irr_set = set(mo_irr[0])
    irr_num = {x:mo_irr[0].count(x) for x in irr_set}
    data_by_irr = {}
    for x in irr_set:
        data_by_irr[x] = {'geometries':R_list}
        data_by_irr[x]['molecular_orbital_energies']  = np.zeros((nR,irr_num[x]))
        data_by_irr[x]['natural_orbital_occupancies'] = np.zeros((nR,irr_num[x]))
        for j,Rj in enumerate(R_list):
            # ---
            idx_x = [i for (i,xi) in enumerate(mo_irr[j]) if xi==x]
            data_by_irr[x]['molecular_orbital_energies'][j,:] = mo_eig[j][idx_x]
            # ---
            idx_x = [i for (i,xi) in enumerate(no_irr[j]) if xi==x]
            data_by_irr[x]['natural_orbital_occupancies'][j,:] = no_occ[j][idx_x]
    # exlusion of A1 core orbital
    if('A1' in irr_set): irrep_core = 'A1'
    else: irrep_core = 'A1g'
    data_by_irr[irrep_core]['molecular_orbital_energies']  = data_by_irr[irrep_core]['molecular_orbital_energies'][:,1:]
    data_by_irr[irrep_core]['natural_orbital_occupancies'] = data_by_irr[irrep_core]['natural_orbital_occupancies'][:,1:]
    return data_by_irr

# ------------------------------------------------------------- 

data_bh = read_scf_fci_data('../../scf_fci/bh_scf_data.txt')
data_bh = organize_by_irrep(data_bh)

for irrep in data_bh.keys():
    R = data_bh[irrep]['geometries']
    E = data_bh[irrep]['molecular_orbital_energies']
    n = data_bh[irrep]['natural_orbital_occupancies']
    for i,symb in zip(range(E.shape[1]),['+','x',r'$*$']):
        ax[0,0].plot(R,E[:,i],color=style[irrep]['color'],ls=style[irrep]['s'],marker=symb)
        ax[0,1].plot(R,n[:,i],color=style[irrep]['color'],ls=style[irrep]['s'],marker=symb)
fill_panel(ax[0,0],'',[0.5,5.5],[0.5,1.5,2.5,3.5,4.5,5.5],[0.5,1.5,2.5,3.5,4.5,5.5],
                   r'$\varepsilon$ [$\mathrm{E_h}$]',[-1.0,1.5],[-1.0,-0.5,0,0.5,1.0,1.5],
                   ['-1.0','-0.5','0.0','0.5','1.0','1.5'],p=20.0,q=10.0)
fill_panel(ax[0,1],'',[0.5,5.5],[0.5,1.5,2.5,3.5,4.5,5.5],[0.5,1.5,2.5,3.5,4.5,5.5],
                   r'$n$',[0.0,2.0],[0.0,0.5,1.0,1.5,2.0],
                   ['0.0','0.5','1.0','1.5','2.0'],p=20.0,q=10.0)
ax[0,0].text(0.5,0.9,'BH',horizontalalignment='center',verticalalignment='center',transform=ax[0,0].transAxes)
ax[0,1].text(0.1,0.5,'BH',horizontalalignment='center',verticalalignment='center',transform=ax[0,1].transAxes)

# ------------------------------------------------------------- 

data_hf = read_scf_fci_data('../../scf_fci/hf_scf_data.txt')
data_hf = organize_by_irrep(data_hf)

for irrep in data_hf.keys():
    R = data_hf[irrep]['geometries']
    E = data_hf[irrep]['molecular_orbital_energies']
    n = data_hf[irrep]['natural_orbital_occupancies']
    for i,symb in zip(range(E.shape[1]),['+','x',r'$*$']):
        ax[1,0].plot(R,E[:,i],color=style[irrep]['color'],ls=style[irrep]['s'],marker=symb)
        ax[1,1].plot(R,n[:,i],color=style[irrep]['color'],ls=style[irrep]['s'],marker=symb)

fill_panel(ax[1,0],'',[0.5,4.5],[0.5,1.5,2.5,3.5,4.5],[0.5,1.5,2.5,3.5,4.5],
                   r'$\varepsilon$ [$\mathrm{E_h}$]',[-2.0,2.0],[-2.0,-1.0,0,1.0,2.0],
                   ['-2.0','-1.0','0.0','1.0','2.0'],p=20.0,q=10.0)
fill_panel(ax[1,1],'',[0.5,4.5],[0.5,1.5,2.5,3.5,4.5],[0.5,1.5,2.5,3.5,4.5],
                   r'$n$',[0.0,2.0],[0.0,0.5,1.0,1.5,2.0],
                   ['0.0','0.5','1.0','1.5','2.0'],p=20.0,q=10.0)
ax[1,0].text(0.5,0.9,'HF',horizontalalignment='center',verticalalignment='center',transform=ax[1,0].transAxes)
ax[1,1].text(0.1,0.5,'HF',horizontalalignment='center',verticalalignment='center',transform=ax[1,1].transAxes)

# ------------------------------------------------------------- 

data_beh2 = read_scf_fci_data('../../scf_fci/beh2_scf_data.txt')
data_beh2 = organize_by_irrep(data_beh2)

for irrep in data_beh2.keys():
    R = data_beh2[irrep]['geometries']
    E = data_beh2[irrep]['molecular_orbital_energies']
    n = data_beh2[irrep]['natural_orbital_occupancies']
    for i,symb in zip(range(E.shape[1]),['+','x',r'$*$']):
        ax[2,0].plot(R,E[:,i],color=style[irrep]['color'],ls=style[irrep]['s'],marker=symb)
        ax[2,1].plot(R,n[:,i],color=style[irrep]['color'],ls=style[irrep]['s'],marker=symb)
fill_panel(ax[2,0],'',[0.5,4.5],[0.5,1.5,2.5,3.5,4.5],[0.5,1.5,2.5,3.5,4.5],
                   r'$\varepsilon$ [$\mathrm{E_h}$]',[-0.5,2.0],[-0.5,0.0,0.5,1.0,1.5,2.0],
                   ['-0.5','0.0','0.5','1.0','1.5','2.0'],p=20.0,q=10.0)
fill_panel(ax[2,1],'',[0.5,4.5],[0.5,1.5,2.5,3.5,4.5],[0.5,1.5,2.5,3.5,4.5],
                   r'$n$',[0.0,2.0],[0.0,0.5,1.0,1.5,2.0],
                   ['0.0','0.5','1.0','1.5','2.0'],p=20.0,q=10.0)
ax[2,0].text(0.5,0.9,r'BeH$_2$',horizontalalignment='center',verticalalignment='center',transform=ax[2,0].transAxes)
ax[2,1].text(0.1,0.5,r'BeH$_2$',horizontalalignment='center',verticalalignment='center',transform=ax[2,1].transAxes)

# ------------------------------------------------------------- 

data_h2o = read_scf_fci_data('../../scf_fci/h2o_scf_data.txt')
data_h2o = organize_by_irrep(data_h2o)

for irrep in data_h2o.keys():
    R = data_h2o[irrep]['geometries']
    E = data_h2o[irrep]['molecular_orbital_energies']
    n = data_h2o[irrep]['natural_orbital_occupancies']
    for i,symb in zip(range(E.shape[1]),['+','x',r'$*$']):
        ax[3,0].plot(R,E[:,i],color=style[irrep]['color'],ls=style[irrep]['s'],marker=symb)
        ax[3,1].plot(R,n[:,i],color=style[irrep]['color'],ls=style[irrep]['s'],marker=symb)
fill_panel(ax[3,0],r'$R$ $[\mathrm{\AA}]$',[0.5,3.5],[0.50,1.25,2.0,2.75,3.5],['0.50','1.25','2.00','2.75','3.50'],
                   r'$\varepsilon$ [$\mathrm{E_h}$]',[-1.5,1.0],[-1.5,-1.0,-0.5,0,0.5,1.0],
                   ['-1.5','-1.0','-0.5','0.0','0.5','1.0'],p=20.0,q=10.0)
fill_panel(ax[3,1],r'$R$ $[\mathrm{\AA}]$',[0.5,3.5],[0.50,1.25,2.0,2.75,3.5],['0.50','1.25','2.00','2.75','3.50'],
                   r'$n$',[0.0,2.0],[0.0,0.5,1.0,1.5,2.0],
                   ['0.0','0.5','1.0','1.5','2.0'],p=20.0,q=10.0)
ax[3,0].text(0.5,0.9,r'H$_2$O',horizontalalignment='center',verticalalignment='center',transform=ax[3,0].transAxes)
ax[3,1].text(0.1,0.5,r'H$_2$O',horizontalalignment='center',verticalalignment='center',transform=ax[3,1].transAxes)

# -------------------------------------------------------------

for s in style:
    ax[0,0].plot([np.nan],[np.nan],color=style[s]['color'],label=style[s]['label'],ls=style[s]['s'])
h,l = ax[0,0].get_legend_handles_labels()

x0L,y0L,dxL,dyL = 0.00,1.05,2.10,0.20
ax[0,0].legend(h,l,fancybox=True,shadow=True,ncol=10,loc=3,
               bbox_to_anchor=(x0L,y0L,dxL,dyL),handlelength=1.25,handletextpad=0.5,columnspacing=1.5)

fname = 'scf_fci.eps'
fig.savefig(fname,format='eps')
os.system('ps2epsi '+fname)
os.system('mv '+fname+'i '+fname)

