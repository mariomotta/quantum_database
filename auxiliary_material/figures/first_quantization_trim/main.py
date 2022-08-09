import numpy as np
import os
import matplotlib.pyplot as plt
import itertools

from matplotlib import rc
rc('font',**{'family':'serif','serif':['cmu serif'],'size':14})
rc('text', usetex=True)

L        = 3.5
fig,ax   = plt.subplots(3,4,figsize=(4*L,0.5*3*L))
fig.subplots_adjust(hspace=0.0,wspace=0.45)

c_list = {'yellow'      : '#DAA520',
          'orange'      : '#FF8856',
          'red'         : '#DC343B',
          'magenta'     : '#B284BE',
          'purple'      : '#888FC7',
          'blue'        : '#0079BF',
          'cyan'        : '#5BA4CF',
          'light-green' : '#B9D146',
          'dark-green'  : '#32826E'}

style = {'1'  :{'color':c_list['red'],        's':'-', 'w':2,'marker':'o','ms':16,'label':r'$n_r=1$'},
         '2'  :{'color':c_list['blue'],       's':'--','w':2,'marker':'P','ms':16,'label':r'$n_r=2$'},
         '3'  :{'color':c_list['cyan'],       's':'--','w':2,'marker':'+','ms':16,'label':r'$n_r=3$'},
         '4'  :{'color':c_list['light-green'],'s':'-', 'w':2,'marker':'x','ms':13,'label':r'$n_r=4$'},
         '5'  :{'color':c_list['dark-green'], 's':'--','w':2,'marker':'3','ms':13,'label':r'$n_r=5$'},
         'fci':{'color':c_list['orange'],     's':'-', 'w':2,'marker':'*','ms':11,'label':r'FCI'},
         'x'  :{'color':c_list['purple'],     's':'-.','w':2,'marker':'^','ms':11,'label':'x'},
         'y'  :{'color':c_list['orange'],     's':'-.','w':2,'marker':'s','ms':16,'label':'y'},
         'z'  :{'color':c_list['yellow'],     's':':', 'w':2,'marker':'D','ms':16,'label':'z'}}

def fill_panel(pan,xlabel,xlim,xticks,xticklabels,ylabel,ylim,yticks,yticklabels,p=20.0,q=7.0):
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

def get_geometries(mol):
    if(mol=='bh'):   return [0.9,1.1,1.3,1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.1,3.3,3.5,3.7,3.9,4.1,4.3,4.5,4.7,4.9,5.1,5.3]
    if(mol=='hf'):   return [0.5,0.7,0.9,1.1,1.3,1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.1,3.3,3.5,3.7,3.9,4.1,4.3,4.5]
    if(mol=='beh2'): return [0.7,0.9,1.1,1.3,1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.1,3.3,3.5,3.7,3.9,4.1,4.3,4.5]
    if(mol=='h2o'):  return [0.7,0.9,1.1,1.3,1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.1,3.3]

# ------------------------------------------------------------- 

for j,mol in enumerate(['bh','hf','beh2','h2o']):
    print(mol)
    data_mol = {'geometries':get_geometries(mol)}
    data_mol['E_fci'] = np.loadtxt('../../../first_quantization/trim/operators/%s_info.txt'%mol)[:,4]
    data_mol['E_fci_trim'] = np.loadtxt('../../../first_quantization/trim/operators/%s_info.txt'%mol)[:,5]
    for depth in [3,4,5]:
        data_mol['E_cascade_%d'%depth] = np.zeros(len(data_mol['geometries']))
        for jj,R in enumerate(data_mol['geometries']):
            data_mol['E_cascade_%d'%depth][jj] = np.loadtxt('../../../first_quantization/trim/circuits/cascade_%d/%s_%s_results.txt'%(depth,mol,str(R)))[0,0]

    d_list = [3,4,5]    
    for d in d_list:
        R = np.array(data_mol['geometries'])
        E = np.array(data_mol['E_cascade_%d'%d])
        X = np.array(data_mol['E_fci_trim'])
        Y = np.array(data_mol['E_fci'])
        for Rj,Dj in zip(R,E-X):
            print(mol,d,Rj,Dj)
        ax[0,j].plot(R,E  ,color=style[str(d)]['color'],ls=style[str(d)]['s'],marker=style[str(d)]['marker'],label=style[str(d)]['label'])
        ax[1,j].plot(R,E-X,color=style[str(d)]['color'],ls=style[str(d)]['s'],marker=style[str(d)]['marker'],label=style[str(d)]['label'])
        d = 'fci'
        ax[2,j].plot(R,X-Y,color=style[str(d)]['color'],ls=style[str(d)]['s'],marker=style[str(d)]['marker']) #,label=style[str(d)]['label'])

# --------------------------------------------------

fill_panel(ax[0,0],'',[0.5,5.5],[0.5,1.5,2.5,3.5,4.5,5.5],['','','','','',''],
                   r'$E$ [$\mathrm{E_h}$]',[-25.05,-24.85],[-25.05,-25.00,-24.95,-24.90,-24.85],['-25.05','-25.00','-24.95','-24.90','-24.85'])
fill_panel(ax[1,0],'',[0.5,5.5],[0.5,1.5,2.5,3.5,4.5,5.5],['','','','','',''],
                   r'$E$-$E_{\mathrm{FCI}}^{\mathrm{trim}}$ [$m\mathrm{E_h}$]',[0.0000,0.0001],[0,2.5e-5,5e-5,7.5e-5,1e-4],['0.000','0.025','0.050','0.075','0.100'])
fill_panel(ax[2,0],r'$R$ [\AA]',[0.5,5.5],[0.5,1.5,2.5,3.5,4.5,5.5],[0.5,1.5,2.5,3.5,4.5,5.5],
                   r'$E_{\mathrm{FCI}}^{\mathrm{trim}}$-$E_{\mathrm{FCI}}$ [$m\mathrm{E_h}$]',[0.0000,0.0004],[0,1e-4,2e-4,3e-4,4e-4],['0.000','0.100','0.200','0.300','0.400'])
ax[0,0].text(0.85,0.9,'BH',horizontalalignment='center',verticalalignment='center',transform=ax[0,0].transAxes)

fill_panel(ax[0,1],'',[0.5,4.5],[0.5,1.5,2.5,3.5,4.5],['','','','',''],
                   r'$E$ [$\mathrm{E_h}$]',[-99.7,-98.5],[-99.7,-99.4,-99.1,-98.8,-98.5],['-99.70','-99.40','-99.10','-98.80','-98.50'])
fill_panel(ax[1,1],'',[0.5,4.5],[0.5,1.5,2.5,3.5,4.5],['','','','',''],
                   r'$E$-$E_{\mathrm{FCI}}^{\mathrm{trim}}$ [$m\mathrm{E_h}$]',[0.0000,0.0001],[0,2.5e-5,5e-5,7.5e-5,1e-4],['0.000','0.025','0.050','0.075','0.100'])
fill_panel(ax[2,1],r'$R$ [\AA]',[0.5,4.5],[0.5,1.5,2.5,3.5,4.5],[0.5,1.5,2.5,3.5,4.5],
                   r'$E_{\mathrm{FCI}}^{\mathrm{trim}}$-$E_{\mathrm{FCI}}$ [$m\mathrm{E_h}$]',[0.0000,0.0004],[0,1e-4,2e-4,3e-4,4e-4],['0.000','0.100','0.200','0.300','0.400'])
ax[0,1].text(0.85,0.9,'HF',horizontalalignment='center',verticalalignment='center',transform=ax[0,1].transAxes)

fill_panel(ax[0,2],'',[0.5,4.5],[0.5,1.5,2.5,3.5,4.5],['','','','',''],
                   r'$E$ [$\mathrm{E_h}$]',[-15.8,-15.0],[-15.8,-15.6,-15.4,-15.2,-15.0],['-15.80','-15.60','-15.40','-15.20','-15.00'])
fill_panel(ax[1,2],'',[0.5,4.5],[0.5,1.5,2.5,3.5,4.5],['','','','',''],
                   r'$E$-$E_{\mathrm{FCI}}^{\mathrm{trim}}$ [$m\mathrm{E_h}$]',[0.0000,0.0001],[0,2.5e-5,5e-5,7.5e-5,1e-4],['0.000','0.025','0.050','0.075','0.100'])
fill_panel(ax[2,2],r'$R$ [\AA]',[0.5,4.5],[0.5,1.5,2.5,3.5,4.5],[0.5,1.5,2.5,3.5,4.5],
                   r'$E_{\mathrm{FCI}}^{\mathrm{trim}}$-$E_{\mathrm{FCI}}$ [$m\mathrm{E_h}$]',[0.0000,0.004],[0,1e-3,2e-3,3e-3,4e-3],['0.000','1.000','2.000','3.000','4.000'])
ax[0,2].text(0.85,0.9,'BeH$_2$',horizontalalignment='center',verticalalignment='center',transform=ax[0,2].transAxes)

fill_panel(ax[0,3],'',[0.5,3.5],[0.5,1.5,2.5,3.5],['','','',''],
                   r'$E$ [$\mathrm{E_h}$]',[-75.75,-75.35],[-75.75,-75.65,-75.55,-75.45,-75.35],['-75.75','-75.65','-75.55','-75.45','-75.35'])
fill_panel(ax[1,3],'',[0.5,3.5],[0.5,1.5,2.5,3.5],['','','',''],
                   r'$E$-$E_{\mathrm{FCI}}^{\mathrm{trim}}$ [$m\mathrm{E_h}$]',[0.0000,0.0001],[0,2.5e-5,5e-5,7.5e-5,1e-4],['0.000','0.025','0.050','0.075','0.100'])
fill_panel(ax[2,3],r'$R$ [\AA]',[0.5,3.5],[0.5,1.5,2.5,3.5],[0.5,1.5,2.5,3.5],
                   r'$E_{\mathrm{FCI}}^{\mathrm{trim}}$-$E_{\mathrm{FCI}}$ [$m\mathrm{E_h}$]',[0.0000,0.0004],[0,1e-4,2e-4,3e-4,4e-4],['0.000','0.100','0.200','0.300','0.400'])
ax[0,3].text(0.85,0.9,'H$_2$O',horizontalalignment='center',verticalalignment='center',transform=ax[0,3].transAxes)

# --------------------------------------------------

h,l = ax[1,0].get_legend_handles_labels()

x0L,y0L,dxL,dyL = 0.30,0.15,0.3,0.5
ax[1,0].legend(h,l,fancybox=True,shadow=True,ncol=1,loc=3,
               bbox_to_anchor=(x0L,y0L,dxL,dyL),handlelength=1.25,handletextpad=0.5,columnspacing=1.5)

fname = 'first_quantization_trim.eps'
fig.savefig(fname,format='eps')
os.system('ps2epsi '+fname)
os.system('mv '+fname+'i '+fname)

