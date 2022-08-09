import numpy as np
import os
import matplotlib.pyplot as plt
import itertools

from matplotlib import rc
rc('font',**{'family':'serif','serif':['cmu serif'],'size':14})
rc('text', usetex=True)

c_list = {'yellow'      : '#DAA520',
          'orange'      : '#FF8856',
          'red'         : '#DC343B',
          'magenta'     : '#B284BE',
          'purple'      : '#888FC7',
          'blue'        : '#0079BF',
          'cyan'        : '#5BA4CF',
          'light-green' : '#B9D146',
          'dark-green'  : '#32826E'}

style = {'trotter_ds':{'color':c_list['red'],        's':'-', 'w':2,'marker':'o','ms':16,'label':'Trotter, d+s'},
         'suzuki_ds' :{'color':c_list['blue'],       's':'-', 'w':2,'marker':'P','ms':16,'label':'Suzuki, d+s'},
         'suzuki_sd' :{'color':c_list['light-green'],'s':'-.','w':2,'marker':'<','ms':13,'label':'Suzuki, s+d'},
         'trotter_sd':{'color':c_list['orange'],     's':'-.','w':2,'marker':'s','ms':16,'label':'Trotter, s+d'}}

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

R_list  = [0.7,0.9,1.1,1.3,1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.1,3.3]
flavor  = ['restricted_closed_shell','unrestricted']
approx  = ['trotter','suzuki']
reverse = ['doubles_singles','singles_doubles']
title   = ['restricted','unrestricted']
dlist   = [1,2]

for depth in dlist:

    L        = 4.5
    fig,ax   = plt.subplots(2,2,figsize=(2*L,0.5*2*L))
    fig.subplots_adjust(hspace=0.0,wspace=0.3)

    res = {}
    for f in flavor:
        for a in approx:
            for x in reverse:
                res['E_vqe_%s_%s_ds_%s'%(f,a,x)] = np.zeros(len(R_list))
                res['S_vqe_%s_%s_ds_%s'%(f,a,x)] = np.zeros(len(R_list))
                res['E_fci_%s_%s_ds_%s'%(f,a,x)] = np.zeros(len(R_list))
                for j,R in enumerate(R_list):
                    DAT = np.loadtxt('../../../second_quantization/circuits/qUCCSD/%s/%s/%s/h2o_reps_%d_%s_results.txt' % (f,x,a,depth,str(R)))
                    res['E_vqe_%s_%s_ds_%s'%(f,a,x)][j] = DAT[0][0]
                    res['S_vqe_%s_%s_ds_%s'%(f,a,x)][j] = DAT[2][0]
                    res['E_fci_%s_%s_ds_%s'%(f,a,x)][j] = DAT[0][1]
    
    for j,f in enumerate(flavor):
        for a in approx:
            for x in reverse:
                E_vqe = np.array(res['E_vqe_%s_%s_ds_%s'%(f,a,x)])
                S_vqe = np.array(res['S_vqe_%s_%s_ds_%s'%(f,a,x)])
                E_fci = np.array(res['E_fci_%s_%s_ds_%s'%(f,a,x)])
                if(x=='doubles_singles'): s = a+'_ds'
                else:                     s = a+'_sd'
                ax[0,j].plot(R_list,E_vqe-E_fci,color=style[s]['color'],label=style[s]['label'],ls=style[s]['s'],marker=style[s]['marker'],ms=5,mec='black',mew=0.5)
                ax[1,j].plot(R_list,S_vqe,color=style[s]['color'],label=style[s]['label'],ls=style[s]['s'],marker=style[s]['marker'],ms=5,mec='black',mew=0.5)
        ax[0,j].text(0.1,0.9,title[j],horizontalalignment='left',verticalalignment='center',transform=ax[0,j].transAxes)
    
    fill_panel(ax[0,0],'',[0.5,3.5],[0.5,1.0,1.5,2.0,2.5,3.0,3.5],['','','','','','',''],
                       r'$E-E_{\mathrm{FCI}}$ [$\mathrm{E_h}$]',[0,0.09],[0.00,0.03,0.06,0.09],
                       ['0.00','0.03','0.06','0.09'],p=20.0,q=10.0)
    fill_panel(ax[1,0],r'$R$ [\AA]',[0.5,3.5],[0.5,1.0,1.5,2.0,2.5,3.0,3.5],[0.5,1.0,1.5,2.0,2.5,3.0,3.5],
                       r'$S^2$ [$\hbar$]',[0,0.008],[0.00,0.002,0.004,0.006,0.008],
                       ['0.000','0.002','0.004','0.006','0.008'],p=20.0,q=10.0)
    
    fill_panel(ax[0,1],'',[0.5,3.5],[0.5,1.0,1.5,2.0,2.5,3.0,3.5],['','','','','','',''],
                       r'$E-E_{\mathrm{FCI}}$ [$\mathrm{E_h}$]',[0,0.009],[0.000,0.003,0.006,0.009],
                       ['0.000','0.003','0.006','0.009'],p=20.0,q=10.0)
    fill_panel(ax[1,1],r'$R$ [\AA]',[0.5,3.5],[0.5,1.0,1.5,2.0,2.5,3.0,3.5],[0.5,1.0,1.5,2.0,2.5,3.0,3.5],
                       r'$S^2$ [$\hbar$]',[0,2],[0.0,0.5,1.0,1.5,2.0],
                       ['0.0','0.5','1.0','1.5','2.0'],p=20.0,q=10.0)
    
    # --------------------------------------------------
    
    h,l = ax[0,0].get_legend_handles_labels()
    x0L,y0L,dxL,dyL = -0.02,1.05,2.20,0.20
    ax[0,0].legend(h,l,fancybox=True,shadow=True,ncol=10,loc=3,
                   bbox_to_anchor=(x0L,y0L,dxL,dyL),handlelength=1.5,handletextpad=0.5,columnspacing=2.0)
    
    fname = 'quccsd_reps_%d.eps' % depth
    fig.savefig(fname,format='eps')
    os.system('ps2epsi '+fname)
    os.system('mv '+fname+'i '+fname)

