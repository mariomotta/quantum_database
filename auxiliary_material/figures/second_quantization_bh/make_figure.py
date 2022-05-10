import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib import rc

rc('font',**{'family':'serif','serif':['cmu serif'],'size':12})
rc('text', usetex=True)

def read_data(fname):
    return np.loadtxt(fname)

def fill_panel(pan,xlabel,xlim,xticks,xticklabels,ylabel,ylim,yticks,yticklabels,p=10.0,q=10.0):
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

def custom_plot(ax,style,d,x,y):
    ax.plot(x,y,color=style[str(d)]['color'],ls=style[str(d)]['s'],marker=style[str(d)]['marker'],label=style[str(d)]['label'],ms=1) #style[str(d)]['ms'])

c_list = {'purple'       : '#B284BE',
          'jacaranda'    : '#888FC7',
          'light_yellow' : '#EEEA62',
          'gold'         : '#FFD300',
          'earwax'       : '#DAA520',
          'brown'        : '#7B3F00',
          'light_blue'   : '#BCD9EA',
          'cerulean'     : '#5BA4CF',
          'cobalt'       : '#0079BF',
          'dark_blue'    : '#055A8C',
          'light_green'  : '#ACDF87',
          'yellow_green' : '#B9D146',
          'mid_green'    : '#68BB59',
          'dark_green'   : '#1E5631',
          'orange'       : '#FF8856',
          'red'          : '#DC343B',
          'light-gray'   : '#adadad',
          'palatinate'   : '#72246C',
          'black'        : 'black'}

style = {'scf':{'color':c_list['light-gray'],  's':':', 'w':2,'marker':'o','ms':3,'label':'SCF'},
         '1'  :{'color':c_list['red'],         's':'-', 'w':2,'marker':'+','ms':5,'label':r'$n_r=1$'},
         '2'  :{'color':c_list['cobalt'],      's':'--','w':2,'marker':'x','ms':7,'label':r'$n_r=2$'},
         '3'  :{'color':c_list['orange'],      's':'-', 'w':2,'marker':'+','ms':5,'label':r'$n_r=3$'},
         '4'  :{'color':c_list['mid_green'],   's':'--','w':2,'marker':'x','ms':7,'label':r'$n_r=4$'},
         '5'  :{'color':c_list['earwax'],      's':'-', 'w':2,'marker':'+','ms':5,'label':r'$n_r=5$'},
         '6'  :{'color':c_list['yellow_green'],'s':'--','w':2,'marker':'x','ms':7,'label':r'$n_r=6$'},
         '7'  :{'color':c_list['gold'],        's':'-', 'w':2,'marker':'+','ms':5,'label':r'$n_r=7$'},
         '8'  :{'color':c_list['jacaranda'],   's':'--','w':2,'marker':'x','ms':7,'label':r'$n_r=8$'},
         'fci':{'color':c_list['black'],       's':'-.','w':2,'marker':'o','ms':1,'label':'FCI'}}

# -----------------------------------------------------------------------------------------------------------

L        = 2.7
fig,ax   = plt.subplots(5,4,figsize=(4*L,0.5*5*L))
fig.subplots_adjust(hspace=0.0,wspace=0.0)

scf_data = np.loadtxt('/Users/mario/Documents/GitHub/VATech/quantum_database/auxiliary_material/scf_fci/bh_scf_energy.txt')

path = '/Users/mario/Documents/GitHub/VATech/quantum_database/second_quantization/circuits/hardware_efficient/'

dist_list = [0.7,0.9,1.1,1.3,1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.1,3.3,3.5,3.7,3.9,4.1,4.3,4.5,4.7,4.9,5.1,5.3]
depth_list = [1,2,3,4,5,6,7,8]
ndist,ndepth = len(dist_list),len(depth_list)

data = np.zeros((ndist,ndepth,4,2))
for i,depth in enumerate(depth_list):
    for j,dist in enumerate(dist_list):
        data[j,i,:,:] = read_data(path+'ry_linear_%d/bh_%s_results.txt' % (depth,dist))

for i,depth in enumerate(depth_list):
    custom_plot(ax[0,0],style,depth,dist_list,data[:,i,0,0])                   # --- energy
    custom_plot(ax[1,0],style,depth,dist_list,data[:,i,0,0]-data[:,i,0,1])     # --- dev from fci
    custom_plot(ax[2,0],style,depth,dist_list,data[:,i,1,0]-data[:,i,1,1])     # --- n
    custom_plot(ax[3,0],style,depth,dist_list,data[:,i,2,0]-data[:,i,2,1])     # --- sz
    custom_plot(ax[4,0],style,depth,dist_list,data[:,i,3,0]-data[:,i,3,1])     # --- s2
    if(depth==8): custom_plot(ax[0,0],style,'scf',scf_data[:,0],scf_data[:,1]) # --- scf
    if(depth==8): custom_plot(ax[0,0],style,'fci',dist_list,data[:,i,0,1])     # --- fci
    if(depth==8): custom_plot(ax[1,0],style,'scf',scf_data[:,0],scf_data[:,1]-data[:,i,0,1])

# -----------------------------------------------------------------------------------------------------------

depth_list = [1,2,3,4]
ndist,ndepth = len(dist_list),len(depth_list)

data = np.zeros((ndist,ndepth,4,2))
for i,depth in enumerate(depth_list):
    for j,dist in enumerate(dist_list):
        data[j,i,:,:] = read_data(path+'ry_full_%d/bh_%s_results.txt' % (depth,dist))

for i,depth in enumerate(depth_list):
    custom_plot(ax[0,1],style,depth,dist_list,data[:,i,0,0])                   # --- energy
    custom_plot(ax[1,1],style,depth,dist_list,data[:,i,0,0]-data[:,i,0,1])     # --- dev from fci
    custom_plot(ax[2,1],style,depth,dist_list,data[:,i,1,0]-data[:,i,1,1])     # --- n
    custom_plot(ax[3,1],style,depth,dist_list,data[:,i,2,0]-data[:,i,2,1])     # --- sz
    custom_plot(ax[4,1],style,depth,dist_list,data[:,i,3,0]-data[:,i,3,1])     # --- s2
    if(depth==4): custom_plot(ax[0,1],style,'scf',scf_data[:,0],scf_data[:,1]) # --- scf
    if(depth==4): custom_plot(ax[0,1],style,'fci',dist_list,data[:,i,0,1])     # --- fci
    if(depth==4): custom_plot(ax[1,1],style,'scf',scf_data[:,0],scf_data[:,1]-data[:,i,0,1])

# -----------------------------------------------------------------------------------------------------------

depth_list = [1,2,3,4]
ndist,ndepth = len(dist_list),len(depth_list)

data = np.zeros((ndist,ndepth,4,2))
for i,depth in enumerate(depth_list):
    for j,dist in enumerate(dist_list):
        data[j,i,:,:] = read_data(path+'cascade_%d/bh_%s_results.txt' % (depth,dist))

for i,depth in enumerate(depth_list):
    custom_plot(ax[0,2],style,depth,dist_list,data[:,i,0,0])                   # --- energy
    custom_plot(ax[1,2],style,depth,dist_list,data[:,i,0,0]-data[:,i,0,1])     # --- dev from fci
    custom_plot(ax[2,2],style,depth,dist_list,data[:,i,1,0]-data[:,i,1,1])     # --- n
    custom_plot(ax[3,2],style,depth,dist_list,data[:,i,2,0]-data[:,i,2,1])     # --- sz
    custom_plot(ax[4,2],style,depth,dist_list,data[:,i,3,0]-data[:,i,3,1])     # --- s2
    if(depth==4): custom_plot(ax[0,2],style,'scf',scf_data[:,0],scf_data[:,1]) # --- scf
    if(depth==4): custom_plot(ax[0,2],style,'fci',dist_list,data[:,i,0,1])     # --- fci
    if(depth==4): custom_plot(ax[1,2],style,'scf',scf_data[:,0],scf_data[:,1]-data[:,i,0,1])

# -----------------------------------------------------------------------------------------------------------

path = '/Users/mario/Documents/GitHub/VATech/quantum_database/second_quantization/circuits/qUCCSD/'

depth_list = [1]
ndist,ndepth = len(dist_list),len(depth_list)

data = np.zeros((ndist,ndepth,4,2))

for i,depth in enumerate(depth_list):
    for j,dist in enumerate(dist_list):
        data[j,i,:,:] = read_data(path+'unrestricted/singles_doubles/trotter/bh_reps_%d_%s_results.txt' % (depth,dist))

for i,depth in enumerate(depth_list):
    custom_plot(ax[0,3],style,depth,dist_list,data[:,i,0,0])                   # --- energy
    custom_plot(ax[1,3],style,depth,dist_list,data[:,i,0,0]-data[:,i,0,1])     # --- dev from fci
    custom_plot(ax[2,3],style,depth,dist_list,data[:,i,1,0]-data[:,i,1,1])     # --- n
    custom_plot(ax[3,3],style,depth,dist_list,data[:,i,2,0]-data[:,i,2,1])     # --- sz
    custom_plot(ax[4,3],style,depth,dist_list,data[:,i,3,0]-data[:,i,3,1])     # --- s2
    if(depth==1): custom_plot(ax[0,3],style,'scf',scf_data[:,0],scf_data[:,1]) # --- scf
    if(depth==1): custom_plot(ax[0,3],style,'fci',dist_list,data[:,i,0,1])     # --- fci
    if(depth==1): custom_plot(ax[1,3],style,'scf',scf_data[:,0],scf_data[:,1]-data[:,i,0,1])

for c in [0,1,2,3]:
    if(c>0): lab = ['']*5
    else:    lab = ['-25.1','-25.0','-24.9','-24.8','-24.7']
    if(c>0): yl = ''
    else:    yl = r'$E$ [$\mathrm{E_h}$]'
    fill_panel(ax[0,c],'',[0.5,5.5],[0.5,1.5,2.5,3.5,4.5,5.5],['','','','','',''],yl,[-25.1,-24.7],[-25.1,-25.0,-24.9,-24.8,-24.7],lab)
    # ---------------- 
    if(c>0): lab = ['']*5
    else:    lab = ['0.00','0.02','0.04','0.06','0.08']
    if(c>0): yl = ''
    else:    yl = r'$E$-$E_{\mathrm{FCI}}$ [$\mathrm{E_h}$]'
    fill_panel(ax[1,c],'',[0.5,5.5],[0.5,1.5,2.5,3.5,4.5,5.5],['','','','','',''],yl,[0,0.08],[0,0.02,0.04,0.06,0.08],lab)
    # ---------------- 
    if(c>0): lab = ['']*5
    else:    lab = ['-0.002','-0.001','0.000','0.001','0.002']
    if(c>0): yl = ''
    else:    yl = r'$N$-$N_{\mathrm{FCI}}$'
    fill_panel(ax[2,c],'',[0.5,5.5],[0.5,1.5,2.5,3.5,4.5,5.5],['','','','','',''],yl,[-0.002,0.002],[-0.002,-0.001,0,0.001,0.002],lab)
    # ---------------- 
    if(c>0): lab = ['']*5
    else:    lab = ['0.00','0.50','1.00','1.50','2.00']
    if(c>0): yl = ''
    else:    yl = r'$S^2$-$S^2_{\mathrm{FCI}}$ $[\hbar]$'
    fill_panel(ax[3,c],'',[0.5,5.5],[0.5,1.5,2.5,3.5,4.5,5.5],['','','','','',''],yl,[0,2],[0,0.5,1.0,1.5,2.0],lab)
    # ---------------- 
    if(c>0): lab = ['']*5
    else:    lab = ['-0.002','-0.001','0.000','0.001','0.002']
    if(c>0): yl = ''
    else:    yl = r'$S_z$-$S_{z,\mathrm{FCI}}$ $[\hbar]$'
    fill_panel(ax[4,c],r'$R$ $[\mathrm{\AA}]$',[0.5,5.5],[0.5,1.5,2.5,3.5,4.5,5.5],
               ['0.5','1.5','2.5','3.5','4.5','5.5'],yl,[-0.002,0.002],[-0.002,-0.001,0,0.001,0.002],lab)
    # ---------------- 
    ax[0,0].text(0.75,0.75,'$R_y$ linear',horizontalalignment='center',verticalalignment='center',transform=ax[0,0].transAxes,fontsize=14)
    ax[0,1].text(0.75,0.75,  '$R_y$ full',horizontalalignment='center',verticalalignment='center',transform=ax[0,1].transAxes,fontsize=14)
    ax[0,2].text(0.75,0.75,     'cascade',horizontalalignment='center',verticalalignment='center',transform=ax[0,2].transAxes,fontsize=14)
    ax[0,3].text(0.75,0.75,     'q-UCCSD',horizontalalignment='center',verticalalignment='center',transform=ax[0,3].transAxes,fontsize=14)

h,l = ax[0,0].get_legend_handles_labels()
x0L,y0L,dxL,dyL = -0.02,1.05,4*L,0.5
ax[0,0].legend(h,l,fancybox=True,shadow=True,ncol=10,loc=3,
               bbox_to_anchor=(x0L,y0L,dxL,dyL),handlelength=1.4,handletextpad=0.5,columnspacing=0.5)

fname = 'second_quantization_bh.eps'
fig.savefig(fname,format='eps')
os.system('ps2epsi '+fname)
os.system('mv '+fname+'i '+fname)

