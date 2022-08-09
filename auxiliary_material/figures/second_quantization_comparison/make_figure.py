import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib import rc
rc('font',**{'family':'serif','serif':['cmu serif'],'size':14})
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
    mk = style[str(d)]['marker']
    if(mk=='o' or mk=='s'): ax.plot(x,y,color=style[str(d)]['color'],ls=style[str(d)]['s'],marker=mk,label=style[str(d)]['label'],ms=style[str(d)]['ms'],mec='black',mew=0.5)
    else:                   ax.plot(x,y,color=style[str(d)]['color'],ls=style[str(d)]['s'],marker=mk,label=style[str(d)]['label'],ms=style[str(d)]['ms'])


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
          'light-gray'   : '#C0C0C0',
          'palatinate'   : '#72246C',
          'black'        : 'black'}

style = {'ry_linear':{'color':c_list['red'],      's':'--', 'w':2,'marker':'o','ms':4,'label':r'R$_y$, linear'},
         'qUCCSD'   :{'color':c_list['cobalt'],   's':'-',  'w':2,'marker':'+','ms':4,'label':r'q-UCCSD'      },
         'ry_full'  :{'color':c_list['orange'],   's':'-.', 'w':2,'marker':'s','ms':4,'label':r'R$_y$, full'  },
         'cascade'  :{'color':c_list['mid_green'],'s':'--', 'w':2,'marker':'x','ms':4,'label':r'cascade'      }}

# -----------------------------------------------------------------------------------------------------------

L        = 3.7
fig,ax   = plt.subplots(2,2,figsize=(2*L,0.5*2*L))
fig.subplots_adjust(hspace=0.0,wspace=0.0,left=0.1,bottom=0.15)

scf_data = np.loadtxt('/Users/mario/Documents/GitHub/VATech/quantum_database/auxiliary_material/scf_fci/bh_scf_energy.txt')
path = '/Users/mario/Documents/GitHub/VATech/quantum_database/second_quantization/circuits/hardware_efficient/'

dist_list = [0.7,0.9,1.1,1.3,1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.1,3.3,3.5,3.7,3.9,4.1,4.3,4.5]
depth_list = [8]
ndist,ndepth = len(dist_list),len(depth_list)

data = np.zeros((ndist,ndepth,4,2))
for i,depth in enumerate(depth_list):
    for j,dist in enumerate(dist_list):
        data[j,i,:,:] = read_data(path+'ry_linear_%d/bh_%s_results.txt' % (depth,dist))

for i,depth in enumerate(depth_list):
    custom_plot(ax[0,0],style,'ry_linear',dist_list,data[:,i,0,0]-data[:,i,0,1])     # --- dev from fci

# -----------------------------------------------------------------------------------------------------------

depth_list = [8]
ndist,ndepth = len(dist_list),len(depth_list)

data = np.zeros((ndist,ndepth,4,2))
for i,depth in enumerate(depth_list):
    for j,dist in enumerate(dist_list):
        data[j,i,:,:] = read_data(path+'ry_full_%d/bh_%s_results.txt' % (depth,dist))

for i,depth in enumerate(depth_list):
    custom_plot(ax[0,0],style,'ry_full',dist_list,data[:,i,0,0]-data[:,i,0,1])     # --- dev from fci

# -----------------------------------------------------------------------------------------------------------

depth_list = [4]
ndist,ndepth = len(dist_list),len(depth_list)

data = np.zeros((ndist,ndepth,4,2))
for i,depth in enumerate(depth_list):
    for j,dist in enumerate(dist_list):
        data[j,i,:,:] = read_data(path+'cascade_%d/bh_%s_results.txt' % (depth,dist))

for i,depth in enumerate(depth_list):
    custom_plot(ax[0,0],style,'cascade',dist_list,data[:,i,0,0]-data[:,i,0,1])     # --- dev from fci

# -----------------------------------------------------------------------------------------------------------

path = '/Users/mario/Documents/GitHub/VATech/quantum_database/second_quantization/circuits/qUCCSD/'

depth_list = [1]
ndist,ndepth = len(dist_list),len(depth_list)

data = np.zeros((ndist,ndepth,4,2))

for i,depth in enumerate(depth_list):
    for j,dist in enumerate(dist_list):
        data[j,i,:,:] = read_data(path+'unrestricted/singles_doubles/trotter/bh_reps_%d_%s_results.txt' % (depth,dist))

for i,depth in enumerate(depth_list):
    custom_plot(ax[0,0],style,'qUCCSD',dist_list,data[:,i,0,0]-data[:,i,0,1])     # --- dev from fci

# ===============================================================

scf_data = np.loadtxt('/Users/mario/Documents/GitHub/VATech/quantum_database/auxiliary_material/scf_fci/beh2_scf_energy.txt')

path = '/Users/mario/Documents/GitHub/VATech/quantum_database/second_quantization/circuits/hardware_efficient/'

dist_list = [0.7,0.9,1.1,1.3,1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.1,3.3,3.5,3.7,3.9,4.1,4.3,4.5]
depth_list = [10]
ndist,ndepth = len(dist_list),len(depth_list)

data = np.zeros((ndist,ndepth,4,2))
for i,depth in enumerate(depth_list):
    for j,dist in enumerate(dist_list):
        data[j,i,:,:] = read_data(path+'ry_linear_%d/beh2_%s_results.txt' % (depth,dist))

for i,depth in enumerate(depth_list):
    custom_plot(ax[0,1],style,'ry_linear',dist_list,data[:,i,0,0]-data[:,i,0,1])     # --- dev from fci

# -----------------------------------------------------------------------------------------------------------

depth_list = [10]
ndist,ndepth = len(dist_list),len(depth_list)

data = np.zeros((ndist,ndepth,4,2))
for i,depth in enumerate(depth_list):
    for j,dist in enumerate(dist_list):
        data[j,i,:,:] = read_data(path+'ry_full_%d/beh2_%s_results.txt' % (depth,dist))

for i,depth in enumerate(depth_list):
    custom_plot(ax[0,1],style,'ry_full',dist_list,data[:,i,0,0]-data[:,i,0,1])     # --- dev from fci

# -----------------------------------------------------------------------------------------------------------

depth_list = [6]
ndist,ndepth = len(dist_list),len(depth_list)

data = np.zeros((ndist,ndepth,4,2))
for i,depth in enumerate(depth_list):
    for j,dist in enumerate(dist_list):
        data[j,i,:,:] = read_data(path+'cascade_%d/beh2_%s_results.txt' % (depth,dist))

for i,depth in enumerate(depth_list):
    custom_plot(ax[0,1],style,'cascade',dist_list,data[:,i,0,0]-data[:,i,0,1])     # --- dev from fci

# -----------------------------------------------------------------------------------------------------------

path = '/Users/mario/Documents/GitHub/VATech/quantum_database/second_quantization/circuits/qUCCSD/'

depth_list = [2]
ndist,ndepth = len(dist_list),len(depth_list)

data = np.zeros((ndist,ndepth,4,2))

for i,depth in enumerate(depth_list):
    for j,dist in enumerate(dist_list):
        data[j,i,:,:] = read_data(path+'unrestricted/singles_doubles/trotter/beh2_reps_%d_%s_results.txt' % (depth,dist))

for i,depth in enumerate(depth_list):
    custom_plot(ax[0,1],style,'qUCCSD',dist_list,data[:,i,0,0]-data[:,i,0,1])     # --- dev from fci

#h,l = ax[0,0].get_legend_handles_labels()
#x0L,y0L,dxL,dyL = 0.50,0.30,0.5,0.5
#ax[0,0].legend(h,l,fancybox=True,shadow=True,ncol=1,loc=3,
#               bbox_to_anchor=(x0L,y0L,dxL,dyL),handlelength=1.5,handletextpad=1.0,columnspacing=1.0,labelspacing=0.1)

h,l = ax[0,0].get_legend_handles_labels()
x0L,y0L,dxL,dyL = -0.03,1.02,L,0.5
ax[0,0].legend(h,l,fancybox=True,shadow=True,ncol=6,loc=3,
               bbox_to_anchor=(x0L,y0L,dxL,dyL),handlelength=2.0,handletextpad=0.5,columnspacing=1.6,labelspacing=0.1)

# ===============================================================

scf_data = np.loadtxt('/Users/mario/Documents/GitHub/VATech/quantum_database/auxiliary_material/scf_fci/hf_scf_energy.txt')

path = '/Users/mario/Documents/GitHub/VATech/quantum_database/second_quantization/circuits/hardware_efficient/'

dist_list = [0.7,0.9,1.1,1.3,1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.1,3.3,3.5,3.7,3.9,4.1,4.3,4.5]
depth_list = [5]
ndist,ndepth = len(dist_list),len(depth_list)

data = np.zeros((ndist,ndepth,4,2))
for i,depth in enumerate(depth_list):
    for j,dist in enumerate(dist_list):
        data[j,i,:,:] = read_data(path+'ry_linear_%d/hf_%s_results.txt' % (depth,dist))

for i,depth in enumerate(depth_list):
    custom_plot(ax[1,0],style,'ry_linear',dist_list,data[:,i,0,0]-data[:,i,0,1])     # --- dev from fci

# -----------------------------------------------------------------------------------------------------------

depth_list = [3]
ndist,ndepth = len(dist_list),len(depth_list)

data = np.zeros((ndist,ndepth,4,2))
for i,depth in enumerate(depth_list):
    for j,dist in enumerate(dist_list):
        data[j,i,:,:] = read_data(path+'ry_full_%d/hf_%s_results.txt' % (depth,dist))

for i,depth in enumerate(depth_list):
    custom_plot(ax[1,0],style,'ry_full',dist_list,data[:,i,0,0]-data[:,i,0,1])     # --- dev from fci

# -----------------------------------------------------------------------------------------------------------

depth_list = [3]
ndist,ndepth = len(dist_list),len(depth_list)

data = np.zeros((ndist,ndepth,4,2))
for i,depth in enumerate(depth_list):
    for j,dist in enumerate(dist_list):
        data[j,i,:,:] = read_data(path+'cascade_%d/hf_%s_results.txt' % (depth,dist))

for i,depth in enumerate(depth_list):
    custom_plot(ax[1,0],style,'cascade',dist_list,data[:,i,0,0]-data[:,i,0,1])     # --- dev from fci

# -----------------------------------------------------------------------------------------------------------

path = '/Users/mario/Documents/GitHub/VATech/quantum_database/second_quantization/circuits/qUCCSD/'

depth_list = [1]
ndist,ndepth = len(dist_list),len(depth_list)

data = np.zeros((ndist,ndepth,4,2))

for i,depth in enumerate(depth_list):
    for j,dist in enumerate(dist_list):
        data[j,i,:,:] = read_data(path+'unrestricted/singles_doubles/trotter/hf_reps_%d_%s_results.txt' % (depth,dist))

for i,depth in enumerate(depth_list):
    custom_plot(ax[1,0],style,'qUCCSD',dist_list,data[:,i,0,0]-data[:,i,0,1])     # --- dev from fci

# ===============================================================================

scf_data = np.loadtxt('/Users/mario/Documents/GitHub/VATech/quantum_database/auxiliary_material/scf_fci/h2o_scf_energy.txt')

path = '/Users/mario/Documents/GitHub/VATech/quantum_database/second_quantization/circuits/hardware_efficient/'

dist_list = [0.7,0.9,1.1,1.3,1.5,1.7,1.9,2.1,2.3,2.5,2.7,2.9,3.1,3.3]
depth_list = [8]
ndist,ndepth = len(dist_list),len(depth_list)

data = np.zeros((ndist,ndepth,4,2))
for i,depth in enumerate(depth_list):
    for j,dist in enumerate(dist_list):
        data[j,i,:,:] = read_data(path+'ry_linear_%d/h2o_%s_results.txt' % (depth,dist))

for i,depth in enumerate(depth_list):
    custom_plot(ax[1,1],style,'ry_linear',dist_list,data[:,i,0,0]-data[:,i,0,1])     # --- dev from fci

# -----------------------------------------------------------------------------------------------------------

depth_list = [7]
ndist,ndepth = len(dist_list),len(depth_list)

data = np.zeros((ndist,ndepth,4,2))
for i,depth in enumerate(depth_list):
    for j,dist in enumerate(dist_list):
        data[j,i,:,:] = read_data(path+'ry_full_%d/h2o_%s_results.txt' % (depth,dist))

for i,depth in enumerate(depth_list):
    custom_plot(ax[1,1],style,'ry_full',dist_list,data[:,i,0,0]-data[:,i,0,1])     # --- dev from fci

# -----------------------------------------------------------------------------------------------------------

depth_list = [4]
ndist,ndepth = len(dist_list),len(depth_list)

data = np.zeros((ndist,ndepth,4,2))
for i,depth in enumerate(depth_list):
    for j,dist in enumerate(dist_list):
        data[j,i,:,:] = read_data(path+'cascade_%d/h2o_%s_results.txt' % (depth,dist))

for i,depth in enumerate(depth_list):
    custom_plot(ax[1,1],style,'cascade',dist_list,data[:,i,0,0]-data[:,i,0,1])     # --- dev from fci

# -----------------------------------------------------------------------------------------------------------

path = '/Users/mario/Documents/GitHub/VATech/quantum_database/second_quantization/circuits/qUCCSD/'

depth_list = [2]
ndist,ndepth = len(dist_list),len(depth_list)

data = np.zeros((ndist,ndepth,4,2))

for i,depth in enumerate(depth_list):
    for j,dist in enumerate(dist_list):
        data[j,i,:,:] = read_data(path+'unrestricted/singles_doubles/trotter/h2o_reps_%d_%s_results.txt' % (depth,dist))

for i,depth in enumerate(depth_list):
    custom_plot(ax[1,1],style,'qUCCSD',dist_list,data[:,i,0,0]-data[:,i,0,1])     # --- dev from fci

for r in [0,1]:
    for c in [0,1]:
        ax[r,c].axhline(0.0016,ls=':',c='black')

yl = r'$E$-$E_{\mathrm{FCI}}$ [$m\mathrm{E_h}$]'
la = [' 0.0','2.0','4.0','6.0']
fill_panel( ax[0,0],                     '',[0.5,4.5],[0.5,1.5,2.5,3.5,4.5],['','','','',''],               yl,[0,0.006],[0,0.002,0.004,0.006],la,p=20.0)
fill_panel( ax[1,0],r'$R$ [$\mathrm{\AA}$]',[0.5,4.5],[0.5,1.5,2.5,3.5,4.5],['0.5','1.5','2.5','3.5','4.5'],yl,[0,0.006],[0,0.002,0.004,0.006],la,p=20.0)
DUMMY1 = ax[0,1].twinx()                                                                                              
fill_panel(  DUMMY1,r'$R$ [$\mathrm{\AA}$]',[0.5,4.5],[0.5,1.5,2.5,3.5,4.5],['','','','',''],               yl,[0,0.006],[0,0.002,0.004,0.006],la,p=20.0)
fill_panel( ax[0,1],                     '',[0.5,4.5],[0.5,1.5,2.5,3.5,4.5],['','','','',''],               '',[0,0.006],[],[],p=20.0)
DUMMY2 = ax[1,1].twinx()                                                                                              
fill_panel(  DUMMY2,r'$R$ [$\mathrm{\AA}$]',[0.5,4.5],[0.5,1.5,2.5,3.5,4.5],['0.5','1.5','2.5','3.5','4.5'],yl,[0,0.006],[0,0.002,0.004,0.006],la,p=20.0)
fill_panel( ax[1,1],r'$R$ [$\mathrm{\AA}$]',[0.5,4.5],[0.5,1.5,2.5,3.5,4.5],['0.5','1.5','2.5','3.5','4.5'],'',[0,0.006],[],[],p=20.0)

ax[0,0].text(0.15,0.85,     'BH',horizontalalignment='center',verticalalignment='center',transform=ax[0,0].transAxes)
ax[0,1].text(0.15,0.85,'BeH$_2$',horizontalalignment='center',verticalalignment='center',transform=ax[0,1].transAxes)
ax[1,0].text(0.15,0.85,     'HF',horizontalalignment='center',verticalalignment='center',transform=ax[1,0].transAxes)
ax[1,1].text(0.15,0.85, 'H$_2$O',horizontalalignment='center',verticalalignment='center',transform=ax[1,1].transAxes)

fname = 'second_quantization_comparison.eps'
fig.savefig(fname,format='eps')
os.system('ps2epsi '+fname)
os.system('mv '+fname+'i '+fname)

