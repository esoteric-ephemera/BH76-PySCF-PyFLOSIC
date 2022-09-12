import yaml
import numpy as np
from os import path
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

from analysis import BH76_analysis,rdir

x_l = np.array([0, 10, 25, 50, 75, 100])
Nx = len(x_l)

scan_dir = './scan_hybs/'
hf_dir = './results_aug-cc-pvqz/'

dfad = {
    'LSDA': 'LSDA', 'PBE': 'PBE', 'SCAN': 'SCAN', 'r2SCAN': 'r$^2$SCAN',
    'LCwPBE': 'LC-$\\omega$PBE', 'SX': 'SX-$x$'
}

lbld = {'LSDA': [0.50,12.25],
    'PBE': [0.10,9.1],
    'LCwPBE': [0.50,.8],
    'SX': [0.10,3.6], 'r2SCAN': [0.35,4.6], 'SCAN': [0.80,5.9,0.75,4.96]
}

cls = ['darkblue','darkorange','tab:green','darkred','k','gray']
lsls = ['-','--','-.',':']

width = 8
aspct = 10./16.#0.55
fsz = 16
fig,ax = plt.subplots(figsize=(width,width*aspct))

ndfa = 1
for dfa in dfad:
    #if dfa == 'SX':
    #    continue
    ndfa += 1


scf_mads = {'SX': 0.0}
for dfa in dfad:
    if dfa == 'SX':
        continue
    tdir = hf_dir + '/' + dfa + '/'+ dfa + '_BH76/'
    tfl = tdir + 'BH76_total.yaml'
    if not path.isfile(tfl):
        BH76_analysis(cdir=tdir)

    td = yaml.load(open(tfl,'r'),Loader=yaml.Loader)
    scf_mads[dfa] = td['Stats']['MAD']

for idfa,dfa in enumerate(dfad):

    tlist = [0.0 for ix in range(Nx)]
    for ix,x in enumerate(x_l):
        if dfa == 'SX':
            if x == 0:
                tdir = scan_dir + '/S00X/SCAN_BH76/'
            else:
                tdir = scan_dir + '/S{:}X/'.format(int(x)) + 'SCAN_{:}_EXX_BH76/'.format(int(x))
        else:
            if x == 0:
                if dfa == 'SCAN':
                    tdir = scan_dir + '/S00X/SCAN_BH76/'
                else:
                    tdir = scan_dir + '/S00X/'+ dfa + '@SCAN_BH76/'
            else:
                tdir = scan_dir + '/S{:}X/'.format(int(x)) + dfa + '@S{:}X_BH76/'.format(int(x))

        tfl = tdir + 'BH76_total.yaml'
        if not path.isfile(tfl):
            BH76_analysis(cdir=tdir)

        td = yaml.load(open(tfl,'r'),Loader=yaml.Loader)
        tlist[ix] = td['Stats']['MAD']

    ax.plot(x_l/100.,tlist,color=cls[idfa],linestyle='-',label=dfad[dfa],linewidth=2)

    if dfa in lbld:
        if len(lbld[dfa]) == 2:
            ax.annotate(dfad[dfa],(lbld[dfa][0],lbld[dfa][1]),\
                color=cls[idfa], fontsize=fsz)
        elif len(lbld[dfa]) == 4:
            ax.annotate(dfad[dfa],(lbld[dfa][2],lbld[dfa][3]), \
                xytext =(lbld[dfa][0],lbld[dfa][1]), color=cls[idfa], \
                fontsize=fsz, arrowprops={'color': cls[idfa], 'arrowstyle':'->', 'lw':2})

ax.set_xlabel('$x$',fontsize=fsz)

xtl = []
for x in x_l:
    if x == 0:
        xtl.append('SCAN')
    elif x == 25:
        xtl.append('SCAN0')
    else:
        xtl.append('{:.2f}'.format(x/100.))

ax.set_xlim(min(x_l)/100.,max(x_l)/100.)
ax.set_ylim(0.0,1.02*max([scf_mads[dfa] for dfa in dfad]))
ax.tick_params(axis='both',labelsize=fsz)
ax.set_xticks(x_l/100.)
ax.set_xticklabels(xtl)

ax.set_ylabel('BH76 MAD (kcal/mol)',fontsize=fsz)
plt.title('DFA@SX-$x$ (solid) and DFA@DFA (dotted)',fontsize=fsz)

for idfa,dfa in enumerate(dfad):
    if scf_mads[dfa] > 0.0:
        ax.hlines(scf_mads[dfa],*ax.get_xlim(),color=cls[idfa],linewidth=2,\
            linestyle=':')

ax.yaxis.set_minor_locator(MultipleLocator(1))
ax.yaxis.set_major_locator(MultipleLocator(5))

#plt.show() ; exit()
plt.savefig('./scan_hyb_mads_one_panel.pdf',dpi=1000,bbox_inches='tight')
