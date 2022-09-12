import yaml
import numpy as np
from os import path
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

from analysis import BH76_analysis,rdir

x_l = [0, 10, 25, 50, 75, 100]
Nx = len(x_l)

scan_dir = './scan_hybs/'
hf_dir = './results_aug-cc-pvqz/'

dfad = {'LSDA': 'LSDA', 'PBE': 'PBE', 'SCAN': 'SCAN', 'r2SCAN': 'r$^2$SCAN',
    'LCwPBE': 'LC-$\\omega$PBE', 'SX': 'SX-$x$'
}

lbld = {'LSDA': [0,0.50,12.7],
    'PBE': [0,0.11,8.7],#[1,40,.55],
    'LCwPBE': [1, 0.65,1.5],
    'SX': [0,0.10,3.6], 'r2SCAN': [0,0.35,4.45], 'SCAN': [1,0.60,.8]#,'PBE':
}

cls = ['darkblue','darkorange','tab:green','darkred','k','gray']
lsls = ['-','--','-.',':']

fig,ax = plt.subplots(2,1,figsize=(6,7))

ndfa = 1
for dfa in dfad:
    #if dfa == 'SX':
    #    continue
    ndfa += 1

ovals = np.zeros((Nx+2,ndfa))
for ix in range(Nx):
    ovals[ix,0] = x_l[ix]/100.

scf_mads = {'SX': 0.0}
for dfa in dfad:
    if dfa == 'SX':
        continue
    tdir = hf_dir + '/' + dfa + '/' + dfa + '_BH76/'
    tfl = tdir + 'BH76_total.yaml'
    if not path.isfile(tfl):
        BH76_analysis(cdir=tdir)

    td = yaml.load(open(tfl,'r'),Loader=yaml.Loader)
    scf_mads[dfa] = td['Stats']['MAD']

at_hf_mads = {'SX': 1.0}
for dfa in dfad:
    if dfa == 'SX':
        continue
    tdir = hf_dir + '/' + dfa + '/' + dfa + '@HF_BH76/'
    tfl = tdir + 'BH76_total.yaml'
    if not path.isfile(tfl):
        BH76_analysis(cdir=tdir)

    td = yaml.load(open(tfl,'r'),Loader=yaml.Loader)
    at_hf_mads[dfa] = td['Stats']['MAD']

for idfa,dfa in enumerate(dfad):

    tlist = [0.0 for ix in range(Nx)]
    for ix,x in enumerate(x_l):
        if dfa == 'SX':
            if x == 0:
                tdir = scan_dir + '/S00X/SCAN_BH76/'
            elif x == 'HF':
                tdir = hf_dir + '/HF_BH76/'
            else:
                tdir = scan_dir + '/S{:}X/'.format(int(x)) + 'SCAN_{:}_EXX_BH76/'.format(int(x))
        else:
            if x == 0:
                if dfa == 'SCAN':
                    tdir = scan_dir + '/S00X/SCAN_BH76/'
                else:
                    tdir = scan_dir + '/S00X/'+ dfa + '@SCAN_BH76/'
            elif x == 'HF':
                tdir = hf_dir + '/' + dfa + '@HF_BH76/'
            else:
                tdir = scan_dir + '/S{:}X/'.format(int(x)) + dfa + '@S{:}X_BH76/'.format(int(x))

        tfl = tdir + 'BH76_total.yaml'
        if not path.isfile(tfl):
            BH76_analysis(cdir=tdir)

        td = yaml.load(open(tfl,'r'),Loader=yaml.Loader)
        tlist[ix] = td['Stats']['MAD']

        #if dfa != 'SX':
        ovals[ix,1+idfa] = td['Stats']['MAD']

    ovals[-2,1+idfa] = scf_mads[dfa]
    if dfa != 'SX':
        ovals[-1,1+idfa] = at_hf_mads[dfa]

    fxl = np.array(x_l)/100.
    ax[0].plot(fxl,tlist,color=cls[idfa],linestyle=lsls[idfa%len(lsls)],label=dfad[dfa])
    if dfa != 'SX':
        ax[1].plot(fxl,np.array(tlist)/scf_mads[dfa],color=cls[idfa],\
            linestyle=lsls[idfa%len(lsls)],label=dfad[dfa])

    if dfa in lbld:
        ax[lbld[dfa][0]].annotate(dfad[dfa],(lbld[dfa][1],lbld[dfa][2]),\
            color=cls[idfa], fontsize=12)

ax[1].set_xlabel('$x$',fontsize=12)

xtl = []
for ix,x in enumerate(x_l):
    if x == 0:
        xtl.append('SCAN')
    elif x == 25:
        xtl.append('0.25/SCAN0')
    else:
        xtl.append('{:}'.format(fxl[ix]))

for i in range(2):
    ax[i].set_xlim(min(fxl),max(fxl))
    ax[i].tick_params(axis='both',labelsize=12)
    ax[i].set_xticks(fxl)
    ax[i].set_xticklabels(xtl)

ax[0].set_ylabel('BH76 MAD \n(kcal/mol)',fontsize=12,rotation=45)
ax[0].yaxis.set_label_coords(-0.085,.65)
ax[1].set_ylabel('$\\frac{\\mathrm{MAD(DFA@SX)}}{\\mathrm{MAD(DFA@DFA)}}$',\
    fontsize=12,rotation=45)
ax[1].yaxis.set_label_coords(-0.05,0.9)
#ax[1].legend(fontsize=12,ncol=2)
ax[1].hlines(1.0,*ax[1].get_xlim(),color='k',linewidth=1,linestyle=':')

ax[0].yaxis.set_minor_locator(MultipleLocator(1))
ax[0].yaxis.set_major_locator(MultipleLocator(5))
ax[1].yaxis.set_minor_locator(MultipleLocator(0.1))
ax[1].yaxis.set_major_locator(MultipleLocator(0.5))

ax[0].annotate('(a)',(0.01,12),fontsize=16)
ax[1].annotate('(b)',(0.01,1.5),fontsize=16)

#plt.show() ; exit()
plt.savefig('./scan_hyb_mads.pdf',dpi=600,bbox_inches='tight')

tmpstr = 'x'
for dfa in dfad:
    if dfa != 'SX':
        tmpstr += ', {:}'.format(dfa)
tmpstr += '\n'
tmpstr += 'SCF' + (', {:}'*(ndfa-1)).format(*ovals[-2,1:]) + '\n'
for ix in range(Nx):
    tmpstr += ('{:}, '*(ndfa-1) + '{:} \n').format(*ovals[ix])
tmpstr += '@HF' + (', {:}'*(ndfa-1)).format(*ovals[-1,1:]) + '\n'

with open('./DFA_at_SXx.csv','w+') as tfl:
    tfl.write(tmpstr)

tmpstr = '$x$'
for dfa in dfad:
    #if dfa != 'SX':
    tmpstr += ' & {:}'.format(dfa)
tmpstr += ' \\\\ \n'
tmpstr += 'SCF' + (' & {:.2f}'*(ndfa-1)).format(*ovals[-2,1:]) + ' \\\\ \hline \n'
for ix in range(Nx):
    lchar = ''
    if ix == Nx-1:
        lchar = '\hline'
    tmpstr += ('{:.2f}' + ' & {:.2f}'*(ndfa-1) + ' \\\\ {:} \n').format(*ovals[ix],lchar)
tmpstr += '@HF' + (' & {:.2f}'*(ndfa-1)).format(*ovals[-1,1:]) + ' \\\\ \n'

with open('./DFA_at_SXx.tex','w+') as tfl:
    tfl.write(tmpstr)
