import yaml
import numpy as np
from os import path
from glob import glob
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

from analysis import BH76_analysis,rdir


sl_name = 'r2SCAN'
sl_tex_name = 'r$^2$SCAN'
abbr_name = 'R2'
hf_dir = './results_aug-cc-pvqz/'

hyb_dir = './{:}_hybrids/'.format(sl_name)

fl_l = glob('./{:}_hybrids/{:}*X'.format(sl_name,sl_name))
x_l = np.sort(np.asarray([float(adir[:-1].split(sl_name)[-1]) for adir in fl_l]))
Nx = len(x_l)

dfad = {'LSDA': 'LSDA', 'PBE': 'PBE', 'BLYP': 'BLYP' ,
    'SCAN': 'SCAN', 'r2SCAN': 'r$^2$SCAN',
    'B3LYP': 'B3LYP', 'LCwPBE': 'LC-$\\omega$PBE',
    abbr_name: '{:}X-$x$'.format(abbr_name)
}

lbld = {'LSDA': [0,0.70,11],
    'PBE': [0,0.60,6.4],
    'BLYP': [1,0.50,.57],
    'B3LYP': [0,0.41,4.1],
    'LCwPBE': [1, 0.65, 1.38],
    'R2': [0,0.80,6.1], 'r2SCAN': [0,0.25,5.4], 'SCAN': [1,0.60,.8]
}

cls = ['darkblue','darkorange','dodgerblue','tab:green','darkred','indigo','k','gray']
lsls = ['-','--','-.',':']

fig,ax = plt.subplots(2,1,figsize=(6,7))

ndfa = 1
for dfa in dfad:
    ndfa += 1

ovals = np.zeros((Nx+2,ndfa))
for ix in range(Nx):
    ovals[ix,0] = x_l[ix]/100.

scf_mads = {abbr_name: 0.0}
for dfa in dfad:
    if dfa == abbr_name:
        continue
    tdir = hf_dir + '/' + dfa + '/' + dfa + '_BH76/'
    tfl = tdir + 'BH76_total.yaml'
    if not path.isfile(tfl):
        BH76_analysis(cdir=tdir)

    td = yaml.load(open(tfl,'r'),Loader=yaml.Loader)
    scf_mads[dfa] = td['Stats']['MAD']

at_hf_mads = {abbr_name: 1.0}
for dfa in dfad:
    if dfa == abbr_name:
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
        if dfa == abbr_name:
            if x == 0:
                tdir = hyb_dir + '/{:}00X/{:}_BH76/'.format(sl_name,sl_name)
            elif x == 'HF':
                tdir = hf_dir + '/HF_BH76/'
            else:
                tdir = hyb_dir + '/{:}{:}X/{:}_{:}_EXX_BH76/'.format(\
                    sl_name,int(x),sl_name,int(x))
        else:
            if x == 0:
                if dfa == sl_name:
                    tdir = hyb_dir + '/{:}00X/{:}_BH76/'.format(sl_name,sl_name)
                else:
                    tdir = hyb_dir + '/{:}00X/{:}@{:}00X_BH76/'.format(\
                        sl_name,dfa,sl_name)
            elif x == 'HF':
                tdir = hf_dir + '/' + dfa + '@HF_BH76/'
            else:
                tdir = hyb_dir + '/{:}{:}X/{:}@{:}{:}X_BH76/'.format(\
                    sl_name,int(x),dfa,sl_name,int(x))

        tfl = tdir + 'BH76_total.yaml'
        if not path.isfile(tfl):
            BH76_analysis(cdir=tdir)

        td = yaml.load(open(tfl,'r'),Loader=yaml.Loader)
        tlist[ix] = td['Stats']['MAD']

        ovals[ix,1+idfa] = td['Stats']['MAD']

    ovals[-2,1+idfa] = scf_mads[dfa]
    if dfa != abbr_name:
        ovals[-1,1+idfa] = at_hf_mads[dfa]

    ax[0].plot(x_l/100., tlist, color=cls[idfa], linestyle=lsls[idfa%len(lsls)],\
        label=dfad[dfa])
    if dfa != abbr_name:
        ax[1].plot(x_l/100.,np.array(tlist)/scf_mads[dfa],color=cls[idfa],\
            linestyle=lsls[idfa%len(lsls)],label=dfad[dfa])

    if dfa in lbld:
        ax[lbld[dfa][0]].annotate(dfad[dfa],(lbld[dfa][1],lbld[dfa][2]),\
            color=cls[idfa], fontsize=12)

ax[1].set_xlabel('$x$',fontsize=12)

xticks = np.array([0,10,25,40,50,60,70,80,90,100])/100.
xtl = []
for x in xticks:
    if x == 0.:
        xtl.append(sl_tex_name)
    elif x == 0.25:
        xtl.append('0.25/{:}0'.format(sl_tex_name))
    else:
        xtl.append('{:}'.format(x))

for i in range(2):
    ax[i].set_xlim(min(x_l)/100.,max(x_l)/100.)
    ax[i].tick_params(axis='both',labelsize=10)
    ax[i].set_xticks(xticks)
    ax[i].set_xticklabels(xtl,rotation=0)

ax[0].set_ylabel('BH76 MAD \n(kcal/mol)',fontsize=12,rotation=45)
ax[0].yaxis.set_label_coords(-0.085,.65)
ax[1].set_ylabel('$\\frac{\\mathrm{MAD(DFA@'+abbr_name+'X)}}{\\mathrm{MAD(DFA@DFA)}}$',\
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
plt.savefig('./{:}_hyb_mads.pdf'.format(sl_name),dpi=600,bbox_inches='tight')

tmpstr = 'x'
for dfa in dfad:
    tmpstr += ', {:}'.format(dfa)
    if dfa == abbr_name:
        tmpstr += 'X-x'
tmpstr += '\n'
tmpstr += 'SCF' + (', {:}'*(ndfa-1)).format(*ovals[-2,1:]) + '\n'
for ix in range(Nx):
    tmpstr += ('{:}, '*(ndfa-1) + '{:} \n').format(*ovals[ix])
tmpstr += '@HF' + (', {:}'*(ndfa-1)).format(*ovals[-1,1:]) + '\n'

with open('./DFA_at_{:}Xx.csv'.format(abbr_name),'w+') as tfl:
    tfl.write(tmpstr)

tmpstr = '$x$'
for dfa in dfad:
    tmpstr += ' & {:}'.format(dfa)
    if dfa == abbr_name:
        tmpstr += 'X-$x$'
tmpstr += ' \\\\ \n'
tmpstr += 'SCF' + (' & {:.2f}'*(ndfa-1)).format(*ovals[-2,1:]) + ' \\\\ \hline \n'
for ix in range(Nx):
    lchar = ''
    if ix == Nx-1:
        lchar = '\hline'
    tmpstr += ('{:.2f}' + ' & {:.2f}'*(ndfa-1) + ' \\\\ {:} \n').format(*ovals[ix],lchar)
tmpstr += '@HF' + (' & {:.2f}'*(ndfa-1)).format(*ovals[-1,1:]) + ' \\\\ \n'

with open('./DFA_at_{:}Xx.tex'.format(abbr_name),'w+') as tfl:
    tfl.write(tmpstr)
