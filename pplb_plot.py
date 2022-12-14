import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import numpy as np
import yaml

from analysis import BH76_analysis
from index_labels import indx_to_name

bh76_labs, _ = indx_to_name()

xl = np.array([0, 10, 25, 50, 75, 100])#, 85, 100])
Nx = xl.shape[0]

cte = {
    # CH_3 + ClF ---> CH_3F + Cl
    11: np.array([-0.0702,-0.0645,-0.0539,-0.0302,-0.0013, 0.0270]),#, 0.0105, 0.0270]),
    # H + HCl ---> H_2 + Cl
    39: np.array([-0.02085,-0.01953,-0.01741,-0.01349,-0.00902, -0.00391]),#,-0.00706,-0.00391]),
    # F + H_2 ---> H + HF
    55: np.array([-0.060588583,-0.050148583,-0.0332,-0.004808583, 0.018791417, 0.035881417]),#, \
#        0.026381417, 0.035881417])
}

scanvals = {}
for arx in cte:
    scanvals[arx] = cte[arx][0]
    cte[arx] /= abs(scanvals[arx])

xbds = [1e20,-1e20]
for arx in cte:
    xbds = [min(xbds[0],cte[arx].min()),max(xbds[1],cte[arx].max())]

wdir = './scan_hybs/'
fnls = ['SCAN']
colors = ['darkblue','darkorange','tab:green']
linestyles = ['-',':']

wd = {}
for arx in cte:
    for irx in range(2):
        wd[arx+irx] = {}
        for afnl in fnls:
            wd[arx+irx][afnl] = np.zeros(Nx)

ybds = [1e20,-1e20,1e20,-1e20]
for ix, ax in enumerate(xl):

    hstr = 'S'+'{:}'.format(ax).zfill(2) + 'X'
    tdir = wdir + hstr + '/'
    if ax == 0:
        atdfa = 'SCAN'
    else:
        atdfa = hstr

    for afnl in fnls:

        #td = BH76_analysis(cdir=tdir + '{:}@{:}_BH76/'.format(afnl,atdfa))
        if ax == 0 and afnl == 'SCAN':
            fnldir = 'SCAN'
        else:
            fnldir = '{:}@{:}'.format(afnl,atdfa)
        td = yaml.load(open(tdir+fnldir + '_BH76/BH76_total.yaml','r'),\
            Loader=yaml.Loader)
        for arx in cte:
            for irx in range(2):
                wd[arx + irx][afnl][ix] = td['RX'][arx + irx]['Error']
                ybds[0] = min(ybds[0],wd[arx + irx][afnl][ix])
                ybds[1] = max(ybds[1],wd[arx + irx][afnl][ix])

fig, ax = plt.subplots(figsize=(6,4))

pstr = 'Reaction Index, DFA, c0, c1, c2, \Delta N -> E_min \n'

for irx, arx in enumerate(cte):

    for jrx in range(2):
        for ifnl, afnl in enumerate(fnls):

            lstr = None
            if ifnl == 0 and jrx == 0:
                lstr = bh76_labs[arx] + '  [{:.4f}]'.format(scanvals[arx])

            amat = np.zeros((Nx,3))
            amat[:,0] = 1.
            amat[:,1] = cte[arx]
            amat[:,2] = cte[arx]**2

            cfs = np.linalg.lstsq(amat, wd[arx + jrx][afnl], rcond=None)[0]
            dnmin = -cfs[1]/(2.*cfs[2])*abs(scanvals[arx])
            print(arx+jrx,afnl,'{:.4f}'.format(dnmin))
            pstr += ('{:}, '*5 + '{:} \n').format(arx+jrx,afnl,*cfs,dnmin)

            itp = np.linspace(cte[arx].min(),cte[arx].max(),2000)
            ax.scatter(cte[arx],wd[arx + jrx][afnl],color=colors[irx])
            ax.plot(itp,cfs[0] + itp*(cfs[1] + itp*cfs[2]),color=colors[irx],\
                linestyle=linestyles[jrx],label=lstr)
            #ax[jrx].plot(cte[arx],wd[arx + jrx][afnl],color=colors[irx],\
            #    linestyle=linestyles[ifnl],marker='o',label=lstr)

with open('./PPLB_pars.csv','w+') as tfl:
    tfl.write(pstr)

ax.set_xlabel('$\\Delta N /|\\Delta N(\\mathrm{SCAN})|$',fontsize=12)
ax.set_ylabel('SCAN@SX-$x$ \nBH error (kcal/mol)',fontsize=12)
#ax[1].set_ylabel('REV BH error (kcal/mol)',fontsize=12)

ax.set_xlim(xbds)
ax.set_ylim(ybds[0],ybds[1])

ax.hlines(0.,*xbds,color='k',linewidth=1)
ax.vlines(0.,ybds[0],ybds[1],color='k',linewidth=1)

ax.xaxis.set_minor_locator(MultipleLocator(0.1))
ax.xaxis.set_major_locator(MultipleLocator(0.2))

ax.yaxis.set_minor_locator(MultipleLocator(1.))
ax.yaxis.set_major_locator(MultipleLocator(2.))

ax.legend(fontsize=10,loc=(0.01,.7),ncol=1,frameon=False,\
    title='Reaction [$\Delta N(\\mathrm{SCAN})$]')
#plt.show() ; exit()
plt.savefig('./scan_hyb_pplb.pdf',dpi=600,bbox_inches='tight')
