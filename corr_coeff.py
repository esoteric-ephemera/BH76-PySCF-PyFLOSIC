import numpy as np
import yaml
from os import listdir, path, system
import matplotlib.pyplot as plt
import matplotlib.cm as mcm
from seaborn import heatmap

from analysis import BH76_analysis
from SP_CSV_to_YAML import conversion

conversion() # ensures FLOSIC errors are recomputed

def str_rep_rules(instr):
    rules = {
        'r2SCAN00X': 'r$^2$SCAN',
        'r2SCAN50X': 'R2X-0.5',
        'r2SCAN': 'r$^2$SCAN',
        'LCwPBE': 'LC-$\\omega$PBE',
        'S50X': 'SX-0.5',
        'S00X': 'SCAN'
    }
    outstr = instr
    for arule in rules:
        outstr = outstr.replace(arule,rules[arule])
    return outstr

wkdir = './corr_coeffs/'
if not path.isdir(wkdir):
    system('mkdir ' + wkdir)

rdir = path.dirname(path.realpath(__file__)) + '/'
bhref = yaml.load(open(rdir + 'BH76_ref_energies.yaml','r'),\
    Loader=yaml.Loader)
Nrx = len(bhref.keys())
ref_l = np.zeros(Nrx)
for irx in bhref:
    ref_l[irx-1] = bhref[irx]['Ref']

def get_cov(xl,yl):
    nelt = 1.*xl.shape[0]

    xm = np.sum(xl)/nelt
    xvar = np.sum(xl**2)/nelt
    xstd = max(0.,xvar - xm**2)**(0.5)

    ym = np.sum(yl)/nelt
    yvar = np.sum(yl**2)/nelt
    ystd = max(0.,yvar - ym**2)**(0.5)

    xycov = np.sum(xl*yl)/nelt - xm*ym

    return xycov/max(1.e-15,xstd*ystd)

dfas = ['LSDA','PBE','SCAN','r2SCAN']
Ndfa = 1.*len(dfas)

scf_dir = './results_aug-cc-pvqz/'
scan_hyb_dir = './scan_hybs/'
r2scan_hyb_dir = './r2SCAN_hybrids/'
flosic_dir = './FLOSIC/'
refs = ['HF','LSDA','PBE','S00X','r2SCAN00X','LCwPBE','S50X','r2SCAN50X','SCAN-FLOSIC']

#for dirl in [scan_hyb_dir,r2scan_hyb_dir]:
#    for adir in listdir(dirl):
#        if path.isdir(dirl + adir) and adir[-1] == 'X':
#            refs.append(adir)

trefl = [str_rep_rules(aref) for aref in refs]
tcmap = mcm.get_cmap('RdBu_r')

def make_plot(mat,flnm):
    fig, ax = plt.subplots(figsize=(8,6))
    heatmap(mat,annot=True,xticklabels=trefl,yticklabels=trefl,cmap=tcmap,\
        vmin=-1., vmax = 1., fmt='.2f')
    plt.yticks(rotation = 60)
    plt.xticks(rotation = 20)


    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    tx = np.linspace(xlim[0],xlim[1],2000)
    slp = (ylim[0] - ylim[1])/(xlim[1] - xlim[0])
    icpt = ylim[1] - xlim[0]*slp
    ax.plot(tx,slp*tx + icpt,color='white',linewidth=3,linestyle='--')
    """
    ax.annotate('FDE',(xlim[0],1.05*ylim[1]),fontsize=12,\
        rotation=-180.*np.arctan(slp)/np.pi)
    ax.annotate('DDE',(xlim[0],0.95*ylim[1]),fontsize=12,\
        rotation=-180.*np.arctan(slp)/np.pi)
    """

    #plt.show() ; exit()
    plt.savefig(flnm,dpi=600,bbox_inches='tight')
    plt.cla()
    plt.clf()
    plt.close()
    return

Nref = len(refs)

avg_covmat = np.zeros((Nref,Nref))
nnorm_coeff = len(dfas)*np.ones((Nref,Nref))

for dfa in dfas:
    covmat = np.zeros((Nref,Nref))

    fde_l = np.zeros((Nref,Nrx))
    dde_l = np.zeros((Nref,Nrx))

    td = scf_dir + '{:}/{:}_BH76'.format(dfa,dfa)
    scf_d = BH76_analysis(cdir=td)

    #trefl = []
    do_ref = [True for iref in range(Nref)]

    for iref, aref in enumerate(refs):

        if aref == dfa or (dfa == 'SCAN' and aref == 'S00X') or \
            (dfa == 'r2SCAN' and aref == 'r2SCAN00X'):
            nscf_d = scf_d.copy()
            do_ref[iref] = False
            #continue

        elif 'FLOSIC' in aref:
            td = flosic_dir + '{:}/{:}@{:}_BH76/'.format(dfa,dfa,aref)
            nscf_d = yaml.load(open(td+'BH76_total.yaml','r'),\
                Loader=yaml.Loader)
        else:

            if aref[0] == 'S' and aref[-1] == 'X':
                if aref == 'S00X':
                    td = scan_hyb_dir + '{:}/{:}@SCAN_BH76/'.format(aref,dfa)
                else:
                    td = scan_hyb_dir + '{:}/{:}@{:}_BH76'.format(aref,dfa,aref)
            elif aref[:6] == 'r2SCAN' and aref[-1] == 'X':
                td = r2scan_hyb_dir + '{:}/{:}@{:}_BH76'.format(aref,dfa,aref)
            else:
                td = scf_dir + '{:}/{:}@{:}_BH76'.format(dfa,dfa,aref)
            nscf_d = BH76_analysis(cdir=td)

        #trefl.append(str_rep_rules(aref))

        for irx in range(Nrx):
            fde_l[iref,irx] = nscf_d['RX'][irx+1]['Energy'] - ref_l[irx]
            dde_l[iref,irx] = scf_d['RX'][irx+1]['Energy'] - nscf_d['RX'][irx+1]['Energy']

    for iref in range(Nref):
        for jref in range(iref,Nref):
            covmat[jref,iref] = get_cov(dde_l[iref,:],dde_l[jref,:])
            covmat[iref,jref] = get_cov(fde_l[iref,:],fde_l[jref,:])

    np.savetxt(wkdir + '{:}_cov_mat.csv'.format(dfa),covmat,delimiter=',')

    make_plot(covmat,wkdir+'{:}_cov_mat.pdf'.format(dfa))

    for iref in range(Nref):
        for jref in range(iref,Nref):
            if iref < jref:
                if do_ref[jref] and do_ref[iref]:
                    avg_covmat[jref,iref] += covmat[jref,iref]
                else:
                    #print(trefl[iref],trefl[jref])
                    nnorm_coeff[jref,iref] -= 1
            avg_covmat[iref,jref] += covmat[iref,jref]

np.savetxt(wkdir + 'norm_coeffs.csv',nnorm_coeff,delimiter=',')
avg_covmat /= nnorm_coeff

np.savetxt(wkdir + 'avg_cov_mat.csv',avg_covmat,delimiter=',')

make_plot(avg_covmat,wkdir+'avg_cov_mat.pdf')
