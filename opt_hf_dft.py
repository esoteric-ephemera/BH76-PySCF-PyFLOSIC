import numpy as np
from os import path
import yaml
from scipy.optimize import minimize,minimize_scalar
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

from dens_sens import dfas, dens_sens_analysis, out_dir, scf_dir

def sigmoid(x,x0,x1):

    y = (x - x0)/(x1 - x0)
    cl = [-35., 84., -70., 20.]

    sig = np.zeros(x.shape)
    tmsk = (0.0 <= y) & (y <= 1.0)
    ym = y[tmsk]
    sig[tmsk] = -ym**4*(cl[0] + ym*(cl[1] + ym*(cl[2] + ym*cl[3])))

    sig[y > 1.0] = 1.0

    return sig

def step(x,x0):

    theta = np.zeros(x.shape)
    theta[x > x0] = 1.0
    return theta

def get_mad_dc_dft(dfa,mxpr_l):

    infl = out_dir + dfa + '_BH76_dens_sens.yaml'
    ds_d = yaml.load(open(infl,'r'), Loader=yaml.Loader)

    infl = scf_dir + dfa + '_BH76/BH76_total.yaml'
    if dfa == 'SCAN0':
        infl = './scan_hybs/S25X/SCAN_25_EXX_BH76/BH76_total.yaml'
    scf_d = yaml.load(open(infl,'r'), Loader=yaml.Loader)

    infl = scf_dir + dfa + '@HF_BH76/BH76_total.yaml'
    athf_d = yaml.load(open(infl,'r'), Loader=yaml.Loader)

    Nmx = mxpr_l.shape[0]
    res_d = {
        'BH76': {'MD': np.zeros(Nmx), 'MAD': np.zeros(Nmx), 'SCF MAD': scf_d['Stats']['MAD']},
        'BH76RC': {'MD': np.zeros(Nmx), 'MAD': np.zeros(Nmx), 'SCF MAD': scf_d['RC Stats']['MAD']}
    }

    for aset in ['BH76','BH76RC']:

        if aset == 'BH76':
            optstr = 'RX'
        elif aset == 'BH76RC':
            optstr = 'RC'

        ds_l = []
        scf_l = []
        hf_l = []

        for indx in ds_d[aset]:
            ds_l.append(abs(ds_d[aset][indx]))
            scf_l.append(scf_d[optstr][indx]['Error'])
            hf_l.append(athf_d[optstr][indx]['Error'])

        ds_l = np.asarray(ds_l)
        scf_l = np.asarray(scf_l)
        hf_l = np.asarray(hf_l)

        npts = 1.0*ds_l.shape[0]

        for i in range(Nmx):
            devs = scf_l + step(ds_l,mxpr_l[i])*(hf_l - scf_l)

            res_d[aset]['MD'][i] = np.sum(devs)/npts
            res_d[aset]['MAD'][i] = np.sum(np.abs(devs))/npts

    return res_d

def get_opt_mix():

    ds_l = []
    scf_l = []
    hf_l = []
    for dfa in dfas:

        infl = out_dir + dfa + '_BH76_dens_sens.yaml'
        if not path.isfile(infl):
            dens_sens_analysis()
        ds_d = yaml.load(open(infl,'r'), Loader=yaml.Loader)

        infl = scf_dir + dfa + '_BH76/BH76_total.yaml'
        if dfa == 'SCAN0':
            infl = './scan_hybs/S25X/SCAN_25_EXX_BH76/BH76_total.yaml'
        scf_d = yaml.load(open(infl,'r'), Loader=yaml.Loader)

        infl = scf_dir + dfa + '@HF_BH76/BH76_total.yaml'
        athf_d = yaml.load(open(infl,'r'), Loader=yaml.Loader)

        for indx in ds_d['BH76']:
            ds_l.append(abs(ds_d['BH76'][indx]))
            scf_l.append(scf_d['RX'][indx]['Error'])
            hf_l.append(athf_d['RX'][indx]['Error'])

        for indx in ds_d['BH76RC']:
            ds_l.append(abs(ds_d['BH76RC'][indx]))
            scf_l.append(scf_d['RC'][indx]['Error'])
            hf_l.append(athf_d['RC'][indx]['Error'])

    ds_l = np.asarray(ds_l)
    scf_l = np.asarray(scf_l)
    hf_l = np.asarray(hf_l)
    npts = 1.0*ds_l.shape[0]
    diff_l = hf_l - scf_l

    def tfunopt(sv):
        sigma = step(ds_l,sv)#sigmoid(ds_l,sv[0],sv[1])
        obj = np.sum(np.abs(scf_l + sigma*diff_l))/npts
        return obj

    res = minimize_scalar(tfunopt,bounds=(1.5,3.0),method='bounded')
    #res = minimize(tfunopt,[2.0,2.5],method='L-BFGS-B',bounds=((0.1,10.0),(0.12,10.0)))
    print(res)

    for dfa in dfas:
        td = get_mad_dc_dft(dfa,0.0)
        print(dfa,td)
    exit()

    ps = [round(y,2) for y in res.x]

    xl = np.linspace(0.0,10.0,5000)
    plt.plot(xl,sigmoid(xl,ps[0],ps[1]))
    plt.show()

    return ps

if __name__ == "__main__":

    fig, ax = plt.subplots(2,figsize=(6,7))
    cl = np.linspace(0.0,15.0,5000)

    label_d = {
        'r2SCAN': 'r$^2$SCAN', 'LCwPBE': 'LC-$\\omega$PBE'
        }

    clist = ['darkblue','darkorange','tab:green','darkred','tab:blue','purple']
    lsl = ['-','--']

    for idfa,dfa in enumerate(dfas):

        j = 0
        if dfa in ['SCAN0','LCwPBE']:
            j = 1

        od = get_mad_dc_dft(dfa,cl)
        lblstr = dfa
        if dfa in label_d:
            lblstr = label_d[dfa]
        ax[0].plot(cl,od['BH76']['MAD']/od['BH76']['SCF MAD'],color=clist[idfa],\
            linestyle=lsl[j])
        ax[1].plot(cl,od['BH76RC']['MAD']/od['BH76RC']['SCF MAD'],color=clist[idfa],\
            label=lblstr,linestyle=lsl[j])

    ax[1].legend()
    ax[1].set_xlabel('Density Sensitivity cutoff (kcal/mol)',fontsize=12)
    for i in range(2):
        ax[i].set_xlim(cl.min(),cl.max())
        ax[i].hlines(1.0,cl.min(),cl.max(),color='k',linewidth=1,linestyle=':')
        ax[i].xaxis.set_minor_locator(MultipleLocator(0.5))
        ax[i].xaxis.set_major_locator(MultipleLocator(2.0))
    ax[0].set_ylabel('BH76',fontsize=12)
    ax[1].set_ylabel('BH76RC',fontsize=12)

    ax[0].yaxis.set_minor_locator(MultipleLocator(0.1))
    ax[0].yaxis.set_major_locator(MultipleLocator(.5))
    ax[1].yaxis.set_minor_locator(MultipleLocator(0.05))
    ax[1].yaxis.set_major_locator(MultipleLocator(.1))

    plt.savefig('./density_sensitivity/ds_cutoff.pdf',dpi=600,bbox_inches='tight')
    exit()

    ps = get_opt_mix()
