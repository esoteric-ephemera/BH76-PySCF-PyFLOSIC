import yaml
from os import path,system
from analysis import BH76_analysis

eH_to_eV = 27.211386245988
eV_to_kcalmol = 23.06054783061903
eH_to_kcalmol = eH_to_eV*eV_to_kcalmol

rdir = path.dirname(path.realpath(__file__)) + '/'

dfas = ['LSDA','PBE', 'BLYP', 'SCAN', 'r2SCAN', 'M06-L', 'MN15-L', 'SCAN0', 'B3LYP','LCwPBE']

scf_dir = './results_aug-cc-pvqz/'
out_dir = './density_sensitivity/'
if not path.isdir(out_dir):
    system('mkdir '+out_dir)

def dens_sens_analysis():

    refs = yaml.load(open(rdir + 'BH76_ref_energies.yaml','r'), Loader=yaml.Loader)
    cmd = yaml.load(open(rdir + 'BH76_chg_2s.yaml','r'), Loader=yaml.Loader)

    ts_l = []
    rp_l = []
    for arx in refs:
        for asys in refs[arx]['Stoich']:
            if (refs[arx]['Stoich'][asys] > 0) and (asys not in ts_l):
                ts_l.append(asys)
            elif (refs[arx]['Stoich'][asys] < 0) and (asys not in rp_l):
                rp_l.append(asys)

    Nts = len(ts_l)
    Nrp = len(rp_l)
    Ntot = Nts + Nrp

    metrics = ['MDS','MADS','STDDEV','AMIN','AMAX']

    subsets = ['RP','TS','Total','BH76','BH76RC']

    tstr = {'RP': 'DFA, MDS, MADS, STDDEV, AMIN, AMAX \n'}
    tex_str = {'RP': 'DFA & ADS & AMIN & AMAX \\\\ \\hline \n'}
    for subset in subsets:
        if subset not in tstr:
            tstr[subset] = ''
        tstr[subset] += subset
        for amet in metrics:
            tstr[subset] += ', '
        tstr[subset] += '\n'

        if subset not in tex_str:
            tex_str[subset] = ''
        tex_str[subset] += '\\multicolumn{4}{c}{\\textit{'+subset+'}} \\\\ \n'

    for idfa,dfa in enumerate(dfas):

        if dfa == 'LSDA':
            lsd_dir = scf_dir + dfa + '/LSDA_BH76/'
        else:
            lsd_dir = scf_dir + dfa + '/{:}@LSDA_BH76/'.format(dfa)
        BH76_analysis(cdir=lsd_dir)
        hf_dir = scf_dir + '{:}/{:}@HF_BH76/'.format(dfa,dfa)
        BH76_analysis(cdir=hf_dir)

        at_lsd = yaml.load(open(lsd_dir+'BH76_total.yaml','r'), Loader=yaml.Loader)
        at_hf = yaml.load(open(hf_dir+'BH76_total.yaml','r'), Loader=yaml.Loader)

        mds = {}
        mads = {}
        var = {}
        amin = {}
        amax = {}
        for subset in subsets:
            mds[subset] = 0.0
            mads[subset] = 0.0
            var[subset] = 0.0
            amin[subset] = 1.e20
            amax[subset] = -1.e20

        dsensd = {'SP': {} , 'BH76': {}, 'BH76RC': {} }
        for asys in ts_l + rp_l:
            s = (at_lsd['SP'][asys] - at_hf['SP'][asys])*eH_to_kcalmol
            dsensd['SP'][asys] = s

            if asys in ts_l:
                tkey = 'TS'
            elif asys in rp_l:
                tkey = 'RP'
            for akey in [tkey,'Total']:
                mds[akey] += s
                mads[akey] += abs(s)
                var[akey] += s**2
                amin[akey] = min(amin[akey],abs(s))
                amax[akey] = max(amax[akey],abs(s))

        for akey in ['RX','RC']:

            if akey == 'RX':
                didx = 'BH76'
            elif akey == 'RC':
                didx = 'BH76RC'

            for asys in at_lsd[akey]:
                s = (at_lsd[akey][asys]['Energy'] - at_hf[akey][asys]['Energy'])
                dsensd[didx][asys] = s

                mds[didx] += s
                mads[didx] += abs(s)
                var[didx] += s**2
                amin[didx] = min(amin[didx],abs(s))
                amax[didx] = max(amax[didx],abs(s))

        norm = {
            'Total': 1.0/Ntot, 'RP': 1.0/Nrp, 'TS': 1.0/Nts,
            'BH76': 1.0/len(at_lsd['RX'].keys()),
            'BH76RC': 1.0/len(at_lsd['RC'].keys())
            }

        dsensd['Stats'] = {}
        for istat,astat in enumerate(subsets):
            ms = mds[astat]*norm[astat]
            mas = mads[astat]*norm[astat]
            vars = var[astat]*norm[astat]

            stddev = max(0.0,vars - mas**2)**(0.5)
            dsensd['Stats'][astat] = { 'MDS': ms, 'MADS': mas, 'VAR': vars, 'STDDEV': stddev,
                'AMIN': amin[astat], 'AMAX': amax[astat]}

            tstr[astat] += ('{:}' + ', {:.2f}'*5 + '\n').format(dfa,ms,mas,stddev,amin[astat],amax[astat])
            lchar = ''
            if idfa == len(dfas)-1:
                lchar = '\\hline'
            tex_str[astat]+= ('{:} & ${:.2f} \\pm {:.2f}$ & {:.2f} & {:.2f} \\\\ {:} \n').format(\
                dfa,mas,stddev,amin[astat],amax[astat],lchar)

        ofl = out_dir + dfa + '_BH76_dens_sens.yaml'
        with open(ofl,'w+') as tfl:
            yaml.dump(dsensd,tfl,Dumper=yaml.Dumper)

    with open(out_dir+'DS_analysis.csv','w+') as tfl:
        for akey in tstr:
            tfl.write(tstr[akey])

    with open(out_dir+'DS_analysis.tex','w+') as tfl:
        for akey in tex_str:
            tfl.write(tex_str[akey])

    return

if __name__ == "__main__":

    dens_sens_analysis()
