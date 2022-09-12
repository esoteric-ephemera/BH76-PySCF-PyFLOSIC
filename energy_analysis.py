import yaml
from os import path
from analysis import BH76_analysis

eH_to_eV = 27.211386245988
eV_to_kcalmol = 23.06054783061903
eH_to_kcalmol = eH_to_eV*eV_to_kcalmol

rdir = path.dirname(path.realpath(__file__)) + '/'
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

def get_energy_shifts(dfa):

    dets = [0., 0.]
    derp = [0., 0.]

    scf_d = yaml.load(open('./{:}/{:}_BH76/BH76_total.yaml'.format(dfa,dfa),'r'),\
        Loader=yaml.Loader)
    hf_d = yaml.load(open('./{:}/{:}@HF_BH76/BH76_total.yaml'.format(dfa,dfa),'r'),\
        Loader=yaml.Loader)

    for isys, asys in enumerate(scf_d['SP']):
        tdiff = hf_d['SP'][asys] - scf_d['SP'][asys]

        if tdiff < 0.:
            wstr = 'Warning: something wrong with system = {:} for DFA = {:}\n'.format(asys,dfa)
            wstr += 'SCF energy higher in energy than HF energy'
            print(wstr)

        if asys in ts_l:
            dets[0] += tdiff
            dets[1] += tdiff**2
        elif asys in rp_l:
            derp[0] += tdiff
            derp[1] += tdiff**2

    dets[0] *= eH_to_kcalmol/Nts
    dets[1] = max(0.,dets[1]*eH_to_kcalmol**2/Nts - dets[0]**2)**(0.5)

    derp[0] *= eH_to_kcalmol/Nrp
    derp[1] = max(0.,derp[1]*eH_to_kcalmol**2/Nrp - derp[0]**2)**(0.5)

    return dets, derp

def make_eshift_table():
    dfas = ['LSDA', 'PBE', 'BLYP', 'SCAN', 'r2SCAN', 'M06-L', 'MN15-L','B3LYP', 'SCAN0', 'LCwPBE']
    conv_d = {'r2SCAN': 'r$^2$SCAN', 'LCwPBE': 'LC-$\\omega$PBE'}

    tstr = 'DFA & TS & RP & Difference \\\\ \\hline \n'
    for dfa in dfas:
        if dfa in conv_d:
            tstr += conv_d[dfa]
        else:
            tstr += dfa

        dEts, dErp = get_energy_shifts(dfa)
        tstr += ' & ${:.2f} \\pm {:.2f}$ & ${:.2f} \\pm {:.2f}$ & {:.2f} \\\\ \n'.format(*dEts, *dErp, dEts[0] - dErp[0])

    with open('./energy_shifts.tex','w+') as tfl:
        tfl.write(tstr)
    return

if __name__ == "__main__":

    make_eshift_table()
