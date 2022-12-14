import yaml
from os import path

rdir = path.dirname(path.realpath(__file__)) + '/'

refs = yaml.load(open(rdir + 'BH76_ref_energies.yaml','r'), Loader=yaml.Loader)
ts_d = {}
for ibh, abh in enumerate(refs):
    for asys in refs[abh]['Stoich']:
        if refs[abh]['Stoich'][asys] > 0:
            if asys not in ts_d:
                ts_d[asys] = [abh]
            else:
                ts_d[asys].append(abh)

def s2_to_s(ssq):
    return (-1. + (1. + 4*ssq)**(0.5))/2.

wpath = './results_aug-cc-pvqz/HF/HF_BH76/'

ostr = 'Transition state, Computed <S^2>, Computed S, Input S , PE (%), BH index 1, BH index 2,...\n'

for asys in ts_d:

    got_two_s = False
    got_res = False
    with open(wpath + asys + '/inp.txt','r') as tfl:
        for arow in tfl:
            if '2S' in arow:
                ideal_s = float(arow.split('=')[-1].strip())/2.
                got_two_s = True
            elif 'restricted' in arow:
                tmp = arow.split('=')[-1].strip().lower()
                if tmp[0] == 't' or tmp == 'true':
                    LR = True
                elif tmp[0] == 'f' or tmp == 'false':
                    LR = False
                got_res = True

            if got_two_s and got_res:
                break

    if LR:
        ostr += '{:}, 0, 0, 0, 0'.format(asys)

    else:
        with open(wpath + asys + '/' + asys + '.txt','r') as tfl:
            for arow in tfl:
                if '<S^2>' in arow:
                    ssq = float(arow.split('<S^2> =')[-1].split('2S+1')[0])
                    comps = s2_to_s(ssq)
                    ostr += '{:}, {:.4f}, {:.4f}, {:}, {:.4f}'\
                        .format(asys,ssq,comps,ideal_s, 100*(comps/ideal_s - 1.))
                    break

    for anindx in ts_d[asys]:
        ostr += ', {:}'.format(anindx)
    ostr += '\n'


with open('./spin_contam.csv','w+') as tfl:
    tfl.write(ostr)
