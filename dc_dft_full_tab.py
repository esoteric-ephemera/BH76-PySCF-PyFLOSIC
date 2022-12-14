import yaml
from os import path,system
from analysis import BH76_analysis
from SP_CSV_to_YAML import conversion
from index_labels import indx_to_name

#yaml.Dumper.ignore_aliases = lambda *args : True

def str_rep_rules(instr):
    rules = {
        'r2SCAN': '\\rrscan{}',
        'LCwPBE': '\\lcwpbe{}',
        'S50X': 'SX-0.5'
    }
    outstr = instr
    for arule in rules:
        outstr = outstr.replace(arule,rules[arule])
    return outstr


tex_str = ''

sets = {}
sets['BH76'], sets['BH76RC'] = indx_to_name()

key_to_set = {'BH76': 'RX', 'BH76RC': 'RC'}

conversion() # ensures FLOSIC errors are recomputed
no_flo_dat = ['LCwPBE','BLYP','B3LYP', 'M06-L','MN15-L']

rdir = path.dirname(path.realpath(__file__)) + '/'

out_dir = './dc_dft_full/'
if not path.isdir(out_dir):
    system('mkdir '+out_dir)

scf_dir = './results_aug-cc-pvqz/'
nrlmol_basis_dir = './results_NRLMOL_cart/'
refs = {
    #'SCAN@HF': './results_aug-cc-pvqz/',
    'LCwPBE': './results_aug-cc-pvqz/',
    'S50X': './scan_hybs/S50X/',
    'SCAN-FLOSIC': './FLOSIC/'
    #'S25X': './scan_hybs/S25X/',
}
Nref = len(refs.keys())

refd = { 'BH76': {}, 'BH76RC': {} }

dfas = ['LSDA','PBE','BLYP','SCAN','r2SCAN','M06-L','MN15-L','B3LYP','LCwPBE','LSDA@HF','PBE@HF','SCAN@HF','r2SCAN@HF']
metrics = ['MFE','MAFE','MDE','MADE']


for aref in refs:
    if aref in ['S25X','S50X']:
        frac = aref[1:-1]
        tdir = refs[aref]+'SCAN_{:}_EXX_BH76/'.format(frac)
    elif aref in ['SCAN@HF']:
        tdir = refs[aref] + 'SCAN/' + aref +'_BH76/'
    else:
        tdir = refs[aref] + aref + '/' + aref +'_BH76/'

    tfl = tdir + 'BH76_total.yaml'
    if aref != 'SCAN-FLOSIC': # conversion() takes care of SCAN-FLOSIC
        BH76_analysis(cdir=tdir)

    td = yaml.load(open(tfl,'r'), Loader=yaml.Loader)
    for aset in key_to_set:
        tkey = key_to_set[aset]
        refd[aset][aref] = {}
        for arx in td[tkey]:
            refd[aset][aref][arx] = td[tkey][arx]['Energy']

for idfa,dfa in enumerate(dfas):

    if idfa > 0:
        tex_str += '\n\\clearpage\n'
    tex_str += '\\subsection{'+str_rep_rules(dfa)+'}\n\n'

    tstr = {'BH76': '', 'BH76RC': ''}
    hstr = {'BH76': '\\\\ \\hline\n & ', 'BH76RC': '\\\\ \\hline\n & '}
    trefs = {}

    wrefs = {
        'BH76': { 'SCF': {}, 'NSCF': {} },
        'BH76RC': { 'SCF': {}, 'NSCF': {} }
    }

    dsplit = dfa.split('@')
    tdfa = dsplit[0]
    if len(dsplit) > 1:
        at_dfa = dsplit[1]

    for aset in hstr:

        for aref in refs:
            if dfa == aref or ((aref == 'SCAN-FLOSIC') and (dfa in no_flo_dat)):
                continue
            sref = str_rep_rules(aref)
            hstr[aset] += ' & \\multicolumn{2}{c}{'+sref+'}'
            trefs[aref] = refs[aref]

        hstr[aset] += ' \\\\ \n'
        hstr[aset] += 'Index & Reaction'

    td_scf = BH76_analysis(cdir=scf_dir+'{:}/{:}_BH76/'.format(tdfa,dfa))
    for aset in key_to_set:
        tkey = key_to_set[aset]
        wrefs[aset]['SCF']['REF'] = {}
        for arx in td_scf[tkey]:
            wrefs[aset]['SCF']['REF'][arx] = td_scf[tkey][arx]['Energy']

    for aref in trefs:

        if aref == 'SCAN-FLOSIC':

            if '@' in dfa:
                tdir_scf = nrlmol_basis_dir+'{:}/{:}_BH76/'.format(tdfa,dfa)
                td_scf = BH76_analysis(cdir=tdir_scf)

            else:
                tdir_scf = './FLOSIC/{:}/{:}_BH76/'.format(dfa,dfa)
                td_nscf = yaml.load(open(tdir_scf+'BH76_total.yaml','r'), Loader=yaml.Loader)

            tdir_nscf = './FLOSIC/{:}/{:}@{:}_BH76/'.format(tdfa,tdfa,aref)
            td_nscf = yaml.load(open(tdir_nscf+'BH76_total.yaml','r'), Loader=yaml.Loader)

        else:

            if aref in ['S25X','S50X']:
                tdir_nscf = refs[aref]+'{:}@{:}_BH76/'.format(tdfa,aref)
            else:
                tdir_nscf = refs[aref]+'{:}/{:}@{:}_BH76/'.format(tdfa,tdfa,aref)
            td_nscf = BH76_analysis(cdir=tdir_nscf)

        for iii in range(2):
            if iii == 0:
                if aref == 'SCAN-FLOSIC':
                    td = td_scf
                else:
                    for aset in key_to_set:
                        wrefs[aset]['SCF'][aref] = wrefs[aset]['SCF']['REF'].copy()
                    continue
                skey = 'SCF'
            elif iii == 1:
                td = td_nscf
                skey = 'NSCF'

            for aset in key_to_set:
                tkey = key_to_set[aset]
                wrefs[aset][skey][aref] = {}
                for arx in td[tkey]:
                    wrefs[aset][skey][aref][arx] = td[tkey][arx]['Energy']

    for aset in tstr:

        nrx = len(sets[aset].keys())

        for aref in trefs:
            hstr[aset] += ' & FE & DE'
        hstr[aset] += ' \\\\ \\hline \n'

        wd = {}
        for aref in trefs:
            wd[aref] = {
                'MDE': 0., 'MADE': 0., 'VARDE': 0.,
                'MFE': 0., 'MAFE': 0., 'VARFE': 0.
            }

        for irx, arx in enumerate(sets[aset]):
            tstr[aset] += '{:} & {:}'.format(arx,sets[aset][arx])

            for aref in trefs:

                fe = wrefs[aset]['NSCF'][aref][arx] - refd[aset][aref][arx]
                de = wrefs[aset]['SCF'][aref][arx] - wrefs[aset]['NSCF'][aref][arx]
                tstr[aset] += ' & {:.2f} & {:.2f}'.format(fe,de)

                wd[aref]['MDE'] += de
                wd[aref]['MADE'] += abs(de)
                wd[aref]['VARDE'] += de**2

                wd[aref]['MFE'] += fe
                wd[aref]['MAFE'] += abs(fe)
                wd[aref]['VARFE'] += fe**2

            lchar = ''
            if irx == nrx - 1:
                lchar = '\\hline'
            tstr[aset] += ' \\\\ {:} \n'.format(lchar)

        for aref in trefs:
            for akey in wd[aref]:
                wd[aref][akey] /= (1.*nrx)

        tstr[aset] += ' & Mean'
        for aref in trefs:
            tstr[aset] += ' & {:.2f} & {:.2f}'\
                .format(wd[aref]['MFE'], wd[aref]['MDE'])
        tstr[aset] += ' \\\\ \n'

        tstr[aset] += ' & Mean Abs.'
        for aref in trefs:
            tstr[aset] += ' & {:.2f} & {:.2f}'\
                .format(wd[aref]['MAFE'], wd[aref]['MADE'])
        tstr[aset] += ' \\\\ \n'

        tstr[aset] += ' & Root-mean-squared'
        for aref in trefs:
            tstr[aset] += ' & {:.2f} & {:.2f}'\
                .format(wd[aref]['VARFE']**(0.5), wd[aref]['VARDE']**(0.5))

        tstr[aset] += ' \\\\ \\hline \n'

        capstr = 'Functional (FE) and density error (DE) metrics for {:} on the {:} set, all in kcal/mol.\n'.format(str_rep_rules(dfa),aset)
        capstr += 'Each set of two columns uses a different proxy for the exact functional as in Table V of the main text.\n'
        capstr += 'Error statistics are reported in the bottom three rows.\n'
        capstr += 'The reaction index used by the GMTKN and the symbolic index are given in the first two columns.'

        sstr = '|rr'*len(trefs.keys())
        tex_str += '\\begin{longtable}{ll'+sstr+'}\n'
        tex_str += '  \\caption{\n'+capstr+'}\n'
        tex_str += '  \\label{tab:'+'dc_dft_{:}_{:}'.format(aset,dfa)+'}\n'
        tex_str += hstr[aset]
        tex_str += '  \\endfirsthead\n'
        tex_str += hstr[aset]
        tex_str += '  \\endhead\n'
        tex_str += tstr[aset]
        tex_str += '\\end{longtable}\n\n'

        """
        with open(out_dir + '/{:}_{:}_DC_DFT.tex'.format(dfa,aset),'w+') as tfl:
            tfl.write('\\begin{longtable}{ll'+sstr+'}\n')
            tfl.write('  \\caption{\n'+capstr+'}\n')
            tfl.write('  \\label{tab:'+'dc_dft_{:}_{:}'.format(aset,dfa)+'}\n')
            tfl.write(hstr[aset])
            tfl.write('  \\endfirsthead\n')
            tfl.write(hstr[aset])
            tfl.write('  \\endhead\n')
            #tfl.write('  \\centering\n')
            tfl.write(tstr[aset])
            tfl.write('\\end{longtable}')
        """

    """
    if idfa > 0:
        tex_str += '\n\\clearpage\n'
    tex_str += '\\subsection{'+str_rep_rules(dfa)+'}\n\n'
    tex_str += '\\input{./supp_mater_tabs/'+'/{:}_BH76_DC_DFT.tex'.format(dfa)+'}\n'
    tex_str += '\\input{./supp_mater_tabs/'+'/{:}_BH76RC_DC_DFT.tex'.format(dfa)+'}\n'
    """

with open(out_dir + '/dc_dft_for_sm.tex','w+') as tfl:
    tfl.write(tex_str)
