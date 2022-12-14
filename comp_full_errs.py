import yaml
from os import path,system

from analysis import BH76_analysis
from SP_CSV_to_YAML import conversion
from index_labels import indx_to_name

rdir = path.dirname(path.realpath(__file__)) + '/'

conversion() # ensures FLOSIC errors are recomputed
sets = {}
sets['BH76'], sets['BH76RC'] = indx_to_name()

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


refd = {}
for aset in ['BH76','BH76RC']:
    refd[aset] = {}
    td = yaml.load(open(rdir + aset + '_ref_energies.yaml','r'), \
        Loader=yaml.Loader)
    for arx in td:
        refd[aset][arx] = td[arx]['Ref']

stat_l = ['MD','MAD','RMSD','MAPE']

scf_dir =  './results_aug-cc-pvqz/'
flosic_dir = './FLOSIC/'
hyb_dir = './scan_hybs/'

mad_div_str = 'DFA & MAD (first 12) & MAD (last 64) \\\\ \hline\n'

dfas = {
    'LSDA': ['SCF','HF','LCwPBE','S50X','SCAN-FLOSIC','LSDA-FLOSIC'],
    'PBE': ['SCF','HF','LCwPBE','S50X','SCAN-FLOSIC','PBE-FLOSIC'],
    'BLYP': ['SCF','HF','LCwPBE','S50X'],
    'SCAN': ['SCF','HF','LCwPBE','S50X','SCAN-FLOSIC'],
    'r2SCAN': ['SCF','HF','LCwPBE','S50X','SCAN-FLOSIC','r2SCAN-FLOSIC'],
    'M06-L': ['SCF','HF','LCwPBE','S50X'],
    'MN15-L': ['SCF','HF','LCwPBE','S50X'],
    'LCwPBE': ['SCF','HF','S50X'],
    'B3LYP': ['SCF','HF','LCwPBE','S50X'],
    'S50X': ['SCF'],
    'LSDA-FLOSIC': ['SCF'],
    'PBE-FLOSIC': ['SCF'],
    'SCAN-FLOSIC': ['SCF'],
    'r2SCAN-FLOSIC': ['SCF']
}

odir = './full_tables/'
if not path.isdir(odir):
    system('mkdir -p ' + odir)

set_to_key = {'BH76': 'RX', 'BH76RC': 'RC'}

otexs = ''

for idfa, dfa in enumerate(dfas):

    if idfa > 0:
        otexs += '\n\\clearpage\n'
    otexs += '\\subsection{'+str_rep_rules(dfa)+'}\n\n'

    tex_str = {'BH76': '', 'BH76RC': ''}
    hstr = {'BH76': '', 'BH76RC': ''}
    for aset in tex_str:
        hstr[aset] += '\\\\ \\hline \nIndex & Reaction & Exact'

    wd = {}
    for aref in dfas[dfa]:
        lbl = '@{:}'.format(aref)

        if 'FLOSIC' in aref or 'FLOSIC' in dfa:
            if 'FLOSIC' in dfa:
                tfl = flosic_dir + '{:}/{:}_BH76/BH76_total.yaml'\
                    .format(dfa,dfa,aref)
            else:
                tfl = flosic_dir + '{:}/{:}@{:}_BH76/BH76_total.yaml'\
                    .format(dfa,dfa,aref)
            wd[aref] = yaml.load(open(tfl,'r'),Loader=yaml.Loader)
        else:
            if aref == 'SCF':
                if dfa[0] == 'S' and dfa[-1] == 'X':
                    frac = dfa[1:-1]
                    tdir = hyb_dir + '{:}/SCAN_{:}_EXX_BH76/'.format(dfa,frac)
                else:
                    tdir = scf_dir + '{:}/{:}_BH76/'.format(dfa,dfa)
                lbl = 'SCF'
            elif aref[0] == 'S' and aref[-1] == 'X':
                tdir = hyb_dir + '{:}/{:}@{:}_BH76/'.format(aref,dfa,aref)
            else:
                tdir = scf_dir + '{:}/{:}@{:}_BH76/'.format(dfa,dfa,aref)
            wd[aref] = BH76_analysis(cdir=tdir)

        for aset in tex_str:
            hstr[aset] += ' & {:}'.format(lbl)

        if dfa in ['LCwPBE','S50X','SCAN-FLOSIC'] and aref == 'SCF':
            mad_12 = 0.
            mad_rest = 0.
            for irx in range(1,13,1):
                mad_12 += abs(wd['SCF']['RX'][irx]['Error'])
            for irx in range(13,77,1):
                mad_rest += abs(wd['SCF']['RX'][irx]['Error'])
            #print(\
            #    '{:} MAD for first 12 reactions = {:.2f} kcal/mol; and {:.2f} kcal/mol for last 64 reactions'\
            #    .format(dfa,mad_12/12.,mad_rest/64.))
            mad_div_str += '{:} & {:.2f} & {:.2f} \\\\ \n'\
                .format(dfa,mad_12/12.,mad_rest/64.)

    for aset in tex_str:

        hstr[aset] += ' \\\\ \\hline \n'

        for irx, arx in enumerate(refd[aset]):
            tex_str[aset] += '{:} & {:} & {:.2f}'.format(arx,sets[aset][arx],\
                refd[aset][arx])

            for aref in dfas[dfa]:
                tex_str[aset] += ' & {:.2f}'.format(wd[aref][set_to_key[aset]][arx]['Error'])

            lchar = ''
            if irx == len(refd[aset])-1:
                lchar = '\\hline'
            tex_str[aset] += ' \\\\ {:} \n'.format(lchar)

        if aset == 'BH76':
            tkey = 'Stats'
        elif aset == 'BH76RC':
            tkey = 'RC Stats'

        for istat, astat in enumerate(stat_l):
            tex_str[aset] += ' & & {:}'.format(astat)
            for aref in dfas[dfa]:
                tex_str[aset] += ' & {:.2f}'.format(wd[aref][tkey][astat])

            lchar = ''
            if istat == len(stat_l)-1:
                lchar = '\\hline'
            tex_str[aset] += ' \\\\ {:} \n'.format(lchar)

        if aset == 'BH76':
            tchar = 'barrier height'
        elif aset == 'BH76RC':
            tchar = 'reaction energy'

        capstr = 'Individual {:} errors for {:} including its self-consistent (SCF) and non-selfconsistent values (indicated with an ``@'' symbol) used throughout.\n'.format(aset,dfa)
        capstr += 'All values are in kcal/mol.\n'
        capstr += "The ``Exact'' column lists the reference values, such that the sum of the ``Exact'' column and any other column would yield the {:} computed with {:} at that density.\n".format(tchar,dfa)
        capstr += 'The reaction index used by the GMTKN and the symbolic name are given in the first two columns.'

        sstr = 'r'*(len(dfas[dfa])+1)
        otexs += '\\begin{longtable}{ll'+sstr+'}\n'
        otexs += '  \\caption{\n'+str_rep_rules(capstr)+'}\n'
        otexs += '  \\label{tab:'+'{:}_errors_{:}'.format(aset,dfa)+'}\n'
        otexs += hstr[aset]
        otexs += '  \\endfirsthead\n'
        otexs += hstr[aset]
        otexs += '  \\endhead\n'
        otexs += tex_str[aset]
        otexs += '\\end{longtable}\n\n'

        """
        with open(odir + '/{:}_{:}_full_errors.tex'.format(dfa,aset),'w+') as tfl:
            tfl.write('\\begin{longtable}{ll'+sstr+'}\n')
            tfl.write('  \\caption{\n'+str_rep_rules(capstr)+'}\n')
            tfl.write('  \\label{tab:'+'{:}_errors_{:}'.format(aset,dfa)+'}\n')
            tfl.write(hstr[aset])
            tfl.write('  \\endfirsthead\n')
            tfl.write(hstr[aset])
            tfl.write('  \\endhead\n')
            #tfl.write('  \\centering\n')
            tfl.write(tex_str[aset])
            tfl.write('\\end{longtable}')
        """

    """
    if idfa > 0:
        otexs += '\n\\clearpage\n'
    otexs += '\\subsection{'+str_rep_rules(dfa)+'}\n\n'
    otexs += '\\input{./SM_error_tables/'+'/{:}_BH76_full_errors.tex'.format(dfa)+'}\n'
    otexs += '\\input{./SM_error_tables/'+'/{:}_BH76RC_full_errors.tex'.format(dfa)+'}\n'
    """

with open(odir + 'full_err_for_sm.tex','w+') as tfl:
    tfl.write(otexs)

with open(odir + 'first_12_mad.tex','w+') as tfl:
    tfl.write(str_rep_rules(mad_div_str))
