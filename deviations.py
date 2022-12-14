import yaml
from os import path,system
from analysis import BH76_analysis
from SP_CSV_to_YAML import conversion
from index_labels import indx_to_name

#yaml.Dumper.ignore_aliases = lambda *args : True

tex_str = ''

sets = {}
sets['BH76'], sets['BH76RC'] = indx_to_name()

key_to_set = {'BH76': 'RX', 'BH76RC': 'RC'}

conversion() # ensures FLOSIC errors are recomputed
no_flo_dat = ['LCwPBE','BLYP','B3LYP', 'M06-L','MN15-L']

rdir = path.dirname(path.realpath(__file__)) + '/'

out_dir = './deviations/'
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

ref_ediff = {'BH76': {}, 'BH76RC': {}}
ref_ediff['BH76'] = yaml.load(open(rdir + 'BH76_ref_energies.yaml','r'),\
    Loader=yaml.Loader)
ref_ediff['BH76RC'] = yaml.load(open(rdir + 'BH76RC_ref_energies.yaml','r'),\
    Loader=yaml.Loader)


refd = { 'BH76': {}, 'BH76RC': {} }

dfas = ['LSDA','PBE','BLYP','SCAN','r2SCAN','M06-L','MN15-L','B3LYP','LCwPBE','SCAN@HF']

metrics = ['MD','MAD','RMSD','MAPE']
ncol = len(refs.keys())+1
summ_str = {'BH76': 'BH76', 'BH76RC': 'BH76RC'}
for aset in summ_str:
    for metric in metrics:
        if metric == 'MAPE':
            continue
        summ_str[aset] += ' & \\multicolumn{'+'{:}'.format(ncol)+'}{c}{'+'{:}}}'.format(metric)
    summ_str[aset] += ' \\\\ \n'
    for metric in metrics:
        if metric == 'MAPE':
            continue
        summ_str[aset] += ' & Exact'
        for aref in refs:
            summ_str[aset] += ' & {:}'.format(aref)
    summ_str[aset] += ' \\\\ \\hline \n'

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
    hstr = {'BH76': '\\\\ \\hline\nIndex & Reaction & Barrier Height & Exact', 'BH76RC': '\\\\ \\hline\nIndex & Reaction & Reaction Energy & Exact'}
    trefs = {}

    wrefs = {}

    tdfa = dfa.split('@')[0]

    for aset in hstr:

        for aref in refs:
            if dfa == aref or ((aref == 'SCAN-FLOSIC') and (dfa in no_flo_dat)):
                continue

            hstr[aset] += ' & {:}'.format(str_rep_rules(aref))
            trefs[aref] = refs[aref]

        hstr[aset] += ' \\\\ \\hline \n'

    wrefs['SCF'] = BH76_analysis(cdir=scf_dir+'{:}/{:}_BH76/'.format(tdfa,dfa))
    for aref in trefs:

        if aref == 'SCAN-FLOSIC':

            if dfa == 'SCAN@HF':
                tdir_scf = nrlmol_basis_dir+'SCAN/SCAN@HF_BH76/'
            else:
                tdir_scf = './FLOSIC/{:}/{:}_BH76/'.format(dfa,dfa)

            wrefs[aref] = yaml.load(open(tdir_scf + 'BH76_total.yaml','r'), Loader=yaml.Loader)

        else:
            wrefs[aref] = wrefs['SCF']

    for iset, aset in enumerate(tstr):

        tkey = key_to_set[aset]
        nrx = len(sets[aset].keys())

        wd = {}
        for aref in trefs:
            wd[aref] = {'MD': 0., 'MAD': 0., 'VAR': 0., 'MAPE': 0.}

        for irx, arx in enumerate(sets[aset]):

            tstr[aset] += '{:} & {:} & {:.2f} & {:.2f}'.format(arx,sets[aset][arx],\
                ref_ediff[aset][arx]['Ref'],wrefs[aref][tkey][arx]['Error'])

            for aref in trefs:

                terr = wrefs[aref][tkey][arx]['Energy'] - \
                    refd[aset][aref][arx]
                tstr[aset] += ' & {:.2f}'.format(terr)

                wd[aref]['MD'] += terr
                wd[aref]['MAD'] += abs(terr)
                wd[aref]['VAR'] += terr**2
                wd[aref]['MAPE'] += 100*abs(terr/refd[aset][aref][arx])

            tlchar = ''
            if irx == len(sets[aset].keys())-1:
                tlchar = '\\hline'

            tstr[aset] += ' \\\\ {:} \n'.format(tlchar)

        for aref in trefs:
            for akey in wd[aref]:
                wd[aref][akey] /= (1.*nrx)
            wd[aref]['RMSD'] = wd[aref]['VAR']**(0.5)

        if aset == 'BH76':
            nskey = 'Stats'
        elif aset == 'BH76RC':
            nskey = 'RC Stats'

        summ_str[aset] += '{:}'.format(dfa)
        for metric in metrics:
            if metric == 'MAPE':
                continue
            summ_str[aset] += ' & {:.2f}'.format(wrefs['SCF'][nskey][metric])
            for aref in refs:
                if aref in trefs:
                    summ_str[aset] += ' & {:.2f}'.format(wd[aref][metric])
                else:
                    summ_str[aset] += ' & '
        lchar = ''
        if idfa == len(dfas)-1 and iset < len(tstr.keys())-1:
            lchar = '\\hline'
        summ_str[aset] += ' \\\\ {:} \n'.format(lchar)

        for imetric, ametric in enumerate(metrics):

            lchar = ''
            if imetric == len(metrics)-1:
                lchar = '\\hline'

            tstr[aset] += ' & & {:} & {:.2f}'.format(ametric,\
                wrefs['SCF'][nskey][ametric])
            for aref in trefs:
                tstr[aset] += ' & {:.2f}'.format(wd[aref][ametric])
            tstr[aset] += ' \\\\ {:} \n'.format(lchar)

        capstr = 'Deviations (in kcal/mol) from the GMTKN reference values (Exact) and various functionals for the {:} set using {:}.\n'.format(aset,str_rep_rules(dfa))
        capstr += "The deviations are approximate errors; the error statistics in the ``Exact'' column correspond to the values in Table \\ref{tab:ak_pyscf_aug-cc-pvqz} of the main text.\n"
        capstr += 'The other columns replace the exact energies with those of a proxy.\n'
        capstr += 'For an approximate functional to be a good proxy for the exact one (in evaluating density- and functional-driven errors), the deviations computed with the high-level reference values and with the approximate functional should be comparable.\n'
        capstr += 'Mean deviations (MDs), mean absolute deviations (MADs), and variances (VARs) are reported in the bottom three rows.\n'
        capstr += 'The reaction index used by the GMTKN and the symbolic index are given in the first two columns.\n'

        if aset == 'BH76':
            tchar = 'Barrier Height'
        elif aset == 'BH76RC':
            tchar = 'Reaction Energy'

        capstr += "The sum of ``{:}'' and ``Exact'' columns yields the corresponding {:} {:}."\
            .format(tchar,dfa,tchar.lower())

        sstr = 'r'*(len(trefs.keys())+2)
        tex_str += '\\begin{longtable}{ll'+sstr+'}\n'
        tex_str += '  \\caption{\n'+capstr+'}\n'
        tex_str += '  \\label{tab:'+'devs_{:}_{:}'.format(aset,dfa)+'}\n'
        tex_str += hstr[aset]
        tex_str += '  \\endfirsthead\n'
        tex_str += hstr[aset]
        tex_str += '  \\endhead\n'
        tex_str += tstr[aset]
        tex_str += '\\end{longtable}\n\n'
        """
        with open(out_dir + '/{:}_{:}_devs.tex'.format(dfa,aset),'w+') as tfl:
            tfl.write('\\begin{longtable}{ll'+sstr+'}\n')
            tfl.write('  \\caption{\n'+capstr+'}\n')
            tfl.write('  \\label{tab:'+'devs_{:}_{:}'.format(aset,dfa)+'}\n')
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
    tex_str += '\\input{./supp_mater_tabs/'+'/{:}_BH76_devs.tex'.format(dfa)+'}\n'
    tex_str += '\\input{./supp_mater_tabs/'+'/{:}_BH76RC_devs.tex'.format(dfa)+'}\n'
    """

with open(out_dir + '/devs_for_sm.tex','w+') as tfl:
    tfl.write(tex_str)

with open(out_dir + '/dev_summ.tex','w+') as tfl:
    for aset in summ_str:
        tfl.write(str_rep_rules(summ_str[aset]))
