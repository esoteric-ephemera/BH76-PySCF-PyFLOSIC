
from setup import setup
from os import system, path

exx_fracs = list(range(0,101,10))
exx_fracs += [25, 75]

sl_dfa_x = 'MGGA_X_R2SCAN'
sl_dfa_c = 'MGGA_C_R2SCAN'
sl_name = 'r2SCAN'

opts = {'basis': 'aug-cc-pvqz', 'write_chkfl':True}

bdir = './' + sl_name + '_hybrids/'
if not path.isdir(bdir):
    system('mkdir -p ' + bdir)

for ax in exx_fracs:
    if ax == 0:
        dfa_dir = bdir + sl_name
        xc_str = sl_dfa_x + ', ' + sl_dfa_c
    elif ax == 100:
        dfa_dir = bdir + sl_name + '_{:}_EXX'.format(ax)
        xc_str = 'HF, ' + sl_dfa_c
    else:
        dfa_dir = bdir + sl_name + '_{:}_EXX'.format(ax)
        xc_str = '{:.2f}*'.format(1.-ax/100.)+ sl_dfa_x \
            + ' + {:.2f}*HF'.format(ax/100.) + ', ' + sl_dfa_c

    setup(xc_str,startclean=True,xclib='LibXC', dfa_dir=dfa_dir, \
        inp_opts=opts)

sub_str = "#!/bin/bash\n\n"
sub_str += "for u in ./*/ ; do\n"
sub_str += "  cd $u\n"
sub_str += "  qsub runjob.sh\n"
sub_str += "  cd ..\n"
sub_str += "done"

with open(bdir + '/sub_all_jobs.sh','w+') as tfl:
    tfl.write(sub_str)
