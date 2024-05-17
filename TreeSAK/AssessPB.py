import os
import argparse


AssessPB_usage = '''
================= AssessPB example commands =================

# Dependency: bpcomp and tracecomp (from PhyloBayes-MPI)

export OMPI_MCA_btl=^openib
TreeSAK AssessPB -c1 chain_1 -c2 chain_2
TreeSAK AssessPB -cdir chain_folder

# This is a wrapper for: 
bpcomp -x 1000 10 chain1 chain2
tracecomp -x 1000 chain1 chain2

=============================================================
'''


def compare2chains(chain_1, chain_2, burn_in, sample_interval, with_bpcomp, with_tracecomp, op_dir, cmd_txt):

    # bpcomp:    -x <burnin> [<every> <until>]. default burnin = 10 percent of the chain
    # tracecomp: -x <burnin> [<every> <until>]. default burnin = 20 percent of the chain

    bpcomp_cmd    = 'bpcomp -o %s/bpcomp -x %s %s %s %s'    % (op_dir, burn_in, sample_interval, chain_1, chain_2)
    tracecomp_cmd = 'tracecomp -o %s/tracecomp -x %s %s %s' % (op_dir, burn_in, chain_1, chain_2)

    cmd_txt_handle = open(cmd_txt, 'a')

    if with_bpcomp is True:
        cmd_txt_handle.write(bpcomp_cmd + '\n')
        os.system(bpcomp_cmd)

    if with_tracecomp is True:
        cmd_txt_handle.write(tracecomp_cmd + '\n')
        os.system(tracecomp_cmd)

    cmd_txt_handle.close()


def AssessPB(args):

    chain_1         = args['c1']
    chain_2         = args['c2']
    chain_dir       = args['cdir']
    burn_in         = args['bi']
    sample_interval = args['si']
    op_dir          = args['o']
    force_overwrite = args['f']
    with_bpcomp     = True
    with_tracecomp  = True

    cmd_txt         = '%s/cmds.txt' % op_dir

    # create output dir
    if os.path.isdir(op_dir) is True:
        if force_overwrite is True:
            os.system('rm -r %s' % op_dir)
        else:
            print('output folder already exist, program exited!')
            exit()
    os.system('mkdir %s' % op_dir)

    if (chain_1 is not None) and (chain_2 is not None) and (chain_dir is None):
        print('Compare two chains')
        compare2chains(chain_1, chain_2, burn_in, sample_interval, with_bpcomp, with_tracecomp, op_dir, cmd_txt)

    elif (chain_1 is None) and (chain_2 is None) and (chain_dir is not None):
        print('Compare multiple chains')
        print('Function to be added!')
        print('Program exited!')

    else:
        print('Please compare either two chains (specified by -c1 and -c2) or multiple chains provided within -cdir')
        print('Program exited!')
        exit()


if __name__ == '__main__':

    AssessPB_parser = argparse.ArgumentParser()
    AssessPB_parser.add_argument('-c1',     required=False, default=None,           help='chain 1')
    AssessPB_parser.add_argument('-c2',     required=False, default=None,           help='chain 2')
    AssessPB_parser.add_argument('-cdir',   required=False, default=None,           help='chain folder')
    AssessPB_parser.add_argument('-bi',     required=False, default=1000,           help='burn-in, default: 1000')
    AssessPB_parser.add_argument('-si',     required=False, default=10,             help='sample interval, default: 10')
    AssessPB_parser.add_argument('-o',      required=True, default=None,            help='output directory')
    AssessPB_parser.add_argument('-f',      required=False, action="store_true",    help='force overwrite')
    args = vars(AssessPB_parser.parse_args())
    AssessPB(args)
