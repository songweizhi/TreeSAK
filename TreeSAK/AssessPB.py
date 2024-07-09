import os
import argparse


AssessPB_usage = '''
====================== AssessPB example commands ======================

# Dependency: bpcomp and tracecomp (from PhyloBayes-MPI)

export OMPI_MCA_btl=^openib

TreeSAK AssessPB -c1 c1dir/c1 -c2 c2dir/c2
TreeSAK AssessPB -c1 c1dir/c1 -c2 c2dir/c2 -c3 c3dir/c3
TreeSAK AssessPB -c1 c1dir/c1 -c2 c2dir/c2 -c3 c3dir/c3 -c4 c4dir/c4
TreeSAK AssessPB -cdir chain_dir

# This is a wrapper for: 
bpcomp -x 1000 10 c1 c2
bpcomp -x 1000 10 c1 c2 c3 c4
tracecomp -x 1000 c1 c2
tracecomp -x 1000 c1 c2 c3 c4

=======================================================================
'''


def compare2chains(chain_1, chain_2, chain_3, chain_4, burn_in, sample_interval, with_bpcomp, with_tracecomp, op_dir, cmd_txt):

    # bpcomp:    -x <burnin> [<every> <until>]. default burnin = 10 percent of the chain
    # tracecomp: -x <burnin> [<every> <until>]. default burnin = 20 percent of the chain

    bpcomp_cmd    = 'bpcomp -o %s/bpcomp -x %s %s %s %s'    % (op_dir, burn_in, sample_interval, chain_1, chain_2)
    tracecomp_cmd = 'tracecomp -o %s/tracecomp -x %s %s %s' % (op_dir, burn_in, chain_1, chain_2)

    if (chain_3 is not None) and (chain_4 is None):
        bpcomp_cmd    = 'bpcomp -o %s/bpcomp -x %s %s %s %s %s'    % (op_dir, burn_in, sample_interval, chain_1, chain_2, chain_3)
        tracecomp_cmd = 'tracecomp -o %s/tracecomp -x %s %s %s %s' % (op_dir, burn_in, chain_1, chain_2, chain_3)

    if (chain_3 is not None) and (chain_4 is not None):
        bpcomp_cmd    = 'bpcomp -o %s/bpcomp -x %s %s %s %s %s %s'    % (op_dir, burn_in, sample_interval, chain_1, chain_2, chain_3, chain_4)
        tracecomp_cmd = 'tracecomp -o %s/tracecomp -x %s %s %s %s %s' % (op_dir, burn_in, chain_1, chain_2, chain_3, chain_4)

    # write out commands
    cmd_txt_handle = open(cmd_txt, 'a')
    cmd_txt_handle.write(bpcomp_cmd + '\n')
    cmd_txt_handle.write(tracecomp_cmd + '\n')
    cmd_txt_handle.close()

    # execute commands
    if with_bpcomp is True:
        print()
        print('====================== bpcomp ======================')
        print()
        print(bpcomp_cmd)
        os.system(bpcomp_cmd)
        print('Guideline')
        print('1. maxdiﬀ < 0.1: good run.')
        print('2. maxdiﬀ < 0.3: acceptable: gives a good qualitative picture of the posterior consensus.')
        print('3. 0.3 < maxdiﬀ < 1: the sample is not yet suﬃciently large and have not converged, but on right track.')
        print('4. if maxdiﬀ = 1 even after 10,000 points: at least one run stuck in a local maximum.')
        print()

    if with_tracecomp is True:
        print('==================== tracecomp ====================')
        print()
        print(tracecomp_cmd)
        print()
        os.system(tracecomp_cmd)
        print()
        print('Guideline')
        print('1. rel diﬀ < 0.1 and minimum eﬀective size > 300: good run.')
        print('2. rel diﬀ < 0.3 and minimum eﬀective size > 50: acceptable run.')
        print()

    print('====================================================')


def AssessPB(args):

    chain_1         = args['c1']
    chain_2         = args['c2']
    chain_3         = args['c3']
    chain_4         = args['c4']
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
        compare2chains(chain_1, chain_2, chain_3, chain_4, burn_in, sample_interval, with_bpcomp, with_tracecomp, op_dir, cmd_txt)

    elif (chain_1 is None) and (chain_2 is None) and (chain_dir is not None):
        print('Compare multiple chains')
        print('Function to be added!')
        print('Program exited!')

    else:
        print('Please compare either no more than four chains (specified by -c1, -c2, -c3 and -c4) or multiple chains provided within -cdir')
        print('Program exited!')
        exit()


if __name__ == '__main__':

    AssessPB_parser = argparse.ArgumentParser()
    AssessPB_parser.add_argument('-c1',     required=False, default=None,           help='chain 1')
    AssessPB_parser.add_argument('-c2',     required=False, default=None,           help='chain 2')
    AssessPB_parser.add_argument('-c3',     required=False, default=None,           help='chain 3')
    AssessPB_parser.add_argument('-c4',     required=False, default=None,           help='chain 4')
    AssessPB_parser.add_argument('-cdir',   required=False, default=None,           help='chain folder')
    AssessPB_parser.add_argument('-bi',     required=False, default=1000,           help='burn-in, default: 1000')
    AssessPB_parser.add_argument('-si',     required=False, default=10,             help='sample interval, default: 10')
    AssessPB_parser.add_argument('-o',      required=True, default=None,            help='output directory')
    AssessPB_parser.add_argument('-f',      required=False, action="store_true",    help='force overwrite')
    args = vars(AssessPB_parser.parse_args())
    AssessPB(args)
