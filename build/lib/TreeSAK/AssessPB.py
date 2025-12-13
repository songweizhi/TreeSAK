import os
import argparse


AssessPB_usage = '''
====================== AssessPB example commands ======================

# Dependency: bpcomp and tracecomp (from PhyloBayes-MPI)
export OMPI_MCA_btl=^openib
TreeSAK AssessPB -c all_chains.txt

# This is a wrapper for (take 4 chains as an example): 
bpcomp -x 1000 10 c1 c2 c3 c4
tracecomp -x 1000 c1 c2 c3 c4

# format of the file provided to -c: directory_path/output_prefix
GTDB_SCG_best50p0_pb_chain1/GTDB_SCG_best50p0_pb_chain1 
GTDB_SCG_best50p0_pb_chain2/GTDB_SCG_best50p0_pb_chain2 
GTDB_SCG_best50p0_pb_chain3/GTDB_SCG_best50p0_pb_chain3 
GTDB_SCG_best50p0_pb_chain4/GTDB_SCG_best50p0_pb_chain4

=======================================================================
'''


def compare2chains(chain_1, chain_2, chain_3, chain_4, burn_in, sample_interval, op_dir, cmd_txt):

    # bpcomp:    -x <burnin> [<every> <until>]. default burnin = 10 percent of the chain
    # tracecomp: -x <burnin> [<every> <until>]. default burnin = 20 percent of the chain

    bpcomp_cmd    = 'bpcomp -o %s/bpcomp -x %s %s %s %s'                % (op_dir, burn_in, sample_interval, chain_1, chain_2)
    tracecomp_cmd = 'tracecomp -o %s/tracecomp -x %s %s %s'             % (op_dir, burn_in, chain_1, chain_2)

    if (chain_3 is not None) and (chain_4 is None):
        bpcomp_cmd    = 'bpcomp -o %s/bpcomp -x %s %s %s %s %s'         % (op_dir, burn_in, sample_interval, chain_1, chain_2, chain_3)
        tracecomp_cmd = 'tracecomp -o %s/tracecomp -x %s %s %s %s'      % (op_dir, burn_in, chain_1, chain_2, chain_3)

    if (chain_3 is not None) and (chain_4 is not None):
        bpcomp_cmd    = 'bpcomp -o %s/bpcomp -x %s %s %s %s %s %s'      % (op_dir, burn_in, sample_interval, chain_1, chain_2, chain_3, chain_4)
        tracecomp_cmd = 'tracecomp -o %s/tracecomp -x %s %s %s %s %s'   % (op_dir, burn_in, chain_1, chain_2, chain_3, chain_4)

    # write out commands
    cmd_txt_handle = open(cmd_txt, 'a')
    cmd_txt_handle.write(bpcomp_cmd + '\n')
    cmd_txt_handle.write(tracecomp_cmd + '\n')
    cmd_txt_handle.close()

    # execute commands
    print('\n================================ bpcomp ================================')
    os.system(bpcomp_cmd)
    print('\nGuideline')
    print('maxdiﬀ < 0.1:     good')
    print('maxdiﬀ < 0.3:     acceptable, gives a good qualitative picture of the posterior consensus.')
    print('0.3 < maxdiﬀ < 1: the sample is not yet suﬃciently large and have not converged, but on right track.')
    print('If maxdiﬀ = 1 even after 10,000 points: at least one run stuck in a local maximum.')
    print('\n============================== tracecomp ==============================\n')
    os.system(tracecomp_cmd)
    print('\nGuideline')
    print('good:       rel diﬀ < 0.1 and minimum eﬀective size > 300')
    print('acceptable: rel diﬀ < 0.3 and minimum eﬀective size > 50')
    print('\n========================================================================\n')


def AssessPB(args):

    chain_file      = args['c']
    burn_in         = args['bi']
    sample_interval = args['si']
    op_dir          = args['o']
    force_overwrite = args['f']

    cmd_txt         = '%s/cmds.txt' % op_dir

    # check is chain_file exist
    if os.path.isfile(chain_file) is False:
        print('%s not found, program exited!' % chain_file)
        exit()

    # check if chains were provided in the file
    chain_list = []
    for each_chain in open(chain_file):
        chain_list.append(each_chain.strip())
    if len(chain_list) < 2:
        print('Provided %s chains, need at least two chains, program exited!' % len(chain_list))
        exit()

    # create output dir
    if os.path.isdir(op_dir) is True:
        if force_overwrite is True:
            os.system('rm -r %s' % op_dir)
        else:
            print('output folder already exist, program exited!')
            exit()
    os.system('mkdir %s' % op_dir)

    if len(chain_list) == 2:
        compare2chains(chain_list[0], chain_list[1], burn_in, sample_interval, op_dir, cmd_txt)
    elif len(chain_list) == 3:
        compare2chains(chain_list[0], chain_list[1], chain_list[2], burn_in, sample_interval, op_dir, cmd_txt)
    elif len(chain_list) == 4:
        compare2chains(chain_list[0], chain_list[1], chain_list[2], chain_list[3], burn_in, sample_interval, op_dir, cmd_txt)


if __name__ == '__main__':

    AssessPB_parser = argparse.ArgumentParser()
    AssessPB_parser.add_argument('-c',      required=False, default=None,           help='a txt file contain all the chains')
    AssessPB_parser.add_argument('-bi',     required=False, default=1000,           help='burn-in, default: 1000')
    AssessPB_parser.add_argument('-si',     required=False, default=10,             help='sample interval, default: 10')
    AssessPB_parser.add_argument('-o',      required=True,  default=None,           help='output directory')
    AssessPB_parser.add_argument('-f',      required=False, action="store_true",    help='force overwrite')
    args = vars(AssessPB_parser.parse_args())
    AssessPB(args)
