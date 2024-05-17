import os
import argparse
from Bio import AlignIO
import multiprocessing as mp


PB_usage = '''
========================================== PB example commands ==========================================

# Dependency: mpirun, pb_mpi and readpb_mpi (from PhyloBayes-MPI)

export OMPI_MCA_btl=^openib
TreeSAK PB -i in.phylip -p chain_name -t 12

# Notes:
1. This is a wrapper for: mpirun -np 12 pb_mpi -d in.phylip -cat -gtr -x 10 -1 -dgam 4 -s chain_name
2. Input MSA need to be in phylip format.
3. To stop a chain, just open the <chain_name>.run ﬁle and replace the 1 by a 0 (echo 0 > <chain_name>.run).
4. Be careful not to restart an already running chain.
5. You can stop a chain and restart it under a diﬀerent degree of parallelization.
6. Generally, PhyloBayes provides good results for a total number of points of 10000-30000.

*  Settings used by Nina Dombrowski: -cat -gtr -x 10 -1 -dgam 4
   For each marker protein family, four parallel chains were run until convergence was reached, unless stated 
   otherwise (maxdiff < 0.3; settings: bpcomp -x 25_burnin chain1 chain2 chain3 chain4). Additionally, we 
   checked for the minimum effective size using tracecomp (minimum effective size > 50; settings: -x 25_burnin 
   chain1 chain2 chain3 chain4).
   
*  Settings used by Fan Lu:
   Four chains were run for each consensus tree, and for each chain over 15,000 cycles (5,000 burn-in) 
   were conducted, until a maxdiff value lower than 0.3 was reached. Otherwise, non-converged chains were 
   continually run to over 20,000 cycles. Posterior predictive tests were conducted using PhyloBayes MPI 
   with the ‘readpb_mpi -x 5000 50 -allppred’ command.

=========================================================================================================
'''


def fa2phy(fasta_in, phy_out):

    alignment = AlignIO.read(fasta_in, 'fasta')

    max_seq_id_len = 0
    for each_seq in alignment:
        seq_id_len = len(each_seq.id)
        if seq_id_len > max_seq_id_len:
            max_seq_id_len = seq_id_len

    with open(phy_out, 'w') as msa_out_handle:
        msa_out_handle.write('%s %s\n' % (len(alignment), alignment.get_alignment_length()))
        for each_seq in alignment:
            seq_id = each_seq.id
            seq_id_with_space = '%s%s' % (seq_id, ' ' * (max_seq_id_len + 2 - len(seq_id)))
            msa_out_handle.write('%s%s\n' % (seq_id_with_space, str(each_seq.seq)))


def PB(args):

    msa_in          = args['i']
    op_dir          = args['o']
    op_prefix       = args['p']
    fa_to_plp       = args['fa2plp']
    num_of_threads  = args['t']
    num_of_chains   = args['chain']
    force_overwrite = args['f']

    ####################################################################################################################

    settings_dombrowski = '-cat -gtr -x 10 -1 -dgam 4'
    setting_to_use      = settings_dombrowski
    msa_in_plp          = '%s/%s.phylip'        % (op_dir, msa_in)
    cmd_txt             = '%s/cmds.txt'         % op_dir

    ####################################################################################################################

    # create output dir
    if os.path.isdir(op_dir) is True:
        if force_overwrite is True:
            os.system('rm -r %s' % op_dir)
        else:
            print('output folder already exist, program exited!')
            exit()
    os.system('mkdir %s' % op_dir)

    # fa_to_phylip
    msa_to_use = msa_in
    if fa_to_plp is True:
        fa2phy(msa_in, msa_in_plp)
        msa_to_use = msa_in_plp

    pb_mpi_cmd_list = []
    jobs_to_run_in_parallel = 0
    if num_of_chains == 1:
        jobs_to_run_in_parallel = 1
        pb_mpi_cmd = 'mpirun -np %s pb_mpi -d %s %s -s %s/%s' % (num_of_threads, msa_to_use, setting_to_use, op_dir, op_prefix)
        pb_mpi_cmd_list.append(pb_mpi_cmd)

    elif num_of_threads <= num_of_chains:
        jobs_to_run_in_parallel = num_of_threads
        for chain_index in range(1, (num_of_chains + 1)):
            current_wd = '%s/%s_chain%s' % (op_dir, op_prefix, chain_index)
            os.mkdir(current_wd)
            pb_mpi_cmd = 'mpirun -np %s pb_mpi -d %s %s -s %s/%s_chain%s' % (1, msa_to_use, setting_to_use, current_wd, op_prefix, chain_index)
            pb_mpi_cmd_list.append(pb_mpi_cmd)

    else:
        jobs_to_run_in_parallel = num_of_chains
        cores_per_run = num_of_threads // num_of_chains
        for chain_index in range(1, (num_of_chains + 1)):
            current_wd = '%s/%s_chain%s' % (op_dir, op_prefix, chain_index)
            os.mkdir(current_wd)
            pb_mpi_cmd = 'mpirun -np %s pb_mpi -d %s %s -s %s/%s_chain%s' % (cores_per_run, msa_to_use, setting_to_use, current_wd, op_prefix, chain_index)
            pb_mpi_cmd_list.append(pb_mpi_cmd)

    # write out commands
    cmd_txt_handle = open(cmd_txt, 'w')
    for cmd in pb_mpi_cmd_list:
        cmd_txt_handle.write(cmd + '\n')
    cmd_txt_handle.close()

    # run chains with mp
    print('Running pb_mpi with multiprocessing')
    pool = mp.Pool(processes=jobs_to_run_in_parallel)
    pool.map(os.system, pb_mpi_cmd_list)
    pool.close()
    pool.join()

    # readpb_cmd = ''
    # bpcomp_cmd = ''
    # print(readpb_cmd)
    # print(bpcomp_cmd)


if __name__ == '__main__':

    PB_parser = argparse.ArgumentParser()
    PB_parser.add_argument('-i',       required=True,                          help='input MSA file')
    PB_parser.add_argument('-o',       required=True,                          help='output directory')
    PB_parser.add_argument('-p',       required=True,                          help='output prefix')
    PB_parser.add_argument('-fa2plp',  required=False, action="store_true",    help='convert MSA format from fasta to phylip')
    PB_parser.add_argument('-chain',   required=False, type=int, default=4,    help='num of chains to run in parallel, default: 4')
    PB_parser.add_argument('-t',       required=False, type=int, default=4,    help='num of cores, default: 4')
    PB_parser.add_argument('-f',       required=False, action="store_true",    help='force overwrite')
    args = vars(PB_parser.parse_args())
    PB(args)
