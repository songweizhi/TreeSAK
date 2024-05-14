import os
import argparse
from Bio import AlignIO


PB_usage = '''
========================================== PB example commands ==========================================

# Dependency: mpirun, pb_mpi and readpb_mpi (from PhyloBayes-MPI)

export OMPI_MCA_btl=^openib
TreeSAK PB -i in.phylip -p chain_name -t 12

# Notes:
1. This is a wrapper for: mpirun -np 12 pb_mpi -d in.phylip -cat -gtr -x 10 -1 -dgam 4 -s chain_name

2. Input MSA need to be in phylip format.

3. Settings used by Nina Dombrowski: -cat -gtr -x 10 -1 -dgam 4
   For each marker protein family, four parallel chains were run until convergence was reached, unless stated 
   otherwise (maxdiff < 0.3; settings: bpcomp -x 25%_burnin chain1 chain2 chain3 chain4). Additionally, we 
   checked for the minimum effective size using tracecomp (minimum effective size > 50; settings: -x 25%_burnin 
   chain1 chain2 chain3 chain4).

4. Settings used by Fan Lu:
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

    msa_in              = args['i']
    op_prefix           = args['p']
    fa_to_plp           = args['fa2plp']
    num_of_threads      = args['t']

    settings_dombrowski = '-cat -gtr -x 10 -1 -dgam 4'
    settings_fan        = ''

    msa_to_use = msa_in
    if fa_to_plp is True:
        msa_in_plp = '%s.phylip' % msa_in
        fa2phy(msa_in, msa_in_plp)
        msa_to_use = msa_in_plp

    pb_mpi_cmd          = 'mpirun -np %s pb_mpi -d %s %s -s %s' % (num_of_threads, msa_to_use, settings_dombrowski, op_prefix)
    print(pb_mpi_cmd)
    os.system(pb_mpi_cmd)

    readpb_cmd          = ''
    print(readpb_cmd)

    bpcomp_cmd          = ''
    print(bpcomp_cmd)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-i',       required=True,                          help='input MSA file')
    parser.add_argument('-p',       required=True,                          help='output prefix')
    parser.add_argument('-fa2plp',  required=False, action="store_true",    help='convert MSA format from fasta to phylip')
    parser.add_argument('-t',       required=False, type=int, default=2,    help='num of threads')
    args = vars(parser.parse_args())
    PB(args)


'''

export OMPI_MCA_btl=^openib
cd /scratch/PI/ocessongwz/Sponge_2023_12_01/MarkerRef2Tree_marker_set_4/s08_iqtree_wd/pb_test
mpirun -np 12 pb_mpi -d marker_pa85_pruner0.BMGE.phylip -cat -gtr -x 10 -1 -dgam 4 -s chain_test

'''

