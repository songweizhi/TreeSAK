import os
import argparse


pb_mpi_usage = '''
======================= pb_mpi example commands =======================

# Dependency: PhyloBayes-MPI

export OMPI_MCA_btl=^openib
TreeSAK pb_mpi -i in.phylip -p chain_name -t 12

# It is a wrapper for:
mpirun -np 2 pb_mpi -d in.phylip -cat -gtr -x 10 -1 -dgam 4 -s chain_name

# Input MSA need to be in phylip format.

# Settings used by Nina Dombrowski: -cat -gtr -x 10 -1 -dgam 4
For each marker protein family, four parallel chains were run until convergence was reached, 
unless stated otherwise (maxdiff < 0.3; settings: bpcomp -x 25%_burnin chain1 chain2 chain3 chain4). 
Additionally, we checked for the minimum effective size using tracecomp 
(minimum effective size > 50; settings: -x 25%_burnin chain1 chain2 chain3 chain4).

=======================================================================
'''


def pb_mpi(args):

    msa_in          = args['i']
    op_prefix       = args['p']
    num_of_threads  = args['t']

    dombrowski_settings = '-cat -gtr -x 10 -1 -dgam 4'

    pb_mpi_cmd = 'mpirun -np %s pb_mpi -d %s %s -s %s' % (num_of_threads, msa_in, dombrowski_settings, op_prefix)
    print(pb_mpi_cmd)
    os.system(pb_mpi_cmd)


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-i',   required=True,                          help='input MSA file')
    parser.add_argument('-p',   required=True,                          help='output prefix')
    parser.add_argument('-t',   required=False, type=int, default=2,    help='num of threads')
    args = vars(parser.parse_args())
    pb_mpi(args)


'''

export OMPI_MCA_btl=^openib
cd /scratch/PI/ocessongwz/Sponge_2023_12_01/MarkerRef2Tree_marker_set_4/s08_iqtree_wd/pb_test
mpirun -np 12 pb_mpi -d marker_pa85_pruner0.BMGE.phylip -cat -gtr -x 10 -1 -dgam 4 -s chain_test

'''

