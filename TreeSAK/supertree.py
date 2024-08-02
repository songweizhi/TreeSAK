import os
import glob
import argparse
from ete3 import Tree
from Bio import AlignIO
import multiprocessing as mp


supertree_usage = '''
====================== supertree example commands ======================

Dependencies: mafft, trimal, bmge and iqtree2

TreeSAK supertree -i best10 -x fa -o best10_astral_tree -bmge -t 12 -f

========================================================================
'''


def sep_path_basename_ext(file_in):

    f_path, f_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'
    f_base, f_ext = os.path.splitext(f_name)

    return f_name, f_path, f_base, f_ext[1:]


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


def PB(msa_in, op_dir, op_prefix, fa_to_plp, num_of_threads, num_of_chains, force_overwrite):

    ####################################################################################################################

    msa_in_name, msa_in_path, msa_in_base, msa_in_ext = sep_path_basename_ext(msa_in)

    settings_dombrowski = '-cat -gtr -x 10 -1 -dgam 4'
    setting_to_use  = settings_dombrowski
    msa_in_plp      = '%s/%s.phylip'        % (op_dir, msa_in_base)
    cmd_txt         = '%s/%s_cmds.txt'      % (op_dir, msa_in_base)

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

    cores_per_chain = 0
    chain_name_list = []
    pb_mpi_cmd_list = []
    jobs_to_run_in_parallel = 0
    if num_of_chains == 1:
        jobs_to_run_in_parallel = 1
        pb_mpi_cmd = 'mpirun -np %s pb_mpi -d %s %s -s %s/%s' % (num_of_threads, msa_to_use, setting_to_use, op_dir, op_prefix)
        chain_name_list.append('%s/%s' % (op_dir, op_prefix))
        pb_mpi_cmd_list.append(pb_mpi_cmd)
        cores_per_chain = num_of_threads

    elif num_of_threads <= num_of_chains:
        jobs_to_run_in_parallel = num_of_threads
        for chain_index in range(1, (num_of_chains + 1)):
            current_wd = '%s/%s_chain%s' % (op_dir, op_prefix, chain_index)
            os.mkdir(current_wd)
            pb_mpi_cmd = 'mpirun -np %s pb_mpi -d %s %s -s %s/%s_chain%s' % (1, msa_to_use, setting_to_use, current_wd, op_prefix, chain_index)
            chain_name_list.append('%s/%s_chain%s' % (current_wd, op_prefix, chain_index))
            pb_mpi_cmd_list.append(pb_mpi_cmd)
            cores_per_chain = 1
    else:
        jobs_to_run_in_parallel = num_of_chains
        cores_per_run = num_of_threads // num_of_chains
        for chain_index in range(1, (num_of_chains + 1)):
            current_wd = '%s/%s_chain%s' % (op_dir, op_prefix, chain_index)
            os.mkdir(current_wd)
            pb_mpi_cmd = 'mpirun -np %s pb_mpi -d %s %s -s %s/%s_chain%s' % (cores_per_run, msa_to_use, setting_to_use, current_wd, op_prefix, chain_index)
            chain_name_list.append('%s/%s_chain%s' % (current_wd, op_prefix, chain_index))
            pb_mpi_cmd_list.append(pb_mpi_cmd)
            cores_per_chain = cores_per_run

    # write out commands
    cmd_txt_handle = open(cmd_txt, 'w')
    for cmd in pb_mpi_cmd_list:
        cmd_txt_handle.write(cmd + '\n')

    cmd_txt_handle.write('\n# To restart a terminated run (e.g., due to walltime limitation)\n')
    for each_chain in chain_name_list:
        cmd_txt_handle.write('mpirun -np %s pb_mpi %s\n' % (cores_per_chain, each_chain))
    cmd_txt_handle.close()

    # run chains with mp
    print('Running pb_mpi with multiprocessing')
    pool = mp.Pool(processes=jobs_to_run_in_parallel)
    pool.map(os.system, pb_mpi_cmd_list)
    pool.close()
    pool.join()

    # assess the results
    if num_of_chains > 1:

        readpb_cmd = 'bpcomp -x 1000 10 %s' % (' '.join(chain_name_list))
        bpcomp_cmd = 'tracecomp -x 1000 %s' % (' '.join(chain_name_list))

        # write out commands
        cmd_txt_handle = open(cmd_txt, 'a')
        cmd_txt_handle.write(readpb_cmd + '\n')
        cmd_txt_handle.write(bpcomp_cmd + '\n')
        cmd_txt_handle.close()

        # report
        print('You may want to use the following commands to assess the results:')
        print(readpb_cmd)
        print(bpcomp_cmd)

    print('Done!')


def supertree(args):

    oma_op_fasta            = args['i']
    fasta_file_ext          = args['x']
    op_dir                  = args['o']
    trim_with_bmge          = args['bmge']
    trim_model              = args['bmge_m']
    entropy_score_cutoff    = args['bmge_esc']
    iqtree_model            = args['iqtree_m']
    force_overwrite         = args['f']
    num_of_threads          = args['t']
    infer_pb_tree           = args['pb']

    # specify path to BMGE.jar
    current_file_path = '/'.join(os.path.realpath(__file__).split('/')[:-1])
    pwd_bmge_jar = '%s/BMGE.jar' % current_file_path

    fa_file_re   = '%s/*.%s' % (oma_op_fasta, fasta_file_ext)
    fa_file_list = ['.'.join(os.path.basename(i).split('.')[:-1]) for i in glob.glob(fa_file_re)]

    if len(fa_file_list) == 0:
        print('No file found in %s, program exited!' % oma_op_fasta)
        exit()

    ################################################################################

    # define file name
    cmd_1_mafft_txt         = '%s/cmd_1_mafft.txt'          % op_dir
    cmd_2_trim_txt          = '%s/cmd_2_trim.txt'           % op_dir
    cmd_3_tree_txt          = '%s/cmd_3_tree.txt'           % op_dir
    cmd_4_astral_txt        = '%s/cmd_4_astral.txt'         % op_dir
    aln_dir                 = '%s/dir_1_msa'                % op_dir
    trimmed_aln_dir         = '%s/dir_2_trimmed_msa'        % op_dir
    tree_dir                = '%s/dir_3_tree'               % op_dir
    combined_gene_tree_file = '%s/combined_trees.txt'       % op_dir
    astral_mapping_txt      = '%s/name_mapping.txt'         % op_dir
    consensus_tree_txt      = '%s/consensus_tree.treefile'  % op_dir

    ################################################################################

    # create output folder
    if os.path.isdir(op_dir) is True:
        if force_overwrite is True:
            os.system('rm -r %s' % op_dir)
        else:
            print('%s exist, program exited!' % op_dir)
            exit()

    os.mkdir(op_dir)
    os.mkdir(aln_dir)
    os.mkdir(trimmed_aln_dir)
    os.mkdir(tree_dir)

    ################################################################################

    cmd_list_mafft = []
    cmd_list_trim  = []
    cmd_list_tree  = []
    cmd_1_mafft_txt_handle = open(cmd_1_mafft_txt, 'w')
    cmd_2_trim_txt_handle = open(cmd_2_trim_txt, 'w')
    cmd_3_tree_txt_handle = open(cmd_3_tree_txt, 'w')
    for each_og in sorted(fa_file_list):

        # define file name
        current_gene_tree_dir   = '%s/dir_3_tree/%s'    % (op_dir, each_og)
        og_fa                   = '%s/%s.%s'            % (oma_op_fasta, each_og, fasta_file_ext)
        og_aln                  = '%s/%s.aln'           % (aln_dir, each_og)
        og_aln_trimmed          = '%s/%s_trimal.aln'    % (trimmed_aln_dir, each_og)
        if trim_with_bmge is True:
            og_aln_trimmed      = '%s/%s_bmge.aln'      % (trimmed_aln_dir, each_og)

        os.system('mkdir %s' % current_gene_tree_dir)

        # prepare commands
        mafft_cmd          = 'mafft-einsi --thread %s --quiet %s > %s'                                             % (1, og_fa, og_aln)

        trim_cmd           = 'trimal -in %s -out %s -automated1'                                                   % (og_aln, og_aln_trimmed)
        if trim_with_bmge is True:
            trim_cmd       = 'java -jar %s -i %s -m %s -t AA -h %s -of %s'                                         % (pwd_bmge_jar, og_aln, trim_model, entropy_score_cutoff, og_aln_trimmed)

        infer_tree_cmd     = 'iqtree2 -s %s --seqtype AA -m %s -B 1000 --wbtl --bnni --prefix %s/%s -T %s --quiet' % (og_aln_trimmed, iqtree_model, current_gene_tree_dir, each_og, num_of_threads)
        if infer_pb_tree is True:
            infer_tree_cmd = 'TreeSAK PB -i %s -o %s -p %s -t %s -n %s -fa2plp'                                    % (og_aln_trimmed, current_gene_tree_dir, each_og, num_of_threads, 4)

        # add commands to list
        cmd_list_mafft.append(mafft_cmd)
        cmd_list_trim.append(trim_cmd)
        cmd_list_tree.append(infer_tree_cmd)

        # write out commands
        cmd_1_mafft_txt_handle.write(mafft_cmd + '\n')
        cmd_2_trim_txt_handle.write(trim_cmd + '\n')
        cmd_3_tree_txt_handle.write(infer_tree_cmd + '\n')

    cmd_1_mafft_txt_handle.close()
    cmd_2_trim_txt_handle.close()
    cmd_3_tree_txt_handle.close()

    # run mafft commands
    print('Running mafft with %s cores for %s OGs' % (num_of_threads, len(fa_file_list)))
    pool = mp.Pool(processes=num_of_threads)
    pool.map(os.system, cmd_list_mafft)
    pool.close()
    pool.join()

    # run trim commands
    print('Trimming with %s cores for %s OGs' % (num_of_threads, len(fa_file_list)))
    pool = mp.Pool(processes=num_of_threads)
    pool.map(os.system, cmd_list_trim)
    pool.close()
    pool.join()

    # run iqtree commands
    if infer_pb_tree is False:
        print('Running iqtree with %s cores for %s OGs' % (num_of_threads, len(fa_file_list)))
        for each_iqtree_cmd in sorted(cmd_list_tree):
            print(each_iqtree_cmd)
            os.system(each_iqtree_cmd)
    else:
        print('Commands for inferring PhyloBayes tree exported to %s' % cmd_3_tree_txt)

    #################################################### run astral ####################################################

    if infer_pb_tree is False:

        # cat gene trees
        os.system('cat %s/*.treefile > %s' % (tree_dir, combined_gene_tree_file))

        gnm_to_gene_dict = dict()
        for each_tree in open(combined_gene_tree_file):
            tree_str = each_tree.strip()
            current_tree = Tree(tree_str, quoted_node_names=True, format=1)
            for node in current_tree.traverse():
                if node.is_leaf():
                    leaf_name = node.name
                    leaf_gnm = '_'.join(leaf_name.split('_')[:-1])
                    if leaf_gnm not in gnm_to_gene_dict:
                        gnm_to_gene_dict[leaf_gnm] = {leaf_name}
                    else:
                        gnm_to_gene_dict[leaf_gnm].add(leaf_name)

        # get the mapping file
        astral_mapping_txt_handle = open(astral_mapping_txt, 'w')
        for each_gnm in gnm_to_gene_dict:
            current_gene_set = gnm_to_gene_dict[each_gnm]
            for each_gene in current_gene_set:
                astral_mapping_txt_handle.write('%s\t%s\n' % (each_gene, each_gnm))
        astral_mapping_txt_handle.close()

        # may need to add more parameters
        astral_cmd = 'astral -i %s -o %s -t %s -a %s' % (combined_gene_tree_file, consensus_tree_txt, num_of_threads, astral_mapping_txt)
        # -r    --round     Integer 4       Number of initial rounds of placements
        # -s    --subsample Integer 4       Number of rounds of subsampling per exploration step

        # write out command
        cmd_4_astral_txt_handle = open(cmd_4_astral_txt, 'w')
        cmd_4_astral_txt_handle.write(astral_cmd + '\n')
        cmd_4_astral_txt_handle.close()

        # run astral
        os.system(astral_cmd)

    else:
        print('Things to do:')
        print('Run PhyloBayes with commands exported to %s' % cmd_3_tree_txt)
        print('Wait until your chains, for each of your protein family, reached convergence')
        print('1. get consensus gene tree for each protein family usingn')
        print('2. get species tree based on the consensus gene trees with trimal')

    ####################################################################################################################


if __name__ == '__main__':

    supertree_parser = argparse.ArgumentParser()
    supertree_parser.add_argument('-i',         required=True,                          help='orthologous gene sequence')
    supertree_parser.add_argument('-x',         required=True,                          help='faa file extension')
    supertree_parser.add_argument('-o',         required=True,                          help='output directory')
    supertree_parser.add_argument('-bmge',      required=False, action="store_true",    help='trim with BMGE, default is trimal')
    supertree_parser.add_argument('-bmge_m',    required=False, default='BLOSUM30',     help='trim model, default: BLOSUM30')
    supertree_parser.add_argument('-bmge_esc',  required=False, default='0.55',         help='entropy score cutoff, default: 0.55')
    supertree_parser.add_argument('-iqtree_m',  required=False, default='LG+G+I',       help='iqtree_model, default: LG+G+I')
    supertree_parser.add_argument('-pb',        required=False, action="store_true",    help='infer tree with PhyloBayes-MPI, default is iqtree')
    supertree_parser.add_argument('-f',         required=False, action="store_true",    help='force overwrite')
    supertree_parser.add_argument('-t',         required=False, type=int, default=1,    help='num of threads, default: 1')
    args = vars(supertree_parser.parse_args())
    supertree(args)
