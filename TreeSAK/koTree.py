import os
import glob
import argparse
from Bio import SeqIO
import multiprocessing as mp


koTree_usage = '''
================================ koTree example commands ================================

TreeSAK koTree -i combined.faa -kegg KEGG_wd -o op_dir -bmge -t 12 -f -fun ko_id.txt
TreeSAK koTree -i combined.faa -kegg KEGG_wd -o op_dir -bmge -t 12 -f -fun K01995
TreeSAK koTree -i combined.faa -kegg KEGG_wd -o op_dir -bmge -t 12 -f -fun K01995,K01996

=========================================================================================
'''


def select_seq(seq_file, seq_id_set, output_file):
    output_file_handle = open(output_file, 'w')
    for seq_record in SeqIO.parse(seq_file, 'fasta'):
        seq_id = seq_record.id
        if seq_id in seq_id_set:
            SeqIO.write(seq_record, output_file_handle, 'fasta-2line')
    output_file_handle.close()


def koTree(args):

    combined_faa            = args['i']
    kegg_annotation_wd      = args['kegg']
    interested_fun_txt      = args['fun']
    op_dir                  = args['o']
    trim_with_bmge          = args['bmge']
    trim_model              = args['bmge_m']
    entropy_score_cutoff    = args['bmge_esc']
    iqtree_model            = args['iqtree_m']
    force_overwrite         = args['f']
    num_of_threads          = args['t']

    # specify path to BMGE.jar
    current_file_path = '/'.join(os.path.realpath(__file__).split('/')[:-1])
    pwd_bmge_jar = '%s/BMGE.jar' % current_file_path

    interested_fun_set = set()
    if os.path.isfile(interested_fun_txt) is False:
        if ',' in interested_fun_txt:
            interested_fun_set = interested_fun_txt.split(',')
        else:
            interested_fun_set.add(interested_fun_txt)
    else:
        for each_fun in open(interested_fun_txt):
            interested_fun_set.add(each_fun.strip().split()[0])

    ################################################################################

    faa_dir                 = '%s/dir_1_faa'                % op_dir
    aln_dir                 = '%s/dir_2_msa'                % op_dir
    trimmed_aln_dir         = '%s/dir_3_trimmed_msa'        % op_dir
    tree_dir                = '%s/dir_4_tree'               % op_dir
    cmd_1_mafft_txt         = '%s/cmd_1_mafft.txt'          % op_dir
    cmd_2_trim_txt          = '%s/cmd_2_trim.txt'           % op_dir
    cmd_3_tree_txt          = '%s/cmd_3_tree.txt'           % op_dir

    ################################################################################

    # create output folder
    if os.path.isdir(op_dir) is True:
        if force_overwrite is True:
            os.system('rm -r %s' % op_dir)
        else:
            print('%s exist, program exited!' % op_dir)
            exit()

    os.mkdir(op_dir)
    os.mkdir(faa_dir)
    os.mkdir(aln_dir)
    os.mkdir(trimmed_aln_dir)
    os.mkdir(tree_dir)

    ################################################################################

    fun_to_gene_dict = dict()
    if kegg_annotation_wd is not None:

        print('Reading in KEGG annotation results')
        file_re = '%s/*KEGG_wd/*_ko_assignment_ABCD.txt' % (kegg_annotation_wd)
        file_list = glob.glob(file_re)

        if len(file_list) == 0:
            print('KEGG annotation file not detected, program exited!')
            exit()

        for each_file in file_list:
            line_index = 0
            for each_line in open(each_file):
                if line_index > 0:
                    each_line_split = each_line.strip().split('\t')
                    if len(each_line_split) == 9:
                        gene_id = each_line_split[0]
                        ko_d_id = each_line_split[4][2:]
                        if ko_d_id in interested_fun_set:
                            if ko_d_id not in fun_to_gene_dict:
                                fun_to_gene_dict[ko_d_id] = set()
                            fun_to_gene_dict[ko_d_id].add(gene_id)
                line_index += 1

    cmd_list_mafft = []
    cmd_list_trim  = []
    cmd_list_tree  = []
    cmd_1_mafft_txt_handle = open(cmd_1_mafft_txt, 'w')
    cmd_2_trim_txt_handle = open(cmd_2_trim_txt, 'w')
    cmd_3_tree_txt_handle = open(cmd_3_tree_txt, 'w')
    for each_fun in sorted(fun_to_gene_dict):

        # define file name
        fun_faa                     = '%s/%s.faa'           % (faa_dir, each_fun)
        current_gene_tree_dir       = '%s/%s'               % (tree_dir, each_fun)
        fun_aln                     = '%s/%s.aln'           % (aln_dir, each_fun)
        fun_aln_trimmed             = '%s/%s_trimal.aln'    % (trimmed_aln_dir, each_fun)
        if trim_with_bmge is True:
            fun_aln_trimmed         = '%s/%s_bmge.aln'      % (trimmed_aln_dir, each_fun)

        # extract sequences
        current_fun_gene_set = fun_to_gene_dict[each_fun]
        select_seq(combined_faa, current_fun_gene_set, fun_faa)

        os.system('mkdir %s' % current_gene_tree_dir)

        # prepare commands
        mafft_cmd      = 'mafft-einsi --thread %s --quiet %s > %s'      % (1, fun_faa, fun_aln)
        trim_cmd       = 'trimal -in %s -out %s -automated1'            % (fun_aln, fun_aln_trimmed)
        if trim_with_bmge is True:
            trim_cmd   = 'java -jar %s -i %s -m %s -t AA -h %s -of %s'  % (pwd_bmge_jar, fun_aln, trim_model, entropy_score_cutoff, fun_aln_trimmed)
        infer_tree_cmd = 'iqtree2 -s %s --seqtype AA -m %s -B 1000 --wbtl --bnni --prefix %s/%s -T %s --quiet' % (fun_aln_trimmed, iqtree_model, current_gene_tree_dir, each_fun, num_of_threads)

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
    print('Running mafft with %s cores for %s commands' % (num_of_threads, len(cmd_list_mafft)))
    pool = mp.Pool(processes=num_of_threads)
    pool.map(os.system, cmd_list_mafft)
    pool.close()
    pool.join()

    # run trim commands
    print('Trimming with %s cores for %s commands' % (num_of_threads, len(cmd_list_trim)))
    pool = mp.Pool(processes=num_of_threads)
    pool.map(os.system, cmd_list_trim)
    pool.close()
    pool.join()

    # run iqtree commands
    print('Running iqtree with %s cores' % num_of_threads)
    for each_iqtree_cmd in sorted(cmd_list_tree):
        print(each_iqtree_cmd)
        os.system(each_iqtree_cmd)


if __name__ == '__main__':

    koTree_parser = argparse.ArgumentParser()
    koTree_parser.add_argument('-i',         required=True,                          help='orthologous gene sequence')
    koTree_parser.add_argument('-fun',       required=True,                          help='interested functions')
    koTree_parser.add_argument('-cog',       required=False, default=None,           help='COG annotation results')
    koTree_parser.add_argument('-o',         required=True,                          help='output directory')
    koTree_parser.add_argument('-bmge',      required=False, action="store_true",    help='trim with BMGE, default is trimal')
    koTree_parser.add_argument('-bmge_m',    required=False, default='BLOSUM30',     help='trim model, default: BLOSUM30')
    koTree_parser.add_argument('-bmge_esc',  required=False, default='0.55',         help='entropy score cutoff, default: 0.55')
    koTree_parser.add_argument('-iqtree_m',  required=False, default='LG+G+I',       help='iqtree_model, default: LG+G+I')
    koTree_parser.add_argument('-f',         required=False, action="store_true",    help='force overwrite')
    koTree_parser.add_argument('-t',         required=False, type=int, default=1,    help='num of threads, default: 1')
    args = vars(koTree_parser.parse_args())
    koTree(args)


'''

cd /scratch/PI/ocessongwz/Sponge_r220/4_OMA_wd/OMA_wd/Output
TreeSAK FunTree -i /scratch/PI/ocessongwz/Sponge_r220/3_combined_genomes_50_5_dRep97_291.faa -fun K01995,K01996,K01997,K01998,K01999 -kegg /scratch/PI/ocessongwz/Sponge_r220/3_combined_genomes_50_5_dRep97_291_KEGG_wd -o interested_fun_tree_branched_chain_aa_transport_system -bmge -t 12 -f

'''
