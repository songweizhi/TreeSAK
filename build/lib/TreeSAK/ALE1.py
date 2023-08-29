import os
import argparse
from Bio import SeqIO
from ete3 import Tree
import multiprocessing as mp
from distutils.spawn import find_executable


ALE1_usage = '''
=============================================== ALE1 example commands ===============================================

TreeSAK ALE1 -i OrthologousGroups.txt -s OrthologousGroups_combined.fasta -p oma -m 50 -t 6 -jt 3 -f -o ALE1_op_dir

=====================================================================================================================
'''


def subset_tree(tree_file_in, leaves_to_keep_list, tree_file_out):

    input_tree = Tree(tree_file_in)
    subset_tree = input_tree.copy()
    subset_tree.prune(leaves_to_keep_list, preserve_branch_length=True)
    if tree_file_out is None:
        return subset_tree.write()
    else:
        subset_tree.write(outfile=tree_file_out)


def get_ortho_to_gene_dict(ortho_groups_txt, og_program):

    ortho_to_gene_dict = dict()
    for each_og in open(ortho_groups_txt):
        if not each_og.startswith('#'):
            og_id = ''
            gene_list = []
            if og_program == 'orthofinder':
                each_og_split = each_og.strip().split(' ')
                og_id = each_og_split[0][:-1]
                gene_list = each_og_split[1:]
            elif og_program == 'oma':
                each_og_split = each_og.strip().split('\t')
                og_id = each_og_split[0]
                group_member_list = each_og_split[1:]
                for each_protein in group_member_list:
                    protein_id = each_protein.split(' ')[0].split(':')[1]
                    gene_list.append(protein_id)
            ortho_to_gene_dict[og_id] = gene_list

    return ortho_to_gene_dict


def prepare_ale_ip_worker(arg_list):

    qualified_og                        = arg_list[0]
    gene_tree_dir                       = arg_list[1]
    ale_wd                              = arg_list[2]
    genome_tree_file_rooted             = arg_list[3]
    gnm_pco_dict                        = arg_list[4]
    gene_tree_ufboot_for_ale            = arg_list[5]
    genome_tree_file_subset_for_ale     = arg_list[6]

    genome_tree_file_subset             = '%s_genome_tree.treefile' % qualified_og
    gene_tree_ufboot                    = '%s.ufboot'               % qualified_og
    gene_tree_treefile                  = '%s.treefile'             % qualified_og
    gene_tree_treefile_subset           = '%s_subset.treefile'      % qualified_og
    pwd_genome_tree_file_subset         = '%s/%s'                   % (gene_tree_dir, genome_tree_file_subset)
    pwd_genome_tree_file_subset_for_ale = '%s/%s'                   % (ale_wd, genome_tree_file_subset_for_ale)
    pwd_gene_tree_ufboot                = '%s/%s'                   % (gene_tree_dir, gene_tree_ufboot)
    pwd_gene_tree_ufboot_for_ale        = '%s/%s'                   % (ale_wd, gene_tree_ufboot_for_ale)
    pwd_gene_tree_treefile              = '%s/%s'                   % (gene_tree_dir, gene_tree_treefile)
    pwd_gene_tree_treefile_subset       = '%s/%s'                   % (gene_tree_dir, gene_tree_treefile_subset)

    # get genomes on gene tree
    gene_gnm_set = set()
    gnm_to_gene_dict = dict()
    for each_gene in Tree(pwd_gene_tree_treefile).get_leaf_names():

        gene_gnm = '_'.join(each_gene.split('_')[:-1])

        gene_gnm_set.add(gene_gnm)
        if gene_gnm not in gnm_to_gene_dict:
            gnm_to_gene_dict[gene_gnm] = {each_gene}
        else:
            gnm_to_gene_dict[gene_gnm].add(each_gene)

    # subset genome tree
    genome_tree_leaf_set = Tree(genome_tree_file_rooted).get_leaf_names()
    gnms_in_both_trees = set(genome_tree_leaf_set).intersection(gene_gnm_set)
    gnm_tree_subset_str = subset_tree(genome_tree_file_rooted, gnms_in_both_trees, None)
    gnm_tree_subset_str_for_ale = gnm_tree_subset_str
    gnm_tree_subset_str_for_ale = gnm_tree_subset_str_for_ale.replace('GCA_', 'GCA').replace('GCF_', 'GCF')

    # write out genome tree subset
    with open(pwd_genome_tree_file_subset, 'w') as pwd_genome_tree_file_subset_handle:
        pwd_genome_tree_file_subset_handle.write(gnm_tree_subset_str)

    # write out genome tree subset for running ALE
    with open(pwd_genome_tree_file_subset_for_ale, 'w') as pwd_genome_tree_file_subset_for_ale_handle:
        pwd_genome_tree_file_subset_for_ale_handle.write(gnm_tree_subset_str_for_ale)

    # get genes to keep in gene tree
    gene_set_to_keep = set()
    for each_gnm in gnms_in_both_trees:
        gene_set_to_keep.update(gnm_to_gene_dict.get(each_gnm, set()))

    # subset gene_tree.treefile
    subset_tree(pwd_gene_tree_treefile, gene_set_to_keep, pwd_gene_tree_treefile_subset)

    # subset gene_tree.ufboot and rename leaves for running ALE
    pwd_gene_tree_ufboot_for_ale_handle = open(pwd_gene_tree_ufboot_for_ale, 'w')
    for each_gene_tree in open(pwd_gene_tree_ufboot):
        gene_tree_str = each_gene_tree.strip()
        gene_tree_str_subset_for_ale = subset_tree(gene_tree_str, gene_set_to_keep, None)
        gene_tree_str_subset_for_ale = gene_tree_str_subset_for_ale.replace('GCA_', 'GCA').replace('GCF_', 'GCF')
        pwd_gene_tree_ufboot_for_ale_handle.write(gene_tree_str_subset_for_ale + '\n')
    pwd_gene_tree_ufboot_for_ale_handle.close()

    # get gene tree leaf name dict (for plot)
    leaf_name_dict = dict()
    for each_gene in Tree(pwd_gene_tree_treefile_subset).get_leaf_names():
        gene_id = each_gene
        gene_genome = '_'.join(gene_id.split('_')[:-1])
        genome_pco = gnm_pco_dict[gene_genome]
        gene_id_with_taxon = '%s_%s' % (genome_pco, gene_id.split('_')[-1])
        leaf_name_dict[gene_id] = gene_id_with_taxon


def ALE1(args):

    orthogroups_op_txt  = args['i']
    combined_faa        = args['s']
    og_program          = args['p']
    min_og_genome_num   = args['m']
    js_num_threads      = args['jst']
    force_create_op_dir = args['f']
    op_dir              = args['o']
    designate_ogs       = []
    to_ignore_ogs_list  = []

    # define output file name
    get_gene_tree_cmds_txt = '%s_cmds.txt' % op_dir

    iqtree_exe = ''
    if find_executable('iqtree2'):
        iqtree_exe = 'iqtree2'
    elif find_executable('iqtree'):
        iqtree_exe = 'iqtree'
    else:
        print('iqtree not detected, program exited!')
        exit()

    if os.path.isdir(op_dir) is True:
        if force_create_op_dir is True:
            os.system('rm -r %s' % op_dir)
        else:
            print('Output folder detected, program exited!')
            exit()

    if force_create_op_dir is True:
        if os.path.isdir(op_dir) is True:
            os.system('rm -r %s' % op_dir)
    os.system('mkdir %s' % op_dir)

    # get ortho_to_gene_dict
    ortho_to_gene_dict = get_ortho_to_gene_dict(orthogroups_op_txt, og_program)

    # get qualified orthogroups
    qualified_og_set = set()
    for each_ortho in ortho_to_gene_dict:
        ortho_gene_set = ortho_to_gene_dict[each_ortho]
        ortho_gnm_set = set()
        for each_gene in ortho_gene_set:
            gene_gnm = '_'.join(each_gene.split('_')[:-1])
            ortho_gnm_set.add(gene_gnm)
        if len(ortho_gnm_set) >= min_og_genome_num:
            qualified_og_set.add(each_ortho)
    print('The total number of identified orthogroups is %s.' % len(ortho_to_gene_dict))
    print('The number of orthogroups spanning >= %s genomes is %s.' % (min_og_genome_num, len(qualified_og_set)))

    # process qualified OG
    og_to_process = sorted([i for i in qualified_og_set])
    if len(designate_ogs) > 0:
        print('The number of designated OGs to process: %s' % len(designate_ogs))
        og_to_process = designate_ogs

    og_to_process_no_ignored = set()
    for each_og in og_to_process:
        if each_og not in to_ignore_ogs_list:
            og_to_process_no_ignored.add(each_og)

    # read sequence into dict
    gene_seq_dict = dict()
    for each_seq in SeqIO.parse(combined_faa, 'fasta'):
        seq_id = each_seq.id
        gene_seq_dict[seq_id] = str(each_seq.seq)

    # extract gene sequences and prepare commands for building gene tree
    print('Preparing commands and sequence files for building gene trees')
    get_gene_tree_cmds_txt_handle = open(get_gene_tree_cmds_txt, 'w')
    for qualified_og in sorted(og_to_process_no_ignored):
        qualified_og_gene_set = ortho_to_gene_dict[qualified_og]
        qualified_og_gene_faa = '%s/%s.faa' % (op_dir, qualified_og)

        # write out commands
        mafft_cmd  = 'mafft-einsi --thread %s --quiet %s.faa > %s.aln'       % (js_num_threads, qualified_og, qualified_og)
        iqtree_cmd = '%s -m LG+G+I -bb 1000 --wbtl -nt %s -s %s.aln -pre %s' % (iqtree_exe, js_num_threads, qualified_og, qualified_og)
        get_gene_tree_cmds_txt_handle.write('%s; %s\n' % (mafft_cmd, iqtree_cmd))

        # write out sequences
        qualified_og_gene_faa_handle = open(qualified_og_gene_faa, 'w')
        for each_gene in qualified_og_gene_set:
            qualified_og_gene_faa_handle.write('>%s\n' % each_gene)
            qualified_og_gene_faa_handle.write('%s\n' % gene_seq_dict[each_gene])
        qualified_og_gene_faa_handle.close()
    get_gene_tree_cmds_txt_handle.close()

    print('Done!')

if __name__ == '__main__':

    ALE1_parser = argparse.ArgumentParser()
    ALE1_parser.add_argument('-i',   required=True,                         help='orthologous groups, either from orthofinder or oma')
    ALE1_parser.add_argument('-s',   required=True,                         help='sequence file, e.g., combined.faa')
    ALE1_parser.add_argument('-p',   required=True,                         help='orthologous identification program, orthofinder or oma')
    ALE1_parser.add_argument('-m',   required=False, type=int, default=50,  help='min_og_genome_num, default: 50')
    ALE1_parser.add_argument('-o',   required=True,                         help='output dir, i.e., OMA working directory')
    ALE1_parser.add_argument('-jst', required=False, type=int, default=3,   help='number of threads for job script, default: 3')
    ALE1_parser.add_argument('-f',   required=False, action="store_true",   help='force overwrite')
    args = vars(ALE1_parser.parse_args())
    ALE1(args)
