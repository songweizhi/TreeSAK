import os
import argparse
from Bio import SeqIO
from ete3 import Tree
from distutils.spawn import find_executable


ALE1_usage = '''
============================================ ALE1 example commands ============================================

TreeSAK ALE1 -i OrthologousGroups.txt -s combined.faa -p oma -m 50 -jst 3 -f -o ALE1_op_dir

===============================================================================================================
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
