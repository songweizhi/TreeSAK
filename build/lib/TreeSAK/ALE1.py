import os
import glob
import argparse
from Bio import SeqIO
from ete3 import Tree
from distutils.spawn import find_executable


ALE1_usage = '''
====================================== ALE1 example commands ======================================

TreeSAK ALE1 -i OrthologousGroups.txt -s combined.faa -p oma -m 50 -jst 3 -f -o ALE1_op_dir -bmge
TreeSAK ALE1 -ms s03_marker_seq -msx fa -p marker_set_1 -m 50 -jst 3 -f -o ALE1_op_dir -bmge

===================================================================================================
'''


def check_dependencies(program_list):

    not_detected_programs = []
    for needed_program in program_list:
        if find_executable(needed_program) is None:
            not_detected_programs.append(needed_program)

    if not_detected_programs != []:
        print('%s not found, program exited!' % ','.join(not_detected_programs))
        exit()


def sep_path_basename_ext(file_in):

    f_path, file_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'

    f_base, f_ext = os.path.splitext(file_name)

    return f_path, f_base, f_ext


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

    orthogroups_op_txt          = args['i']
    combined_faa                = args['s']
    og_program                  = args['p']
    marker_seq_dir              = args['ms']
    marker_seq_ext              = args['msx']
    min_og_genome_num           = args['m']
    js_num_threads              = args['jst']
    force_create_op_dir         = args['f']
    op_dir                      = args['o']
    trim_with_bmge              = args['bmge']
    bmge_trim_model             = args['bmge_m']
    bmge_entropy_score_cutoff   = args['bmge_esc']
    designate_ogs               = []
    to_ignore_ogs_list          = []

    # check dependencies
    check_dependencies(['java', 'blastp', 'mafft-einsi'])

    # specify path to BMGE.jar
    current_file_path = '/'.join(os.path.realpath(__file__).split('/')[:-1])
    pwd_bmge_jar      = '%s/BMGE.jar' % current_file_path

    # define output file name
    get_gene_tree_cmds_txt = '%s_cmds.txt' % op_dir

    # determine the version of iqtree available on the system
    if find_executable('iqtree2'):
        iqtree_exe = 'iqtree2'
    elif find_executable('iqtree'):
        iqtree_exe = 'iqtree'
    else:
        print('iqtree not detected, program exited!')
        exit()

    # create op_dir
    if os.path.isdir(op_dir) is True:
        if force_create_op_dir is True:
            os.system('rm -r %s' % op_dir)
        else:
            print('Output folder detected, program exited!')
            exit()
    os.system('mkdir %s' % op_dir)

    if (orthogroups_op_txt is not None) and (combined_faa is not None) and (marker_seq_dir is None):

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

            og_aln          = '%s.aln'          % qualified_og
            og_aln_trimmed  = '%s_trimmed.aln'  % qualified_og

            # write out commands
            mafft_cmd       = 'mafft-einsi --thread %s --quiet %s.faa > %s'         % (js_num_threads, qualified_og, og_aln)
            trim_cmd        = 'java -jar %s -i %s -m %s -t AA -h %s -of %s'         % (pwd_bmge_jar, og_aln, bmge_trim_model, bmge_entropy_score_cutoff, og_aln_trimmed)
            iqtree_cmd      = '%s -m LG+G+I -bb 1000 --wbtl -nt %s -s %s -pre %s'   % (iqtree_exe, js_num_threads, og_aln, qualified_og)
            if trim_with_bmge is True:
                iqtree_cmd  = '%s -m LG+G+I -bb 1000 --wbtl -nt %s -s %s -pre %s'   % (iqtree_exe, js_num_threads, og_aln_trimmed, qualified_og)

            if trim_with_bmge is True:
                get_gene_tree_cmds_txt_handle.write('%s; %s; %s\n' % (mafft_cmd, trim_cmd, iqtree_cmd))
            else:
                get_gene_tree_cmds_txt_handle.write('%s; %s\n' % (mafft_cmd, iqtree_cmd))

            # write out sequences
            qualified_og_gene_faa_handle = open(qualified_og_gene_faa, 'w')
            for each_gene in qualified_og_gene_set:
                qualified_og_gene_faa_handle.write('>%s\n' % each_gene)
                qualified_og_gene_faa_handle.write('%s\n' % gene_seq_dict[each_gene])
            qualified_og_gene_faa_handle.close()
        get_gene_tree_cmds_txt_handle.close()

    elif (orthogroups_op_txt is None) and (combined_faa is None) and (marker_seq_dir is not None):

        marker_seq_re = '%s/*.%s' % (marker_seq_dir, marker_seq_ext)
        marker_seq_list = glob.glob(marker_seq_re)

        marker_to_gene_dict = dict()
        for each_file in marker_seq_list:
            _, f_base, _ = sep_path_basename_ext(each_file)
            marker_to_gene_dict[f_base] = set()
            for each_seq in SeqIO.parse(each_file, 'fasta'):
                marker_to_gene_dict[f_base].add(each_seq.id)

        # get qualified orthogroups
        qualified_og_set = set()
        for each_ortho in marker_to_gene_dict:
            ortho_gene_set = marker_to_gene_dict[each_ortho]
            ortho_gnm_set = set()
            for each_gene in ortho_gene_set:
                gene_gnm = '_'.join(each_gene.split('_')[:-1])
                ortho_gnm_set.add(gene_gnm)
            if len(ortho_gnm_set) >= min_og_genome_num:
                qualified_og_set.add(each_ortho)
        print('The total number of identified orthogroups is %s.'       % len(marker_to_gene_dict))
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

        # extract gene sequences and prepare commands for building gene tree
        print('Preparing commands for building gene trees')
        get_gene_tree_cmds_txt_handle = open(get_gene_tree_cmds_txt, 'w')
        for qualified_og in sorted(og_to_process_no_ignored):

            # copy sequence file into output directory
            os.system('cp %s/%s.%s %s/' % (marker_seq_dir, qualified_og, marker_seq_ext, op_dir))

            qualified_og_aln          = '%s.aln'          % qualified_og
            qualified_og_aln_trimmed  = '%s_trimmed.aln'  % qualified_og

            # write out commands
            mafft_cmd  = 'mafft-einsi --thread %s --quiet %s.%s > %s'               % (js_num_threads, qualified_og, marker_seq_ext, qualified_og_aln)
            trim_cmd   = 'java -jar %s -i %s -m %s -t AA -h %s -of %s'              % (pwd_bmge_jar, qualified_og_aln, bmge_trim_model, bmge_entropy_score_cutoff, qualified_og_aln_trimmed)
            iqtree_cmd = '%s -m LG+G+I -bb 1000 --wbtl -nt %s -s %s -pre %s'        % (iqtree_exe, js_num_threads, qualified_og_aln, qualified_og)
            if trim_with_bmge is True:
                iqtree_cmd = '%s -m LG+G+I -bb 1000 --wbtl -nt %s -s %s -pre %s'    % (iqtree_exe, js_num_threads, qualified_og_aln_trimmed, qualified_og)

            if trim_with_bmge is False:
                get_gene_tree_cmds_txt_handle.write('%s; %s\n' % (mafft_cmd, iqtree_cmd))
            else:
                get_gene_tree_cmds_txt_handle.write('%s; %s; %s\n' % (mafft_cmd, trim_cmd, iqtree_cmd))

        get_gene_tree_cmds_txt_handle.close()

    print('Sequece files exported to %s.'                       % op_dir)
    print('Commands for inferring gene tree exported to %s.'    % get_gene_tree_cmds_txt)
    print('Done!')


if __name__ == '__main__':

    ALE1_parser = argparse.ArgumentParser()
    ALE1_parser.add_argument('-i',          required=False, default=None,          help='orthologous groups, either from orthofinder or oma')
    ALE1_parser.add_argument('-s',          required=False, default=None,          help='sequence file, e.g., combined.faa')
    ALE1_parser.add_argument('-ms',         required=False, default=None,          help='input is a folder holds the sequence of each marker')
    ALE1_parser.add_argument('-msx',        required=False, default='fa',          help='file extension of marker sequence file, default: fa')
    ALE1_parser.add_argument('-p',          required=True,                         help='orthologous identification program, orthofinder or oma')
    ALE1_parser.add_argument('-m',          required=False, type=int, default=50,  help='min_og_genome_num, default: 50')
    ALE1_parser.add_argument('-bmge',       required=False, action="store_true",   help='trim MSA with BMGE, default no trimming')
    ALE1_parser.add_argument('-bmge_m',     required=False, default='BLOSUM30',    help='BMGE trim model, default: BLOSUM30')
    ALE1_parser.add_argument('-bmge_esc',   required=False, default='0.55',        help='BMGE entropy score cutoff, default: 0.55')
    ALE1_parser.add_argument('-o',          required=True,                         help='output dir, i.e., OMA working directory')
    ALE1_parser.add_argument('-jst',        required=False, type=int, default=3,   help='number of threads specified in job script, default: 3')
    ALE1_parser.add_argument('-f',          required=False, action="store_true",   help='force overwrite')
    args = vars(ALE1_parser.parse_args())
    ALE1(args)
