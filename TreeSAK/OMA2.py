import os
import glob
import argparse
from Bio import SeqIO


OMA2_usage = '''
============================== OMA2 example commands ==============================

TreeSAK OMA2 -i OrthologousGroups.txt -s OrthologousGroupsFasta -o op_dir -f -n 3
TreeSAK OMA2 -i OrthologousGroups.txt -s OrthologousGroupsFasta -o op_dir -f -c 85

===================================================================================
'''


def sep_path_basename_ext(file_in):

    f_path, f_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'
    f_base, f_ext = os.path.splitext(f_name)

    return f_name, f_path, f_base, f_ext[1:]


def get_gnm_og_cov(og_dir, og_ext, og_cov_txt):

    og_file_re   = '%s/*.%s' % (og_dir, og_ext)
    og_file_list = glob.glob(og_file_re)

    gnm_to_og_dict = dict()
    for og_file in og_file_list:
        _, _, og_id, _ = sep_path_basename_ext(og_file)
        for each_seq in SeqIO.parse(og_file, 'fasta'):
            seq_id = each_seq.id
            gnm_id = '_'.join(seq_id.split('_')[:-1])
            if gnm_id not in gnm_to_og_dict:
                gnm_to_og_dict[gnm_id] = set()
            gnm_to_og_dict[gnm_id].add(og_id)

    og_cov_txt_handle = open(og_cov_txt, 'w')
    for each_gnm in sorted(list(gnm_to_og_dict.keys())):
        gnm_og_set = gnm_to_og_dict[each_gnm]
        og_cov = len(gnm_og_set)*100/len(og_file_list)
        og_cov = float("{0:.2f}".format(og_cov))
        og_cov_txt_handle.write('%s\t%s\n' % (each_gnm, og_cov))
    og_cov_txt_handle.close()


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


def select_seq(seq_in, seq_id_list, seq_out):
    output_file_handle = open(seq_out, 'w')
    for seq_record in SeqIO.parse(seq_in, 'fasta'):
        if seq_record.id in seq_id_list:
            output_file_handle.write('>%s\n' % seq_record.id)
            output_file_handle.write('%s\n' % seq_record.seq)
    output_file_handle.close()


def OMA2(args):

    og_txt              = args['i']
    og_seq_dir          = args['s']
    gnm_txt             = args['g']
    op_dir              = args['o']
    force_overwrite     = args['f']
    min_gene_num        = args['n']
    min_gene_cov        = args['c']

    if (min_gene_num is None) and (min_gene_cov is None):
        print('Please specify either -c or -n, program exited!')
        exit()
    elif (min_gene_num is not None) and (min_gene_cov is not None):
        print('-c and -n are not compatible, program exited!')
        exit()

    og_txt_out = ''
    gnm_og_num_txt = ''
    filtered_seq_dir = ''
    if min_gene_num is not None:
        og_txt_out       = '%s/OrthologousGroups_num%s.txt'             % (op_dir, min_gene_num)
        gnm_og_num_txt   = '%s/OrthologousGroups_num%s_per_genome.txt'  % (op_dir, min_gene_num)
        filtered_seq_dir = '%s/OrthologousGroupsFasta_num%s'            % (op_dir, min_gene_num)
    if min_gene_cov is not None:
        og_txt_out       = '%s/OrthologousGroups_cov%s.txt'             % (op_dir, min_gene_cov)
        gnm_og_num_txt   = '%s/OrthologousGroups_cov%s_per_genome.txt'  % (op_dir, min_gene_cov)
        filtered_seq_dir = '%s/OrthologousGroupsFasta_cov%s'            % (op_dir, min_gene_cov)

    # check genome files
    interested_gnm_set = set()
    if gnm_txt is not None:
        if os.path.isfile(gnm_txt) is True:
            for each_gnm in open(gnm_txt):
                gnm_id = each_gnm.strip().split()[0]
                interested_gnm_set.add(gnm_id)
        else:
            print('%s not found, program exited!' % gnm_txt)
            exit()

    # create dir
    if os.path.isdir(op_dir) is True:
        if force_overwrite is True:
            os.system('rm -r %s' % op_dir)
        else:
            print('output folder detected, program exited!')
            exit()
    os.system('mkdir %s' % op_dir)
    os.system('mkdir %s' % filtered_seq_dir)

    # get overall genome set
    overall_gnm_set = set()
    for each_line in open(og_txt):
        if not each_line.startswith('#'):
            each_line_split = each_line.strip().split('\t')
            gene_list = each_line_split[1:]
            for each_gene in gene_list:
                gene_gnm = each_gene.split(':')[0]
                overall_gnm_set.add(gene_gnm)

    qualified_og_set = set()
    id_to_name_dict = dict()
    gene_to_extract_dict = dict()
    og_txt_out_handle = open(og_txt_out, 'w')
    for each_line in open(og_txt):
        if not each_line.startswith('#'):
            each_line_split = each_line.strip().split('\t')
            og_id = each_line_split[0]
            filename = 'OG%s' % int(og_id[3:])
            id_to_name_dict[og_id] = filename
            gene_list = each_line_split[1:]
            filtered_gene_set = set()
            for each_gene in gene_list:
                gene_gnm = each_gene.split(':')[0]
                gene_id = each_gene.split(':')[1].split(' ')[0]
                if len(interested_gnm_set) == 0:
                    filtered_gene_set.add(gene_id)
                else:
                    if gene_gnm in interested_gnm_set:
                        filtered_gene_set.add(gene_id)

            qualified_og = False
            if min_gene_num is not None:
                if len(filtered_gene_set) >= float(min_gene_num):
                    qualified_og = True
            if min_gene_cov is not None:
                if len(interested_gnm_set) == 0:
                    gnm_cov = len(filtered_gene_set)*100/len(overall_gnm_set)
                else:
                    gnm_cov = len(filtered_gene_set)*100/len(interested_gnm_set)

                if gnm_cov >= float(min_gene_cov):
                    qualified_og = True

            if qualified_og is True:
                qualified_og_set.add(og_id)
                og_txt_out_handle.write('%s\t%s\n' % (filename, ','.join(sorted(list(filtered_gene_set)))))
                gene_to_extract_dict[og_id] = filtered_gene_set
    og_txt_out_handle.close()

    for each_og in gene_to_extract_dict:
        seq_file_name = id_to_name_dict[each_og]
        source_file   = '%s/%s.fa' % (og_seq_dir, seq_file_name)
        filtered_file = '%s/%s.fa' % (filtered_seq_dir, seq_file_name)
        select_seq(source_file, gene_to_extract_dict[each_og], filtered_file)

    # get_gnm_og_cov
    get_gnm_og_cov(filtered_seq_dir, 'fa', gnm_og_num_txt)

    # report
    if min_gene_num is not None:
        print('The number of OG with genes >= %s is %s' % (min_gene_num, len(qualified_og_set)))
    if min_gene_cov is not None:
        print('The number of OG with coverage >= %s is %s' % (min_gene_cov, len(qualified_og_set)))

    print('Done!')


if __name__ == '__main__':

    OMA2_parser = argparse.ArgumentParser()
    OMA2_parser.add_argument('-i',   required=True,                         help='OrthologousGroups.txt')
    OMA2_parser.add_argument('-s',   required=True,                         help='sequence dir, OrthologousGroupsFasta')
    OMA2_parser.add_argument('-g',   required=False, default=None,          help='interested genomes')
    OMA2_parser.add_argument('-o',   required=True,  default=None,          help='output directory')
    OMA2_parser.add_argument('-n',   required=False, default=None,          help='minimal number of gene in a OG, not compatible with -c')
    OMA2_parser.add_argument('-c',   required=False, default=None,          help='minimal genome coverage cutoff, not compatible with -n')
    OMA2_parser.add_argument('-f',   required=False, action="store_true",   help='force overwrite')
    args = vars(OMA2_parser.parse_args())
    OMA2(args)
