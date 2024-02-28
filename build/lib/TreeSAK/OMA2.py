import os
import argparse
from Bio import SeqIO


OMA2_usage = '''
============================ OMA2 example commands ============================

TreeSAK OMA2 -i OrthologousGroups.txt -s OrthologousGroupsFasta -o op_dir -f 

===============================================================================
'''


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

    og_txt_out          = '%s/OrthologousGroups_min%s.txt'  % (op_dir, min_gene_num)
    filtered_seq_dir    = '%s/OrthologousGroupsFasta_min%s' % (op_dir, min_gene_num)

    # check genome files
    interested_gnm_set = set()
    if gnm_txt is not None:
        if os.path.isfile(gnm_txt) is True:
            for each_gnm in open(gnm_txt):
                interested_gnm_set.add(each_gnm.strip())
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
            if len(filtered_gene_set) > min_gene_num:
                og_txt_out_handle.write('%s\t%s\n' % (filename, ','.join(sorted(list(filtered_gene_set)))))
                gene_to_extract_dict[og_id] = filtered_gene_set
    og_txt_out_handle.close()

    for each_og in gene_to_extract_dict:
        seq_file_name = id_to_name_dict[each_og]
        source_file   = '%s/%s.fa' % (og_seq_dir, seq_file_name)
        filtered_file = '%s/%s.fa' % (filtered_seq_dir, seq_file_name)
        select_seq(source_file, gene_to_extract_dict[each_og], filtered_file)

    print('Done!')


if __name__ == '__main__':

    OMA2_parser = argparse.ArgumentParser()
    OMA2_parser.add_argument('-i',   required=True,                          help='OrthologousGroups.txt')
    OMA2_parser.add_argument('-s',   required=True,                          help='sequence dir, OrthologousGroupsFasta')
    OMA2_parser.add_argument('-g',   required=False, default=None,           help='id of interested genomes')
    OMA2_parser.add_argument('-o',   required=True,  default=None,           help='output directory')
    OMA2_parser.add_argument('-n',   required=False, type=int, default=3,    help='minimal number of gene in a OG, deafult: 3')
    OMA2_parser.add_argument('-f',   required=False, action="store_true",    help='force overwrite')
    args = vars(OMA2_parser.parse_args())
    OMA2(args)
