import os
import argparse
from ete3 import Tree


iTOL_gene_tree_usage = '''
====================== iTOL_gene_tree example commands ======================

TreeSAK iTOL_gene_tree -tree genes.tree -i gnm_taxon.txt -o gene_taxon.txt
TreeSAK iTOL_gene_tree -txt gene_id.txt -i gnm_taxon.txt -o gene_taxon.txt

=============================================================================
'''


def iTOL_gene_tree(args):

    input_tree_file = args['tree']
    input_txt_file  = args['txt']
    meta_txt        = args['i']
    op_txt          = args['o']
    include_na      = args['na']

    if (input_tree_file is None) and (input_txt_file is None):
        print('Please provide gene id with at least one approach, program exited!')
        exit()

    if os.path.isfile(meta_txt) is False:
        print('Metadata file not found, program exited!')
        exit()

    metadata_dict = dict()
    for each_gnm in open(meta_txt):
        each_gnm_split = each_gnm.strip().split('\t')
        if len(each_gnm_split) == 2:
            gnm_id = each_gnm_split[0]
            meta_value = each_gnm_split[1]
            metadata_dict[gnm_id] = meta_value

    gene_id_set = set()
    if input_tree_file is not None:
        if os.path.isfile(input_tree_file) is False:
            print('Tree file not found, program exited!')
            exit()
        else:
            for leaf in Tree(input_tree_file, format=1):
                gene_id_set.add(leaf.name)

    if input_txt_file is not None:
        if os.path.isfile(input_txt_file) is False:
            print('Txt file not found, program exited!')
            exit()
        else:
            for each_id in open(input_txt_file):
                gene_id_set.add(each_id.strip())

    op_txt_handle = open(op_txt, 'w')
    for gene_id in gene_id_set:
        gnm_id = '_'.join(gene_id.split('_')[:-1])
        gnm_meta = metadata_dict.get(gnm_id, 'na')
        if include_na is True:
            op_txt_handle.write('%s\t%s\n' % (gene_id, gnm_meta))
        else:
            if gnm_meta != 'na':
                op_txt_handle.write('%s\t%s\n' % (gene_id, gnm_meta))
    op_txt_handle.close()

    print('Done!')


if __name__ == '__main__':

    iTOL_gene_tree_parser = argparse.ArgumentParser(usage=iTOL_gene_tree_usage)
    iTOL_gene_tree_parser.add_argument('-i',     required=True,                         help='input metadata')
    iTOL_gene_tree_parser.add_argument('-tree',  required=False, default=None,          help='gene id, in tree file')
    iTOL_gene_tree_parser.add_argument('-txt',   required=False, default=None,          help='gene id, in txt file')
    iTOL_gene_tree_parser.add_argument('-o',     required=True,                         help='output metadata')
    iTOL_gene_tree_parser.add_argument('-na',    required=False, action='store_true',   help='include leaves with na values')
    args = vars(iTOL_gene_tree_parser.parse_args())
    iTOL_gene_tree(args)
