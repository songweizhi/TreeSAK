import os
import argparse


ALE7_usage = '''
=============================== ALE7 example commands ===============================

# get presence/absence of interested functions in reconstructed ancestral genomes

TreeSAK ALE7 -6 ALE6_op_dir -fun ko.txt -node internal_node.txt -o Fun_PA.txt
TreeSAK ALE7 -6 ALE6_op_dir -fun K01995,K01995 -node 359,466,470 -o Fun_PA.txt
TreeSAK ALE7 -6 ALE6_op_dir -fun arCOG07811,K01995 -node 359,466,470 -o Fun_PA.txt

# needed input files:
-6: annotation_COG/annotation_KEGG

=====================================================================================
'''


def ALE7(args):

    ale6_op_dir        = args['6']
    interested_fun_txt = args['fun']
    interested_gnm_txt = args['node']
    op_txt             = args['o']

    if os.path.isdir(ale6_op_dir) is False:
        print('%s not found, program exited!' % ale6_op_dir)
        exit()

    interested_fun_set = set()
    if os.path.isfile(interested_fun_txt) is False:
        if ',' in interested_fun_txt:
            interested_fun_set = interested_fun_txt.split(',')
        else:
            interested_fun_set.add(interested_fun_txt)
    else:
        for each_fun in open(interested_fun_txt):
            interested_fun_set.add(each_fun.strip().split()[0])

    interested_node_set = set()
    if os.path.isfile(interested_gnm_txt) is False:
        if ',' in interested_gnm_txt:
            interested_node_set = interested_gnm_txt.split(',')
        else:
            interested_node_set.add(interested_gnm_txt)
    else:
        for each_node in open(interested_gnm_txt):
            interested_node_set.add(each_node.strip().split()[0])

    interested_fun_list_sorted = sorted(list(interested_fun_set))
    interested_gnm_list_sorted = sorted(list(interested_node_set))

    op_txt_handle = open(op_txt, 'w')
    op_txt_handle.write('ID\t' + '\t'.join(interested_fun_list_sorted) + '\n')
    for each_node in interested_gnm_list_sorted:
        node_annotation_cog  = '%s/annotation_COG/%s_COG.txt'   % (ale6_op_dir, each_node)
        node_annotation_kegg = '%s/annotation_KEGG/%s_KEGG.txt' % (ale6_op_dir, each_node)

        node_fun_set = set()
        if os.path.isfile(node_annotation_cog):
            for each_line in open(node_annotation_cog):
                node_fun_set.add(each_line.strip().split()[0])
        if os.path.isfile(node_annotation_kegg):
            for each_line in open(node_annotation_kegg):
                node_fun_set.add(each_line.strip().split()[0])

        fun_pa_list = [each_node]
        for each_fun in interested_fun_list_sorted:
            if each_fun in node_fun_set:
                fun_pa_list.append('1')
            else:
                fun_pa_list.append('0')

        op_txt_handle.write('\t'.join(fun_pa_list) + '\n')
    op_txt_handle.close()

    print('Results exported to %s' % op_txt)
    print('Done!')


if __name__ == '__main__':

    ALE7_parser = argparse.ArgumentParser()
    ALE7_parser.add_argument('-6',      required=True,                          help='ALE6 output directory')
    ALE7_parser.add_argument('-fun',    required=True,                          help='interested functions')
    ALE7_parser.add_argument('-node',   required=True,                          help='interested internal nodes')
    ALE7_parser.add_argument('-o',      required=True,                          help='output directory')
    args = vars(ALE7_parser.parse_args())
    ALE7(args)
