import argparse
from ete3 import Tree


TaxonTree_usage = '''
================================ TaxonTree example commands ================================

TreeSAK TaxonTree -i ar53_r220.tree -tax o__Nitrososphaerales -o o__Nitrososphaerales.tree

============================================================================================
'''


def TaxonTree(args):

    tree_file_in     = args['i']
    interested_taxon = args['tax']
    tree_file_out    = args['o']

    input_tree = Tree(tree_file_in, quoted_node_names=True, format=1)

    matched_node_list = []
    for node in input_tree.traverse():
        if (node.name == interested_taxon) or (interested_taxon in node.name):
            matched_node_list.append(node.name)

    if len(matched_node_list) == 1:
        for node in input_tree.traverse():
            if node.name in matched_node_list:
                node.write(outfile=tree_file_out)
    else:
        print('There are multiple matched nodes. program exited!')
        print('Matched nodes: %s' % ','.join(matched_node_list))
        exit()

    print('Subset tree exported to: %s' % tree_file_out)
    print('Done!')


if __name__ == '__main__':

    TaxonTree_parser = argparse.ArgumentParser()
    TaxonTree_parser.add_argument('-i',   required=True,  help='input tree file')
    TaxonTree_parser.add_argument('-tax', required=True,  help='interested taxon')
    TaxonTree_parser.add_argument('-o',   required=True,  help='output tree file')
    args = vars(TaxonTree_parser.parse_args())
    TaxonTree(args)
