import os
import argparse
from ete3 import Tree
from Bio import SeqIO


guide_tree_usage = '''
============================== guide_tree example command ==============================

TreeSAK guide_tree -i in.tree -g gnm_group.txt -o out.tree
TreeSAK guide_tree -i in.tree -g gnm_group.txt -o out.tree -id interested_gnm.txt
TreeSAK guide_tree -i in.tree -g gnm_group.txt -o out.tree -id representative_gene.fa
TreeSAK guide_tree -i in.tree -g gnm_group.txt -o out.tree -id representative_gene.aln

# Note
Topology tree leaves that do not have assigned genomes will be omitted from the output tree.

========================================================================================
'''

def sep_path_basename_ext(file_in):

    f_path, f_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'
    f_base, f_ext = os.path.splitext(f_name)
    f_ext = f_ext[1:]

    return f_name, f_path, f_base, f_ext


def guide_tree(args):

    tree_in                     = args['i']
    tree_in_leaves_to_ignore    = args['x']
    gnm_group_txt               = args['g']
    interested_gene_gnm_txt_fa  = args['id']
    tree_out                    = args['o']

    if tree_in_leaves_to_ignore is not None:
        if os.path.isfile(tree_in_leaves_to_ignore) is False:
            pass
        else:
            pass

    if os.path.isfile(tree_in) is False:
        print('%s not exist, program exited!' % tree_in)
        exit()

    count_all_gene_genome = True
    interested_gene_genome_set = set()
    if interested_gene_gnm_txt_fa is not None:

        if os.path.isfile(interested_gene_gnm_txt_fa) is False:
            print('%s not exist, program exited!' % interested_gene_gnm_txt_fa)
            exit()

        count_all_gene_genome = False
        _, _, _, id_f_ext = sep_path_basename_ext(interested_gene_gnm_txt_fa)
        if id_f_ext in ['fa', 'fasta','fas','ffn','faa','fna', 'aln']:
            for each_seq in SeqIO.parse(open(interested_gene_gnm_txt_fa), 'fasta'):
                interested_gene_genome_set.add(each_seq.id)
        else:
            for each_line in open(interested_gene_gnm_txt_fa):
                interested_gene_genome_set.add(each_line.strip().split('')[0])

    # define output/tmp file name
    _, _, f_in_base, f_in_ext            = sep_path_basename_ext(tree_in)
    _, f_out_path, f_out_base, f_out_ext = sep_path_basename_ext(tree_out)
    tree_in_subset_file                  = '%s/%s_subset.%s'                  % (f_out_path, f_in_base, f_in_ext)
    seqs_not_on_guide_tree_txt           = '%s/%s_seqs_not_on_guide_tree.txt' % (f_out_path, f_out_base)

    t_in = Tree(tree_in, format=0)
    topo_tree_leaves = [i.name for i in t_in.get_leaves()]

    seq_to_tax_dict = dict()
    grp_to_seq_dict = dict()
    seq_to_grp_dict = dict()
    seqs_not_on_topo_tree = set()
    for each_seq in open(gnm_group_txt):
        each_seq_split = each_seq.strip().split('\t')
        seq_id   = each_seq_split[0]

        count_current_seq = False
        if count_all_gene_genome is True:
            count_current_seq = True
        else:
            if seq_id in interested_gene_genome_set:
                count_current_seq = True

        if count_current_seq is True:
            seq_tax  = each_seq_split[1]
            seq_tax_split = seq_tax.strip().split(';')
            seq_to_tax_dict[seq_id] = seq_tax

            seq_grp = ''
            for each_r in seq_tax_split:
                if each_r in topo_tree_leaves:
                    seq_grp = each_r

            seq_to_grp_dict[seq_id] = seq_grp
            if seq_grp not in grp_to_seq_dict:
                grp_to_seq_dict[seq_grp] = set()
            grp_to_seq_dict[seq_grp].add(seq_id)

            if seq_grp == '':
                seqs_not_on_topo_tree.add(seq_id)

    if len(seqs_not_on_topo_tree) > 0:
        print('%s sequences can not be placed on guide tree, for details: %s'    % (len(seqs_not_on_topo_tree), seqs_not_on_guide_tree_txt))
        seqs_not_on_guide_tree_txt_handle = open(seqs_not_on_guide_tree_txt, 'w')
        for each_seq in sorted(list(seqs_not_on_topo_tree)):
            seqs_not_on_guide_tree_txt_handle.write('%s\t%s\n' % (each_seq, seq_to_tax_dict[each_seq]))
        seqs_not_on_guide_tree_txt_handle.close()

    # subset tree (if necessary)
    topo_tree_grp_with_seq = set()
    topo_tree_grp_without_seq = set()
    for leaf in topo_tree_leaves:
        if leaf in grp_to_seq_dict:
            topo_tree_grp_with_seq.add(leaf)
        else:
            topo_tree_grp_without_seq.add(leaf)

    # subset input topo tree
    topo_tree_to_use = tree_in
    if len(topo_tree_grp_without_seq) > 0:
        print('Removing the following leaves from %s, as no sequence in %s belongs to these groups: %s' % (tree_in, gnm_group_txt, ','.join(topo_tree_grp_without_seq)))
        subset_tree = t_in.copy()
        subset_tree.prune(topo_tree_grp_with_seq)
        subset_tree.write(outfile=tree_in_subset_file, format=9)
        topo_tree_to_use = tree_in_subset_file

    t_in = Tree(topo_tree_to_use, format=0)

    for leaf in t_in:
        leaf_name = leaf.name

        leaf_member = grp_to_seq_dict[leaf_name]
        leaf_member_str = ','.join(leaf_member)
        if len(leaf_member) >= 2:
            leaf_member_str = '(' + leaf_member_str + ')'
        if len(leaf_member) == 1:
            leaf.name = leaf_member_str
        else:
            leaf_p = leaf.up
            leaf_p.add_child(Tree((leaf_member_str + ';'), format=0))
            leaf_p.remove_child(leaf)
    t_in.write(outfile=tree_out, format=9)


if __name__ == '__main__':

    guide_tree_parser = argparse.ArgumentParser(usage=guide_tree_usage)
    guide_tree_parser.add_argument('-i',    required=True,                  help='input topo tree')
    guide_tree_parser.add_argument('-g',    required=True,                  help='gene/genome group/taxonomy')
    guide_tree_parser.add_argument('-id',   required=False, default=None,   help='interested genes/genomes')
    guide_tree_parser.add_argument('-o',    required=True,                  help='output tree')
    args = vars(guide_tree_parser.parse_args())
    guide_tree(args)
