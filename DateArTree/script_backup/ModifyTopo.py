from ete3 import Tree
from os.path import join, dirname, exists


def get_topology_without_alpha(st, two_nodes_1, two_nodes_2):

    _st = st.copy()
    n_replaced = _st.get_common_ancestor(two_nodes_1)

    if len(n_replaced.get_leaf_names()) == 10:
        n_p = n_replaced.up
        n_p_p = n_p.up   # All proteobacteria
        n_p.remove_child(n_replaced)
        return _st, n_p_p, n_p
    else:
        beta_gamma = _st.get_common_ancestor(two_nodes_2)
        for _ in n_replaced.children[::]:
            n_replaced.remove_child(_)
        n_replaced.add_child(beta_gamma)
        return _st, all_pro, n_replaced


def read_tree(in_tree, format=None):

    if isinstance(in_tree, str) and exists(in_tree):
        if format=='auto':
            for f in [0,1,2,3,4,5]:
                try:
                    t = Tree(in_tree, format=f)
                    return t
                except:
                    pass
        else:
            t = Tree(open(in_tree).read(), format=format)
    elif isinstance(in_tree, Tree):
        t = in_tree
    else:
        raise IOError('unknown input')
    return t


def erase_name(in_tree_file, format=0):

    t = read_tree(in_tree_file, format=format)
    for n in t.traverse():
        if not n.is_leaf():
            n.name = ''
    return t


def get_mito(dataset):

    euk_tree = Tree(euk_reference_tree, 8)
    euk_tree.prune(dataset)
    et = erase_name(euk_tree)
    return {"Mito": et}


########################################################################################################################

intree                  = f'./dating/topology/mixture_models/deno100/phy/deno100.final_TP1.newick'
euk_reference_tree      = '/mnt/home-backup/thliao/AOB/analysis/update_calibrations/mito_dating/phylo/manual_topology/euk.tre'
euk_list_txt            = '/mnt/home-backup/thliao/AOB/analysis/update_calibrations/mito_dating/euk.list'
base_odir               = "./dating/topology"

Rickettsiales_lineage   = ['GCA_008189685.1', 'GCA_003015145.1']
Magneto_lineage         = ["GCA_000014865.1", "GCA_002109495.1"]
remaining_alpha         = ['GCA_000264455.2', 'GCA_002924445.1']
Holo                    = "GCA_000469665.2"
two_nodes_1             = ['GCA_000014865.1', 'GCA_000264455.2']
two_nodes_2             = ['GCA_018655245.1', 'GCA_002356115.1']

TP_dict                 = {"TP1": "((((Holo,other_alpha),Rick),Mito),Magneto);",
                           "TP2": "(((Holo,other_alpha),(Rick,Mito)),Magneto);",
                           "TP3": "((((Holo,Rick),other_alpha),Mito),Magneto);",
                           "TP4": "((((Mito,Rick),Holo),other_alpha),Magneto);"}

# copy from '/home-user/sswang/project/Mito/results/euk_tree/euk.tre'. The Fig S2A in wang 2021 NC
# rephrase Porphyra purpurea into Porphyra umbilicalis
# add Polysphondylium pallidum manually

'''

1. provide a input tree
2. provide a set of tree skeleton
3. use two leaves to determine a clade.

'''

########################################################################################################################


for tp_name, TP in TP_dict.items():

    print('%s\t%s' % (tp_name, TP))

    # get internal node to tree string dict

    st = Tree(intree, 8)
    _st, all_pro, n_replaced = get_topology_without_alpha(st.copy(), two_nodes_1, two_nodes_2)

    mito_usage = [_ for _ in open(euk_list_txt).read().split('\n')]
    m_dict = get_mito(mito_usage)

    nodes_dict = {"Rick"        : st.get_common_ancestor(Rickettsiales_lineage),
                  "other_alpha" : st.get_common_ancestor(remaining_alpha),
                  "Magneto"     : st.get_common_ancestor(Magneto_lineage),
                  "Holo"        : [_ for _ in st.traverse() if _.name == Holo][0],
                  "Mito"        : m_dict['Mito']}

    for k, n in nodes_dict.items():
        TP = TP.replace(k, n.write(format=3).strip(';'))
    n_replaced.add_child(Tree(TP, format=3))

