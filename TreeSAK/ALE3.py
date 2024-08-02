import os
import argparse
import pandas as pd


ALE3_usage = '''
================= ALE3 example commands =================

TreeSAK ALE3 -2 ALE2_op_dir -o ALE3_op_dir_30 -f -c 30
TreeSAK ALE3 -2 ALE2_op_dir -o ALE3_op_dir_75 -f -c 75

# Needs the uml_rec files

=========================================================
'''


def ale_parser(rec_folder, SpeciesTreeRef_newick, TableInfo_tsv, TableEvents_tsv, GeneTrees_nwk):

    rec_files = [x for x in os.listdir(rec_folder) if x.endswith("uml_rec")]

    table_info = list()
    table_events = list()
    for rec_file in rec_files:
        with open(os.path.join(rec_folder, rec_file)) as f:
            fam = rec_file.replace(".ale.uml_rec", "")
            lines = f.readlines()
            stree = lines[2].strip()
            ll = lines[6].strip().split()[-1]
            dp, tp, lp = lines[8].strip().split("\t")[1:]
            n_reconciled_trees = int(lines[9].strip().split()[0])
            reconciled_trees = lines[11:n_reconciled_trees + 11]
            de, te, le, se = lines[11 + n_reconciled_trees + 1].split("\t")[1:]
            table = lines[11 + n_reconciled_trees + 3:]

        table_info.append((fam, ll, dp, tp, lp, de, te, le, se))
        table_events.append((fam, table))

    # write out SpeciesTreeRef.newick
    with open(SpeciesTreeRef_newick, "w") as f:
        f.write(stree.split("\t")[-1])

    # write out TableInfo.tsv
    with open(TableInfo_tsv, "w") as f:
        head = "\t".join(["Family", "LL", "Dp", "Tp", "Lp", "De", "Te", "Le", "Se"]) + "\n"
        f.write(head)
        for info in table_info:
            f.write("\t".join(info))

    # write out TableEvents.tsv
    with open(TableEvents_tsv, "w") as f:
        header = "Family\tBranchType\t" + table[0].replace("# of", "Branch")
        f.write(header)
        for fam, events in table_events:
            for b in events[1:]:
                f.write(fam + "\t" + b)

    # write out GeneTrees.nwk
    with open(GeneTrees_nwk, "w") as f:
        for t in reconciled_trees:
            f.write(t)


def get_verticality_and_transfer_propensity(TableEvents_tsv, verticality_txt, transfer_propensity_txt, fun_des_dict):

    df = pd.read_csv(TableEvents_tsv, sep="\t")
    dfb = df.groupby("Branch", as_index=False).sum()
    dff = df.groupby("Family").sum()

    dfb["Verticality"] = dfb["singletons"] / (dfb["singletons"] + dfb["Originations"] + dfb["Transfers"])
    dff["TransferPropensity"] = dff["Transfers"] / (dff["singletons"] + dff["Transfers"])

    verticality_dict = dfb.to_dict()['Verticality']
    transfer_propensity_dict = dff.to_dict()['TransferPropensity']

    with open(verticality_txt, 'w') as verticality_txt_handle:
        verticality_txt_handle.write('Branch\tVerticality\n')
        for each_key in sorted(list(verticality_dict.keys())):
            verticality_txt_handle.write('%s\t%s\n' % (each_key, verticality_dict[each_key]))

    with open(transfer_propensity_txt, 'w') as transfer_propensity_txt_handle:

        # write out header
        if len(fun_des_dict) == 0:
            transfer_propensity_txt_handle.write('OG\tTransfer_propensity\n')
        else:
            transfer_propensity_txt_handle.write('OG\tTransfer_propensity\tDescription\n')

        for each_key in sorted(list(transfer_propensity_dict.keys())):
            transfer_propensity = transfer_propensity_dict[each_key]
            transfer_propensity = float("{0:.3f}".format(transfer_propensity))
            each_key   = each_key.replace(('genome_tree.newick_'), '')
            each_key   = each_key.replace('.ufboot', '')
            if len(fun_des_dict) == 0:
                transfer_propensity_txt_handle.write('%s\t%s\n' % (each_key, transfer_propensity))
            else:
                transfer_propensity_txt_handle.write('%s\t%s\t%s\n' % (each_key, transfer_propensity, fun_des_dict.get(each_key, 'na')))


def ALE3(args):

    uml_rec_dir          = args['2']
    gene_presence_cutoff = args['c']
    fun_des_txt          = args['a']
    op_dir               = args['o']
    force_create_op_dir  = args['f']

    # read in fun_des_txt
    fun_des_dict = dict()
    if fun_des_txt is not None:
        if os.path.isfile(fun_des_txt) is True:
            for each_line in open(fun_des_txt):
                each_line_split = each_line.strip().split('\t')
                fun_des_dict[each_line_split[0]] = each_line_split[1]
        else:
            print('Specified "-a" not found, program exited!')
            exit()

    SpeciesTreeRef_newick   = '%s/SpeciesTreeRef.newick'    % op_dir
    TableInfo_tsv           = '%s/TableInfo.tsv'            % op_dir
    TableEvents_tsv         = '%s/TableEvents.tsv'          % op_dir
    GeneTrees_nwk           = '%s/GeneTrees.nwk'            % op_dir
    gene_content_dir        = '%s/GeneContent'              % op_dir
    gene_content_txt        = '%s/GeneContent.txt'          % op_dir
    verticality_txt         = '%s/Verticality.txt'          % op_dir
    transfer_propensity_txt = '%s/Transfer_propensity.txt'  % op_dir

    if os.path.isdir(op_dir) is True:
        if force_create_op_dir is True:
            os.system('rm -r %s' % op_dir)
        else:
            print('Output folder detected, program exited!')
            exit()
    os.system('mkdir %s' % op_dir)
    os.system('mkdir %s' % gene_content_dir)

    # parsing ALE2 outputs
    print('Parsing ALE2 outputs')
    ale_parser(uml_rec_dir, SpeciesTreeRef_newick, TableInfo_tsv, TableEvents_tsv, GeneTrees_nwk)

    # get_verticality_and_transfer_propensity
    print('Getting verticality and transfer propensity')
    get_verticality_and_transfer_propensity(TableEvents_tsv, verticality_txt, transfer_propensity_txt, fun_des_dict)

    # get genome content
    og_set = set()
    branch_to_og_dict = dict()
    col_index = {}
    for each_line in open(TableEvents_tsv):
        each_line_split = each_line.strip().split('\t')
        if each_line.startswith('Family'):
            col_index = {key: i for i, key in enumerate(each_line_split)}
        else:
            gene_family   = each_line_split[col_index['Family']]
            gene_family   = gene_family.replace(('genome_tree.newick_'), '')
            gene_family   = gene_family.replace('.ufboot', '')
            gene_branch   = each_line_split[col_index['Branch']]
            gene_presence = float(each_line_split[col_index['presence']])
            if gene_presence >= (gene_presence_cutoff/100):
                og_set.add(gene_family)
                if gene_branch not in branch_to_og_dict:
                    branch_to_og_dict[gene_branch] = set()
                branch_to_og_dict[gene_branch].add(gene_family)

    # write out gene content for each branch
    for each_branch in branch_to_og_dict:
        branch_gene_content = branch_to_og_dict[each_branch]

        current_gene_content_txt = '%s/%s.txt' % (gene_content_dir, each_branch)
        current_gene_content_txt_handle = open(current_gene_content_txt, 'w')
        for each_gene in sorted(list(branch_gene_content)):
            if len(fun_des_dict) == 0:
                current_gene_content_txt_handle.write('%s\n' % each_gene)
            else:
                current_gene_content_txt_handle.write('%s\t%s\n' % (each_gene, fun_des_dict.get(each_gene, 'na')))
        current_gene_content_txt_handle.close()

    og_list_sorted = sorted(list(og_set))

    gene_content_txt_handle = open(gene_content_txt, 'w')
    gene_content_txt_handle.write('Branch\t' + '\t'.join(og_list_sorted) + '\n')
    for each_gnm in sorted(list(branch_to_og_dict.keys())):
        og_pa_list = [each_gnm]
        for each_og in og_list_sorted:
            if each_og in branch_to_og_dict[each_gnm]:
                og_pa_list.append('1')
            else:
                og_pa_list.append('0')
        gene_content_txt_handle.write('\t'.join(og_pa_list) + '\n')
    gene_content_txt_handle.close()

    print('Protein families in GeneContent.txt: %s' % len(og_list_sorted))
    print('Genomes/branches in GeneContent.txt: %s' % len(branch_to_og_dict))


if __name__ == '__main__':

    ALE3_parser = argparse.ArgumentParser()
    ALE3_parser.add_argument('-2',   required=True,                             help='Folder with the uml_rec files')
    ALE3_parser.add_argument('-c',   required=False, type=float, default=75,   help='gene family presence cutoff in percentage, default: 75')
    ALE3_parser.add_argument('-a',   required=False, default=None,              help='OG functional description')
    ALE3_parser.add_argument('-o',   required=True,                             help='output dir')
    ALE3_parser.add_argument('-f',   required=False, action="store_true",       help='force overwrite')
    args = vars(ALE3_parser.parse_args())
    ALE3(args)
