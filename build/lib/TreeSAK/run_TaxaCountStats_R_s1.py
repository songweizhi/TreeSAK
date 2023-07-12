from __future__ import print_function
import os
import re
import sys
import numpy
import argparse
from ete3 import Tree
import multiprocessing as mp
from operator import itemgetter
from collections import defaultdict


def get_rename_dict(tree_str_in, mag_rename_dict, mag_cluster_dict, sponge_mag_tax_dict, gtdb_gnm_tax_dict):

    # rename dict: {'old_name':'new_name'}

    leaf_rename_dict = {}
    for leaf in Tree(tree_str_in, format=1):

        leaf_name_gnm = '_'.join(leaf.name.split('_')[:-1])
        leaf_name_gnm = mag_rename_dict.get(leaf_name_gnm, leaf_name_gnm)
        leaf_cluster = mag_cluster_dict.get(leaf_name_gnm, 'cluster_0')

        leaf_name_gnm_no_source = leaf_name_gnm
        if '.gtdb' in leaf_name_gnm_no_source:
            leaf_name_gnm_no_source = leaf_name_gnm[:-5]
        if '.ncbi' in leaf_name_gnm:
            leaf_name_gnm_no_source = leaf_name_gnm[:-5]

        # get mag_taxon_str
        gnm_taxon_str = 'NA'
        if leaf_name_gnm_no_source in sponge_mag_tax_dict:
            gnm_taxon_str = sponge_mag_tax_dict[leaf_name_gnm_no_source]
        if leaf_name_gnm_no_source in gtdb_gnm_tax_dict:
            gnm_taxon_str = gtdb_gnm_tax_dict[leaf_name_gnm_no_source]

        # get mag_taxon_str (GCA GCF things)
        if gnm_taxon_str == 'NA':
            mag_id_no_ext_no_source_GCF = leaf_name_gnm_no_source.replace('GCA', 'GCF')
            if mag_id_no_ext_no_source_GCF in gtdb_gnm_tax_dict:
                gnm_taxon_str = gtdb_gnm_tax_dict[mag_id_no_ext_no_source_GCF]

        gnm_taxon_str_no_space = gnm_taxon_str.replace(' ', '_')
        gnm_taxon_str_no_space = gnm_taxon_str_no_space.replace(';', '|')
        leaf_name_new = '%s|%s|strain__%s' % (leaf_cluster, gnm_taxon_str_no_space, '_'.join(leaf.name.split('_')[:-1]))

        leaf_rename_dict[leaf.name] = leaf_name_new

    return leaf_rename_dict


def rename_tree(tree_str_in, rename_dict):

    t_in = Tree(tree_str_in, format=1)
    for leaf in t_in:
        leaf_name = leaf.name
        leaf_name_new = rename_dict.get(leaf_name, leaf_name)
        leaf.name = leaf_name_new

    return t_in.write()


def gtdb_gnm_metadata_parser(gtdb_genome_metadata):

    genome_to_taxon_dict = {}
    genome_to_completeness_dict = {}
    genome_to_contamination_dict = {}
    genome_to_biosample_dict = {}

    col_index = {}
    for each_ref in open(gtdb_genome_metadata):
        each_ref_split = each_ref.strip().split('\t')
        if each_ref.startswith('accession'):
            col_index = {key: i for i, key in enumerate(each_ref_split)}
        else:
            ref_accession = each_ref_split[0][3:]
            gnm_completeness = float(each_ref_split[2])
            gnm_contamination = float(each_ref_split[3])
            gtdb_taxon = each_ref_split[col_index['gtdb_taxonomy']]
            ncbi_biosample = each_ref_split[col_index['ncbi_biosample']]
            genome_to_taxon_dict[ref_accession] = gtdb_taxon
            genome_to_completeness_dict[ref_accession] = gnm_completeness
            genome_to_contamination_dict[ref_accession] = gnm_contamination
            genome_to_biosample_dict[ref_accession] = ncbi_biosample

    return genome_to_completeness_dict, genome_to_contamination_dict, genome_to_taxon_dict, genome_to_biosample_dict


def parse_taxonomy(taxon_name):  # given a taxon name, try to return whatever taxonomic info is available as a list starting with the highest level classification and going lower (or a map?)
    #name_elements = re.split("\|", taxon_name)
    name_elements = taxon_name.split('|')
    #print('name_elements')
    #print(name_elements)

    if (len(name_elements) < 8) or (len(name_elements) > 9):
        print("Nonstandard!")
        quit()

    name_map = {}
    name_map['cluster'] = name_elements[0]
    name_map['domain'] = name_elements[1]
    name_map['phylum'] = name_elements[2]
    name_map['class'] = name_elements[3]
    name_map['order'] = name_elements[4]
    name_map['family'] = name_elements[5]
    name_map['genus'] = name_elements[6]
    name_map['species'] = name_elements[7]
    if len(name_elements) == 9:
        name_map['ncbi_id'] = name_elements[8]
    return name_map


def summarize_taxonomy(name_list, tax_level, name_to_tax_dict):  # take a list of names from a clade and summarize taxonomic info (labels and their frequencies)
    total_size = len(name_list)  # it perhaps makes sense to normalize by the size of the clade
    breakdown = {}
    for name in name_list:
        info = name_to_tax_dict[name]
        if info[tax_level] in breakdown:
            breakdown[info[tax_level]] += 1.0 / float(total_size)
        else:
            breakdown[info[tax_level]] = 1.0 / float(total_size)
    return breakdown


def count_sister_taxa(target_label, tree_in_ml, tree_in_bs, output_file):

    # read in argument
    # target_label = args['l']
    # tree_in_ml   = args['ml']
    # tree_in_bs   = args['bs']
    # output_file  = args['out']

    # edit target_label to make the comparisons at a desired taxonomic level

    # compute the most frequent sister group of each (monophyletic?) group on the tree, to identify trends in gene transfers, "unstable" taxa, etc.

    # read the ML tree, set up the taxonomy stuff, and calculate the number of clades per label, and the sizes of those clades (to report at the end)
    labels = {}
    name_to_tax_info = defaultdict(dict)
    all_tree_leaf_names = []
    ml_tree = Tree(tree_in_ml)  # note that ete3 treats this input tree as rooted
    for leaf in ml_tree:
        taxonomy = parse_taxonomy(leaf.name)
        name_to_tax_info[leaf.name] = taxonomy
        all_tree_leaf_names.append(leaf.name)
        leaf.add_feature("tax", taxonomy[target_label])
        labels[taxonomy[target_label]] = 1
    groups = labels.keys()

    # compute the number of clades (weizhi: monophyletic group) per label in the ML tree, and their sizes
    ML_groups = defaultdict(list)  # the list is the size of each clade, len(list) is the number of clades for that label in the ML tree
    # ML_groups: the number of leaves in each monophyletic groups of the corresponding target_label (e.g. genus)
    for label in groups:
        node_num = 0
        for monophyletic_clade in ml_tree.get_monophyletic(values=[label], target_attr="tax"):  # get monophyletic groups for each target_label (e.g. genus)
            # print('node')
            # print(node)
            size_clade = 0  # get the number of leaf (size_clade) in the monophyletic group
            for leaf in monophyletic_clade:
                size_clade += 1
            ML_groups[label].append(size_clade)
            node_num += 1
        # print('The number of monophyletic clade (node_num) of %s (label):\t%s' % (label, node_num))

    # print()
    # print('Dict holds the number of monophyletic clades per taxon, and their sizes')
    # print('ML_groups:\t %s' % ML_groups)
    # print()

    summary = defaultdict(dict)
    clades_per_group = defaultdict(list)
    treeNum = -1
    for line in open(tree_in_bs):  # read in each bootstrap tree

        treeNum += 1
        tree = Tree(line.rstrip())
        for leaf in tree:
            tax = name_to_tax_info[leaf.name]  # this should set up taxonomy correctly...
            leaf.add_feature("tax", tax[target_label])  # this adds a feature called tax to the leaf, with the attribute of the phylum name
        for label in groups:
            clades_per_group[label].append(0.0)  # setup the clade counting for this particular tree
        tree.unroot()  # Weizhi: why is this

        # iterate over groups that are monophyletic for the taxon label of choice.
        # Choose the smallest sister branch for the comparison. (Assume root is within the larger sister clade (Weizhi:why?))
        for label in groups:

            # print('tree.get_monophyletic(values=[label], target_attr="tax")')
            # print(tree.get_monophyletic(values=[label], target_attr="tax"))
            # print('---------------------------------------------------------------------------------------------------v')
            # print('label: %s' % label)
            monophyletic_clade_index = 1
            for monophyletic_clade in tree.get_monophyletic(values=[label], target_attr="tax"):  # node: monophyletic clade
                clades_per_group[label][treeNum] += 1.0
                # print node.get_ascii()
                sister_clades = monophyletic_clade.get_sisters()

                # print('--------------------v')
                # print('monophyletic clade %s in %s (label)' % (monophyletic_clade_index, label))
                monophyletic_clade_index += 1
                #print(monophyletic_clade)
                # print(len(sisters))
                # for leaf in sisters[0]:
                #     print(leaf.name)
                # print(sisters)
                # print('sisters of current monophyletic clade')
                sister_index = 1
                for each_sister in sister_clades:
                    current_sister_leaf_list = []
                    for leaf in each_sister:
                        current_sister_leaf_list.append(leaf.name)
                    # print('sister %s has %s leaves: %s' % (sister_index, len(current_sister_leaf_list), ','.join([])))
                    sister_index += 1

                if monophyletic_clade.is_root():  # monophyletic clade is root
                    continue

                # Weizhi: bifurcation
                elif len(sister_clades) == 1:  # not at the trifurcation. Do something a bit hacky to find the bigger sister clade

                    taxa_in_sister = []
                    for leaf in sister_clades[0]:
                        taxa_in_sister.append(leaf.name)

                    size_sister = len(taxa_in_sister)

                    taxa_in_group = []
                    for leaf in monophyletic_clade:
                        taxa_in_group.append(leaf.name)

                    taxa_in_other_groups = []  # what does OG mean? (other groups?)
                    for leaf_name in all_tree_leaf_names:
                        if leaf_name in taxa_in_sister:
                            continue
                        elif leaf_name in taxa_in_group:
                            continue
                        else:
                            taxa_in_other_groups.append(leaf_name)
                    size_other_groups = len(taxa_in_other_groups)

                    sister_tax = {}  # taxa in the smaller groups (either the sister group or the OG)
                    if size_other_groups > size_sister:
                        sister_tax = summarize_taxonomy(taxa_in_sister, target_label, name_to_tax_info)
                    else:
                        sister_tax = summarize_taxonomy(taxa_in_other_groups, target_label, name_to_tax_info)

                    # print('size_sister: %s' % size_sister)
                    # print('size_other_groups: %s' % size_other_groups)
                    # print('sister_tax (not really, actually taxa in the smaller one (either the sister group or the OG))')
                    # print(sister_tax)

                    # store the tax info of the sister group
                    for element in sister_tax:
                        # print('element: %s' % element)
                        #print('summary[label]: %s' % summary[label])
                        if element in summary[label]:
                            summary[label][element] += sister_tax[element]
                            #print('summary (in): %s' % summary)
                        else:
                            summary[label][element] = sister_tax[element]
                            #print('summary (not in): %s' % summary)

                else:  # trifurcation in tree. Just treat the two sisters in the same way.

                    taxa_in_sisters_1 = []
                    for leaf in sister_clades[0]:
                        taxa_in_sisters_1.append(leaf.name)

                    taxa_in_sisters_2 = []
                    for leaf in sister_clades[1]:
                        taxa_in_sisters_2.append(leaf.name)

                    # get the size of two sisters
                    size_s1 = len(taxa_in_sisters_1)
                    size_s2 = len(taxa_in_sisters_2)

                    # print('size_s1: %s' % size_s1)
                    # print('size_s2: %s' % size_s2)

                    # get taxa in the smaller sister group
                    sister_tax = {}
                    if size_s1 > size_s2:
                        sister_tax = summarize_taxonomy(taxa_in_sisters_2, target_label, name_to_tax_info)
                    else:
                        sister_tax = summarize_taxonomy(taxa_in_sisters_1, target_label, name_to_tax_info)

                    # print('sister_tax (taxa in the smaller sister group)')
                    # print(sister_tax)

                    for element in sister_tax:
                        if element in summary[label]:
                            summary[label][element] += sister_tax[element]
                        else:
                            summary[label][element] = sister_tax[element]

            #     print('--------------------^')
            #     print()
            # print('---------------------------------------------------------------------------------------------------^')

    # now print out some kind of summary. For each label, the sorted list of sister taxa and their frequencies?
    outh = open(output_file, "w")
    for label in summary:
        num_groups = len(ML_groups[label])
        size_str = ''
        if num_groups == 1:
            size_str = ML_groups[label][0]
        else:
            size_str = ','.join(str(x) for x in (sorted(ML_groups[label], reverse=True)))

        avg_num_clades = float("{0:.4f}".format(numpy.mean(clades_per_group[label])))
        total_num_clades = numpy.sum(clades_per_group[label])
        sorted_sisters = sorted(summary[label].items(), key=itemgetter(1), reverse=True)

        # if label == 'g__TA-20':
        #     print('ML_groups[label]:\t%s'               % ML_groups[label])
        #     print('clades_per_group[label]:\t%s'        % clades_per_group[label])
        #     print('avg_num_clades (mean of list):\t%s'  % numpy.mean(clades_per_group[label]))
        #     print('total_num_clades (sum of list):\t%s' % total_num_clades)
        #     print('summary[label]:\t%s'                 % summary[label])
        #     print('sorted_sisters\t%s'                  % sorted_sisters)

        for tup in sorted_sisters:
            double_normalize = float(tup[1]) / float(total_num_clades)  # normalize the frequencies by the total number of clades, to account for different bootstrap numbers/MCMC sample numbers
            double_normalize = float("{0:.4f}".format(double_normalize))

            str_to_write = '%s\t%s\t%s\t%s\t%s\t%s' % (label, tup[0], float("{0:.4f}".format(tup[1])), avg_num_clades, double_normalize, size_str)
            outh.write(str_to_write + '\n')
            # if label == 'g__TA-20':
            #     print(str_to_write)
    outh.close()


def count_sister_taxa_worker(arg_list):

    mag_rename_dict                 = arg_list[0]
    mag_cluster_dict                = arg_list[1]
    sponge_archaeal_MAG_tax_dict    = arg_list[2]
    gtdb_ar_gnm_tax_dict            = arg_list[3]
    tree_ml                         = arg_list[4]
    ufboot_file                     = arg_list[5]
    target_label                    = arg_list[6]
    tree_ml_renamed                 = arg_list[7]
    ufboot_file_renamed             = arg_list[8]
    count_sister_taxa_op_txt        = arg_list[9]
    gene_id                         = arg_list[10]
    renamed_gnm_to_cluster_dir      = arg_list[11]

    # rename ml tree
    tree_ml_renamed_handle = open(tree_ml_renamed, 'w')
    current_tree_rename_dict = get_rename_dict(tree_ml, mag_rename_dict, mag_cluster_dict, sponge_archaeal_MAG_tax_dict, gtdb_ar_gnm_tax_dict)
    tree_ml_str_renamed = rename_tree(tree_ml, current_tree_rename_dict)
    tree_ml_renamed_handle.write(tree_ml_str_renamed + '\n')
    tree_ml_renamed_handle.close()

    current_renamed_gnm_to_cluster_txt = '%s/%s.txt' % (renamed_gnm_to_cluster_dir, gene_id)
    current_renamed_gnm_to_cluster_txt_handle = open(current_renamed_gnm_to_cluster_txt, 'w')
    for each_leaf in current_tree_rename_dict:
        renamed_leaf = current_tree_rename_dict[each_leaf]
        cluster_id = renamed_leaf.split('|')[0]
        current_renamed_gnm_to_cluster_txt_handle.write('%s\t%s\n' % (renamed_leaf, cluster_id))
    current_renamed_gnm_to_cluster_txt_handle.close()

    # rename ufboot trees
    ufboot_file_renamed_handle = open(ufboot_file_renamed, 'w')
    for each_tree in open(ufboot_file):
        tree_str = each_tree.strip()
        current_tree_rename_dict = get_rename_dict(tree_str, mag_rename_dict, mag_cluster_dict, sponge_archaeal_MAG_tax_dict, gtdb_ar_gnm_tax_dict)
        tree_str_renamed = rename_tree(tree_str, current_tree_rename_dict)
        ufboot_file_renamed_handle.write(tree_str_renamed + '\n')
    ufboot_file_renamed_handle.close()

    # count_sister_taxa
    count_sister_taxa(target_label, tree_ml_renamed, ufboot_file_renamed, count_sister_taxa_op_txt)


def run_count_sister_taxa(genome_metadata_ar53_r207_Mac, sponge_MAG_GTDB_archaea, hog_id_txt, contree_and_ufboot_dir, archaeal_mags_renamed_for_prokka_txt, gnm_cluster_txt, target_label, num_threads, output_dir):

    # define file name
    renamed_gnm_to_cluster_dir              = '%s/renamed_genome_to_cluster'            % output_dir
    renamed_gnm_to_cluster_tmp_txt          = '%s/renamed_genome_to_cluster_tmp.txt'    % output_dir
    renamed_gnm_to_cluster_txt              = '%s/renamed_genome_to_cluster.txt'        % output_dir
    renamed_gnm_to_cluster_iTOL_txt         = '%s/renamed_genome_to_cluster_iTOL.txt'   % output_dir
    renamed_contree_dir                     = '%s/renamed_contree'                      % output_dir
    renamed_ufboot_dir                      = '%s/renamed_ufboot'                       % output_dir
    count_sister_taxa_op_dir                = '%s/count_sister_taxa_op'                 % output_dir

    os.mkdir(output_dir)
    os.mkdir(renamed_contree_dir)
    os.mkdir(renamed_ufboot_dir)
    os.mkdir(count_sister_taxa_op_dir)
    os.mkdir(renamed_gnm_to_cluster_dir)

    _, _, gtdb_ar_gnm_tax_dict, _ = gtdb_gnm_metadata_parser(genome_metadata_ar53_r207_Mac)

    sponge_archaeal_MAG_tax_dict = {}
    for each in open(sponge_MAG_GTDB_archaea):
        if not each.startswith('user_genome'):
            each_split = each.strip().split('\t')
            sponge_archaeal_MAG_tax_dict[each_split[0]] = each_split[1]

    hog_list = []
    for each_hog in open(hog_id_txt):
        hog_list.append(each_hog.strip())

    mag_cluster_dict = {}
    for each_gnm in open(gnm_cluster_txt):
        each_gnm_split = each_gnm.strip().split('\t')
        mag_cluster_dict[each_gnm_split[1]] = each_gnm_split[0]

    mag_rename_dict = {}
    for each_mag in open(archaeal_mags_renamed_for_prokka_txt):
        each_mag_split = each_mag.strip().split('\t')
        before_rename = each_mag_split[0]
        after_rename = each_mag_split[1]
        mag_rename_dict[after_rename] = before_rename

    argument_lol = []
    for og_id in hog_list:

        # define file name
        tree_ml                  = '%s/%s_iqtree.contree'               % (contree_and_ufboot_dir, og_id)
        ufboot_file              = '%s/%s_iqtree.ufboot'                % (contree_and_ufboot_dir, og_id)
        tree_ml_renamed          = '%s/%s_iqtree_renamed.contree'       % (renamed_contree_dir, og_id)
        ufboot_file_renamed      = '%s/%s_iqtree_renamed.ufboot'        % (renamed_ufboot_dir, og_id)
        count_sister_taxa_op_txt = '%s/%s_iqtree_count_sister_taxa.txt' % (count_sister_taxa_op_dir, og_id)

        current_arg_list = [mag_rename_dict, mag_cluster_dict, sponge_archaeal_MAG_tax_dict, gtdb_ar_gnm_tax_dict, tree_ml, ufboot_file, target_label, tree_ml_renamed, ufboot_file_renamed, count_sister_taxa_op_txt, og_id, renamed_gnm_to_cluster_dir]
        argument_lol.append(current_arg_list)

    # run with multiprocessing
    pool = mp.Pool(processes=num_threads)
    pool.map(count_sister_taxa_worker, argument_lol)
    pool.close()
    pool.join()

    # combine renamed_gnm_to_cluster files
    os.system('cat %s/*.txt > %s' % (renamed_gnm_to_cluster_dir, renamed_gnm_to_cluster_tmp_txt))
    os.system('cat %s | sort | uniq > %s' % (renamed_gnm_to_cluster_tmp_txt, renamed_gnm_to_cluster_txt))
    BioSAK_iTOL_cmd = 'BioSAK iTOL -ColorRange -lg %s -lt Cluster -out %s' % (renamed_gnm_to_cluster_txt, renamed_gnm_to_cluster_iTOL_txt)
    os.system(BioSAK_iTOL_cmd)


####################################################### Test 2023-07-11 ########################################################

genome_metadata_ar53_r207_Mac           = '/Users/songweizhi/DB/GTDB_r207/ar53_metadata_r207.tsv'
sponge_MAG_GTDB_archaea                 = '/Users/songweizhi/Documents/Research/Sponge_Hologenome/0_backup/0_metadata_mac/Sponge_MAGs_1677.ar53.summary.tsv'
hog_id_txt                              = '/Users/songweizhi/Documents/Research/Sponge_Hologenome/0_backup/5_Archaeal_tree_50_5_Markers_by_split_wd/HOG_id.txt'
contree_and_ufboot_dir                  = '/Users/songweizhi/Documents/Research/Sponge_Hologenome/0_backup/5_Archaeal_tree_50_5_Markers_by_split_wd/contree_and_ufboot_files'
archaeal_mags_renamed_for_prokka_txt    = '/Users/songweizhi/Documents/Research/Sponge_Hologenome/0_backup/0_metadata_mac/Archaeal_mags_renamed_for_prokka.txt'
gnm_cluster_txt                         = '/Users/songweizhi/Documents/Research/Sponge_Hologenome/0_backup/5_Archaeal_tree_50_5_Markers_by_split_wd/genome_clusters_v1.txt'
target_label                            = 'cluster'
num_threads                             = 10
output_dir                              = '/Users/songweizhi/Desktop/count_sister_taxa_op_Test_2023_07_11'

run_count_sister_taxa(genome_metadata_ar53_r207_Mac, sponge_MAG_GTDB_archaea, hog_id_txt, contree_and_ufboot_dir, archaeal_mags_renamed_for_prokka_txt, gnm_cluster_txt, target_label, num_threads, output_dir)
