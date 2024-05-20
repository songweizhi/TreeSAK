from __future__ import print_function
import os
import glob
import numpy
import argparse
import subprocess
from ete3 import Tree
import multiprocessing as mp
from operator import itemgetter
from collections import defaultdict


SplitScore2_usage = '''
======================== SplitScore2 example commands ========================

TreeSAK SplitScore2 -i step1_op_dir -g gnm_cluster.tsv -k gnm_taxon.txt -f -t 10 -o step_2_op_dir

# format of gnm_cluster.tsv (tab separated)
GCA_013330055.1 c01_UBA8516
GCA_023251795.1 c01_UBA8516
GCA_023251295.1 c01_UBA8516
GCA_005877305.1 c02_TA-20
GCA_013287585.1 c02_TA-20

# gnm_taxon.txt: GTDB format

# install R packages
install.packages("optparse")
install.packages("plyr")
install.packages("dbplyr")
install.packages("dplyr")
install.packages("tidyr")
install.packages("ggplot2")
install.packages("data.table")
install.packages("RColorBrewer")
install.packages("gplots")
install.packages("ape")

=============================================================================
'''


def sep_path_basename_ext(file_in):
    f_path, file_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'
    f_base, f_ext = os.path.splitext(file_name)
    return f_path, f_base, f_ext


def check_executables(program_list):

    not_detected_programs = []
    for needed_program in program_list:

        if subprocess.call(['which', needed_program], stdout=open(os.devnull, 'wb')) != 0:
            not_detected_programs.append(needed_program)

    if not_detected_programs != []:
        print('%s not detected, program exited!' % ','.join(not_detected_programs))
        exit()


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


def get_rename_dict(tree_str_in, mag_cluster_dict, gtdb_gnm_tax_dict):

    # rename dict: {'old_name':'new_name'}

    leaf_rename_dict = {}
    for leaf in Tree(tree_str_in, format=1):

        leaf_name_gnm = '_'.join(leaf.name.split('_')[:-1])
        leaf_cluster = mag_cluster_dict.get(leaf_name_gnm, 'cluster_0')

        # get mag_taxon_str
        gnm_taxon_str = 'NA'
        if leaf_name_gnm in gtdb_gnm_tax_dict:
            gnm_taxon_str = gtdb_gnm_tax_dict[leaf_name_gnm]

        # get mag_taxon_str (GCA GCF things)
        if gnm_taxon_str == 'NA':
            mag_id_no_ext_no_source_GCF = leaf_name_gnm.replace('GCA', 'GCF')
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


def parse_taxonomy(taxon_name):  # given a taxon name, try to return whatever taxonomic info is available as a list starting with the highest level classification and going lower (or a map?)

    name_elements = taxon_name.split('|')
    if (len(name_elements) < 8) or (len(name_elements) > 9):
        print("Nonstandard!")
        quit()

    name_map = dict()
    name_map['cluster'] = name_elements[0]
    name_map['domain']  = name_elements[1]
    name_map['phylum']  = name_elements[2]
    name_map['class']   = name_elements[3]
    name_map['order']   = name_elements[4]
    name_map['family']  = name_elements[5]
    name_map['genus']   = name_elements[6]
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
            size_clade = 0  # get the number of leaf (size_clade) in the monophyletic group
            for leaf in monophyletic_clade:
                size_clade += 1
            ML_groups[label].append(size_clade)
            node_num += 1

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
            monophyletic_clade_index = 1
            for monophyletic_clade in tree.get_monophyletic(values=[label], target_attr="tax"):  # node: monophyletic clade
                clades_per_group[label][treeNum] += 1.0
                sister_clades = monophyletic_clade.get_sisters()
                monophyletic_clade_index += 1
                sister_index = 1
                for each_sister in sister_clades:
                    current_sister_leaf_list = []
                    for leaf in each_sister:
                        current_sister_leaf_list.append(leaf.name)
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

                    # store the tax info of the sister group
                    for element in sister_tax:
                        if element in summary[label]:
                            summary[label][element] += sister_tax[element]
                        else:
                            summary[label][element] = sister_tax[element]

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

                    # get taxa in the smaller sister group
                    sister_tax = {}
                    if size_s1 > size_s2:
                        sister_tax = summarize_taxonomy(taxa_in_sisters_2, target_label, name_to_tax_info)
                    else:
                        sister_tax = summarize_taxonomy(taxa_in_sisters_1, target_label, name_to_tax_info)

                    for element in sister_tax:
                        if element in summary[label]:
                            summary[label][element] += sister_tax[element]
                        else:
                            summary[label][element] = sister_tax[element]

    # now print out some kind of summary. For each label, the sorted list of sister taxa and their frequencies?
    outh = open(output_file, "w")
    for label in summary:
        num_groups = len(ML_groups[label])
        size_str = ''
        if num_groups == 1:
            size_str = ML_groups[label][0]
        else:
            size_str = ','.join(str(x) for x in (sorted(ML_groups[label], reverse=True)))

        avg_num_clades   = float("{0:.4f}".format(numpy.mean(clades_per_group[label])))
        total_num_clades = numpy.sum(clades_per_group[label])
        sorted_sisters   = sorted(summary[label].items(), key=itemgetter(1), reverse=True)

        for tup in sorted_sisters:
            double_normalize = float(tup[1]) / float(total_num_clades)  # normalize the frequencies by the total number of clades, to account for different bootstrap numbers/MCMC sample numbers
            double_normalize = float("{0:.4f}".format(double_normalize))
            str_to_write     = '%s\t%s\t%s\t%s\t%s\t%s' % (label, tup[0], float("{0:.4f}".format(tup[1])), avg_num_clades, double_normalize, size_str)
            outh.write(str_to_write + '\n')
    outh.close()


def count_sister_taxa_worker(arg_list):

    mag_cluster_dict                = arg_list[0]
    gnm_tax_dict                    = arg_list[1]
    tree_ml                         = arg_list[2]
    ufboot_file                     = arg_list[3]
    target_label                    = arg_list[4]
    tree_ml_renamed                 = arg_list[5]
    ufboot_file_renamed             = arg_list[6]
    count_sister_taxa_op_txt        = arg_list[7]
    gene_id                         = arg_list[8]
    renamed_gnm_to_cluster_dir      = arg_list[9]

    # rename ml tree
    tree_ml_renamed_handle = open(tree_ml_renamed, 'w')
    current_tree_rename_dict = get_rename_dict(tree_ml, mag_cluster_dict, gnm_tax_dict)
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
        current_tree_rename_dict = get_rename_dict(tree_str, mag_cluster_dict, gnm_tax_dict)
        tree_str_renamed = rename_tree(tree_str, current_tree_rename_dict)
        ufboot_file_renamed_handle.write(tree_str_renamed + '\n')
    ufboot_file_renamed_handle.close()

    # count_sister_taxa
    count_sister_taxa(target_label, tree_ml_renamed, ufboot_file_renamed, count_sister_taxa_op_txt)


def run_count_sister_taxa(gtdb_classification_txt, hog_list, contree_dir, ufboot_dir, gnm_cluster_txt, target_label, num_threads, output_dir, force_overwrite):

    # define file name
    renamed_gnm_to_cluster_dir      = '%s/renamed_genome_to_cluster'            % output_dir
    renamed_gnm_to_cluster_tmp_txt  = '%s/renamed_genome_to_cluster_tmp.txt'    % output_dir
    renamed_gnm_to_cluster_txt      = '%s/renamed_genome_to_cluster.txt'        % output_dir
    renamed_gnm_to_cluster_iTOL_txt = '%s/renamed_genome_to_cluster_iTOL.txt'   % output_dir
    renamed_contree_dir             = '%s/renamed_contree'                      % output_dir
    renamed_ufboot_dir              = '%s/renamed_ufboot'                       % output_dir
    count_sister_taxa_op_dir        = '%s/count_sister_taxa_op'                 % output_dir

    if os.path.isdir(output_dir) is True:
        if force_overwrite is True:
            os.system('rm -r %s' % output_dir)
        else:
            print('%s exist, program exited!' % output_dir)
            exit()

    os.mkdir(output_dir)
    os.mkdir(renamed_contree_dir)
    os.mkdir(renamed_ufboot_dir)
    os.mkdir(count_sister_taxa_op_dir)
    os.mkdir(renamed_gnm_to_cluster_dir)

    ####################################################################################################################

    gnm_tax_dict = {}
    for each in open(gtdb_classification_txt):
        if not each.startswith('user_genome'):
            each_split = each.strip().split('\t')
            gnm_tax_dict[each_split[0]] = each_split[1]

    ####################################################################################################################

    mag_cluster_dict = {}
    for each_gnm in open(gnm_cluster_txt):
        each_gnm_split = each_gnm.strip().split('\t')
        mag_cluster_dict[each_gnm_split[0]] = each_gnm_split[1]

    argument_lol = []
    for og_id in hog_list:

        # define file name
        tree_ml                  = '%s/%s.contree'               % (contree_dir, og_id)
        ufboot_file              = '%s/%s.ufboot'                % (ufboot_dir, og_id)
        tree_ml_renamed          = '%s/%s_renamed.contree'       % (renamed_contree_dir, og_id)
        ufboot_file_renamed      = '%s/%s_renamed.ufboot'        % (renamed_ufboot_dir, og_id)
        count_sister_taxa_op_txt = '%s/%s_count_sister_taxa.txt' % (count_sister_taxa_op_dir, og_id)

        if os.path.isfile(tree_ml) is False:
            print('%s not found!' % tree_ml)
            exit()

        current_arg_list = [mag_cluster_dict, gnm_tax_dict, tree_ml, ufboot_file, target_label, tree_ml_renamed, ufboot_file_renamed, count_sister_taxa_op_txt, og_id, renamed_gnm_to_cluster_dir]
        argument_lol.append(current_arg_list)

    # run with multiprocessing
    pool = mp.Pool(processes=num_threads)
    pool.map(count_sister_taxa_worker, argument_lol)
    pool.close()
    pool.join()

    # combine renamed_gnm_to_cluster files
    os.system('cat %s/*.txt > %s'                                           % (renamed_gnm_to_cluster_dir, renamed_gnm_to_cluster_tmp_txt))
    os.system('cat %s | sort | uniq > %s'                                   % (renamed_gnm_to_cluster_tmp_txt, renamed_gnm_to_cluster_txt))
    BioSAK_iTOL_cmd = 'BioSAK iTOL -ColorRange -lg %s -lt Cluster -o %s'    % (renamed_gnm_to_cluster_txt, renamed_gnm_to_cluster_iTOL_txt)
    os.system(BioSAK_iTOL_cmd)


def get_taxa_count_stats(step_1_op_dir, hog_list_sorted, get_taxa_count_stats_wd, force_overwrite, TaxaCountStats_Rscript):

    # define input files to R script
    combined_contree_file           = '%s/combined.contree'                     % get_taxa_count_stats_wd
    genes_to_remove_txt             = '%s/Genes_to_remove.txt'                  % get_taxa_count_stats_wd
    list_of_trees_txt               = '%s/List_of_trees.txt'                    % get_taxa_count_stats_wd
    mapping_txt                     = '%s/mapping.txt'                          % get_taxa_count_stats_wd
    marker_list_txt                 = '%s/MarkerList.txt'                       % get_taxa_count_stats_wd
    combined_count_sister_taxa_op   = '%s/combined_count_sister_taxa_op.txt'    % get_taxa_count_stats_wd
    TaxaCountStats_op               = '%s/TaxaCountStats_output.txt'            % get_taxa_count_stats_wd

    if os.path.isdir(get_taxa_count_stats_wd) is True:
        if force_overwrite is True:
            os.system('rm -r %s' % get_taxa_count_stats_wd)
        else:
            print('%s exist, program exited!' % get_taxa_count_stats_wd)
            exit()
    os.mkdir(get_taxa_count_stats_wd)

    cluster_to_domain_dict = {}
    marker_list_txt_handle = open(marker_list_txt, 'w')
    marker_list_txt_handle.write('MarkerID\n')
    list_of_trees_txt_handle = open(list_of_trees_txt, 'w')
    combined_contree_file_handle = open(combined_contree_file, 'w')
    combined_count_sister_taxa_op_handle = open(combined_count_sister_taxa_op, 'w')
    combined_count_sister_taxa_op_handle.write('MarkerID\tGroup_of_interest\tSister_taxa\tNormalized_sum_of_occurances\tsplits\tNormalized2_sum_of_occurances\tClusters\n')
    for each_hog in hog_list_sorted:

        # write out to combined_count_sister_taxa_op
        pwd_count_sister_taxa_op_txt = '%s/count_sister_taxa_op/%s_count_sister_taxa.txt' % (step_1_op_dir, each_hog)
        with open(pwd_count_sister_taxa_op_txt) as count_sister_taxa_op_txt_opened:
            for each_line in count_sister_taxa_op_txt_opened:
                combined_count_sister_taxa_op_handle.write('%s\t%s' % (each_hog, each_line))

        # write out to combined_contree_file
        pwd_renamed_contree_file = '%s/renamed_contree/%s_renamed.contree' % (step_1_op_dir, each_hog)
        with open(pwd_renamed_contree_file, 'r') as pwd_renamed_contree_file_opened:
            combined_contree_file_handle.write(pwd_renamed_contree_file_opened.readline())

        # add to cluster_to_domain_dict
        t_in = Tree(pwd_renamed_contree_file, format=1)
        for leaf in t_in:
            leaf_name_split = leaf.name.split('|')
            cluster_to_domain_dict[leaf_name_split[0]] = leaf_name_split[1]

        # write out to marker_list_txt
        marker_list_txt_handle.write(each_hog + '\n')
        list_of_trees_txt_handle.write(each_hog + '\n')

    marker_list_txt_handle.close()
    list_of_trees_txt_handle.close()
    combined_contree_file_handle.close()
    combined_count_sister_taxa_op_handle.close()

    # prepare mapping_txt
    mapping_txt_handle = open(mapping_txt, 'w')
    mapping_txt_handle.write('Cluster\tDomain\n')
    for each_cluster in cluster_to_domain_dict:
        mapping_txt_handle.write('%s\t%s\n' % (each_cluster, cluster_to_domain_dict[each_cluster]))
    mapping_txt_handle.close()

    # prepare genes_to_remove_txt
    genes_to_remove_txt_handle = open(genes_to_remove_txt, 'w')
    genes_to_remove_txt_handle.write('MarkerID\n')
    genes_to_remove_txt_handle.close()

    # run TaxaCountStats.R
    get_TaxaCountStats_cmd = 'Rscript %s -t %s -l %s -g %s -x %s -s %s -r %s -o %s > /dev/null' % (TaxaCountStats_Rscript, combined_contree_file, list_of_trees_txt, mapping_txt, marker_list_txt, combined_count_sister_taxa_op, genes_to_remove_txt, TaxaCountStats_op)
    print('Running: ' + get_TaxaCountStats_cmd)
    os.system(get_TaxaCountStats_cmd)


def group_marker(taxa_counts_tats_op_txt, marker_seq_dir, op_dir):

    # define file name
    marker_set_top_25_txt       = '%s/top25.txt'    % (op_dir)
    marker_set_top_50_txt       = '%s/top50.txt'    % (op_dir)
    marker_set_top_75_txt       = '%s/top75.txt'    % (op_dir)
    marker_set_top_100_txt      = '%s/top100.txt'   % (op_dir)
    marker_set_top_25_seq_dir   = '%s/top25'        % (op_dir)
    marker_set_top_50_seq_dir   = '%s/top50'        % (op_dir)
    marker_set_top_75_seq_dir   = '%s/top75'        % (op_dir)
    marker_set_top_100_seq_dir  = '%s/top100'       % (op_dir)

    os.system('mkdir %s' % marker_set_top_25_seq_dir)
    os.system('mkdir %s' % marker_set_top_50_seq_dir)
    os.system('mkdir %s' % marker_set_top_75_seq_dir)
    os.system('mkdir %s' % marker_set_top_100_seq_dir)

    marker_set_top_25 = set()
    marker_set_top_50 = set()
    marker_set_top_75 = set()
    marker_set_top_100 = set()
    header_index_dict = {}
    for each_marker in open(taxa_counts_tats_op_txt):
        each_marker_split = each_marker.replace('\n', '').split('\t')
        if each_marker.startswith('MarkerID\t'):
            header_index_dict = {k: v for v, k in enumerate(each_marker_split)}
        else:
            marker_id = each_marker_split[header_index_dict['MarkerID']]
            best_25perc = each_marker_split[header_index_dict['best_25perc']]
            best_50perc = each_marker_split[header_index_dict['best_50perc']]
            worst_25perc = each_marker_split[header_index_dict['worst_25perc']]
            if best_25perc != '':
                marker_set_top_25.add(marker_id)
                os.system('cp %s/%s.fa %s/' % (marker_seq_dir, marker_id, marker_set_top_25_seq_dir))
            if best_50perc != '':
                marker_set_top_50.add(marker_id)
                os.system('cp %s/%s.fa %s/' % (marker_seq_dir, marker_id, marker_set_top_50_seq_dir))
            if worst_25perc == '':
                marker_set_top_75.add(marker_id)
                os.system('cp %s/%s.fa %s/' % (marker_seq_dir, marker_id, marker_set_top_75_seq_dir))
            marker_set_top_100.add(marker_id)
            os.system('cp %s/%s.fa %s/' % (marker_seq_dir, marker_id, marker_set_top_100_seq_dir))

    with open(marker_set_top_25_txt, 'w') as marker_set_top_25_txt_handle:
        marker_set_top_25_txt_handle.write('\n'.join(sorted([i for i in marker_set_top_25])))
    with open(marker_set_top_50_txt, 'w') as marker_set_top_50_txt_handle:
        marker_set_top_50_txt_handle.write('\n'.join(sorted([i for i in marker_set_top_50])))
    with open(marker_set_top_75_txt, 'w') as marker_set_top_75_txt_handle:
        marker_set_top_75_txt_handle.write('\n'.join(sorted([i for i in marker_set_top_75])))
    with open(marker_set_top_100_txt, 'w') as marker_set_top_100_txt_handle:
        marker_set_top_100_txt_handle.write('\n'.join(sorted([i for i in marker_set_top_100])))


def SplitScore2(args):

    step_1_op_dir               = args['i']
    gnm_group_txt               = args['g']
    gtdb_classification_txt     = args['k']
    force_overwrite             = args['f']
    num_of_threads              = args['t']
    step_2_op_dir               = args['o']
    target_label                = 'cluster'

    check_executables(['Rscript'])


    current_file_path           = '/'.join(os.path.realpath(__file__).split('/')[:-1])
    TaxaCountStats_Rscript      = '%s/TaxaCountStats.R'                                     % current_file_path
    qualified_og_seq_dir        = '%s/qualified_OGs'                                        % step_1_op_dir
    contree_file_re             = '%s/*.contree'                                            % qualified_og_seq_dir
    ufboot_file_re              = '%s/*.ufboot'                                             % qualified_og_seq_dir
    count_sister_taxa_op_dir    = '%s/count_sister_taxa_wd'                                 % step_2_op_dir
    get_taxa_count_stats_op_dir = '%s/get_taxa_count_stats_wd'                              % step_2_op_dir
    TaxaCountStats_output_txt   = '%s/get_taxa_count_stats_wd/TaxaCountStats_output.txt'    % step_2_op_dir

    contree_file_set_base = set()
    for each_contree_file in glob.glob(contree_file_re):
        _, f_base, _ = sep_path_basename_ext(each_contree_file)
        contree_file_set_base.add(f_base)

    ufboot_file_set_base = set()
    for each_ufboot_file in glob.glob(ufboot_file_re):
        _, f_base, _ = sep_path_basename_ext(each_ufboot_file)
        ufboot_file_set_base.add(f_base)

    contree_ufboot_shared = set(contree_file_set_base).intersection(ufboot_file_set_base)
    contree_ufboot_shared_sorted = sorted([i for i in contree_ufboot_shared])

    # create output folder
    if os.path.isdir(step_2_op_dir) is True:
        if force_overwrite is True:
            os.system('rm -r %s' % step_2_op_dir)
        else:
            print('%s exist, program exited!' % step_2_op_dir)
            exit()
    os.mkdir(step_2_op_dir)

    print('Counting sister taxa with %s cores' % num_of_threads)
    run_count_sister_taxa(gtdb_classification_txt, contree_ufboot_shared_sorted, qualified_og_seq_dir, qualified_og_seq_dir, gnm_group_txt, target_label, num_of_threads, count_sister_taxa_op_dir, force_overwrite)

    print('Summarising sister taxa')
    get_taxa_count_stats(count_sister_taxa_op_dir, contree_ufboot_shared_sorted, get_taxa_count_stats_op_dir, force_overwrite, TaxaCountStats_Rscript)

    print('Exporting markers by split score')
    group_marker(TaxaCountStats_output_txt, qualified_og_seq_dir, step_2_op_dir)

    print('Done!')


if __name__ == '__main__':

    SplitScore2_parser = argparse.ArgumentParser()
    SplitScore2_parser.add_argument('-i', required=True,                        help='output dir from SplitScore1')
    SplitScore2_parser.add_argument('-g', required=True,                        help='genome group')
    SplitScore2_parser.add_argument('-k', required=True,                        help='genome taxon, GTDB format')
    SplitScore2_parser.add_argument('-f', required=False, action="store_true",  help='force overwrite')
    SplitScore2_parser.add_argument('-t', required=False, type=int, default=1,  help='num of threads, default: 1')
    SplitScore2_parser.add_argument('-o', required=True,                        help='output directory')
    args = vars(SplitScore2_parser.parse_args())
    SplitScore2(args)
