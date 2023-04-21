import os

# file in
gnm_id_txt                          = '/Users/songweizhi/Desktop/DateArTree/03_decorate_tree/combined_gnm_id.txt'
gnm_taxon_txt                       = '/Users/songweizhi/Desktop/DateArTree/01_genome_selection/0_genome_metadata/gtdbtk.ar53.summary.tsv'
root_order_txt                      = '/Users/songweizhi/Desktop/DateArTree/03_decorate_tree/root_orders.txt'

# file out
itol_label_out                      = '/Users/songweizhi/Desktop/DateArTree/03_decorate_tree/iTOL_gnm_label.txt'
itol_color_range_root_order         = '/Users/songweizhi/Desktop/DateArTree/03_decorate_tree/iTOL_gnm_ColorRange_root_order.txt'
itol_color_range_root_order_iTOL    = '/Users/songweizhi/Desktop/DateArTree/03_decorate_tree/iTOL_gnm_ColorRange_root_order_iTOL.txt'

# read in gnm_taxon_txt
max_tax_str_len = 0
gnm_taxon_dict = dict()
for each_gnm in open(gnm_taxon_txt):
    each_gnm_split = each_gnm.strip().split('\t')
    gnm_id = each_gnm_split[0]
    acc_id = gnm_id
    if '.gtdb' in gnm_id:
        acc_id = gnm_id.split('.gtdb')[0]
    elif '.1_' in gnm_id:
        acc_id = gnm_id.split('.1_')[0] + '.1'
    elif '.2_' in gnm_id:
        acc_id = gnm_id.split('.2_')[0] + '.2'
    elif '.3_' in gnm_id:
        acc_id = gnm_id.split('.3_')[0] + '.3'
    else:
        pass
        #print(gnm_id)

    gnm_taxon = each_gnm_split[1]
    gnm_taxon_split = gnm_taxon.split(';')
    gnm_taxon_pco = ';'.join(gnm_taxon_split[1:4])
    if len(gnm_taxon_pco) > max_tax_str_len:
        max_tax_str_len = len(gnm_taxon_pco)
    gnm_taxon_dict[acc_id] = gnm_taxon_pco


# write out
itol_label_out_handle = open(itol_label_out, 'w')
itol_label_out_handle.write('LABELS\nSEPARATOR TAB\nDATA\n')
for each_gnm in open(gnm_id_txt):
    gnm_id = each_gnm.strip()
    gnm_taxon_str = gnm_taxon_dict.get(gnm_id, 'NA')
    itol_label_out_handle.write('%s\t%s__%s\n' % (gnm_id, gnm_id, gnm_taxon_str))
itol_label_out_handle.close()


root_order_list = []
for each_order in open(root_order_txt):
    root_order_list.append(each_order.strip())

itol_color_range_root_order_handle = open(itol_color_range_root_order, 'w')
for each_gnm in open(gnm_id_txt):
    gnm_id = each_gnm.strip()
    gnm_taxon_str = gnm_taxon_dict.get(gnm_id, 'NA')
    gnm_taxon_str_split = gnm_taxon_str.split(';')
    tax_order = ''
    for each_rank in gnm_taxon_str_split:
        if each_rank.startswith('o__'):
            tax_order = each_rank
    if tax_order in root_order_list:
        itol_color_range_root_order_handle.write('%s\troot_order\n' % (gnm_id))
    # else:
    #     itol_color_range_root_order_handle.write('%s\tnonroot_order\n' % (gnm_id))
itol_color_range_root_order_handle.close()

itol_cmd = 'BioSAK iTOL -ColorRange -lg %s -lt RootOrder -out %s' % (itol_color_range_root_order, itol_color_range_root_order_iTOL)
os.system(itol_cmd)

