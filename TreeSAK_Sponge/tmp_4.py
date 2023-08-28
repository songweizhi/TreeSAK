
op_txt = '/Users/songweizhi/Desktop/Sponge_2023_08_25/s1_get_genome/Cdb_dRep97_rename_iTOL.txt'


representative_gnm_set = set()
for each_gnm in open('/Users/songweizhi/Desktop/Cdb_dRep97_representative_gnms.txt'):
    representative_gnm_set.add(each_gnm.strip())

gnm_to_cluster_dict = dict()
for each in open('/Users/songweizhi/Desktop/Cdb_dRep97.txt'):
    each_split = each.strip().split('\t')
    gnm_to_cluster_dict[each_split[0]] = each_split[1]

gnm_to_cpl_dict = dict()
gnm_to_ctm_dict = dict()
for each in open('/Users/songweizhi/Desktop/Sponge_2023_08_25/s1_get_genome/Nitrosopumilaceae_50_5_GTDB_NCBI_367_Sponge_quality.txt'):
    each_split = each.strip().split(',')
    gnm_to_cpl_dict[each_split[0]] = each_split[1]
    gnm_to_ctm_dict[each_split[0]] = each_split[2]

gnm_host_dict = dict()
gnm_host_genus_dict = dict()
for each_gnm in open('/Users/songweizhi/Documents/Research/Sponge_Hologenome/1_metadata/0_metadata_final.txt'):
    each_gnm_split = each_gnm.strip().split('\t')
    gnm_id = each_gnm_split[0]
    host_taxon = each_gnm_split[9]
    host_genus = 'na'
    if host_taxon not in ['nonsponge', 'na']:
        host_genus = host_taxon.split(';')[-1]
    gnm_host_dict[gnm_id] = host_taxon
    gnm_host_genus_dict[gnm_id] = host_genus

op_txt_handle = open(op_txt, 'w')
op_txt_handle.write('LABELS\nSEPARATOR TAB\nDATA\n')
for each_gnm in gnm_to_cluster_dict:
    gnm_id = each_gnm.replace('.sponge', '').replace('.gtdb', '').replace('.ncbi', '')
    if each_gnm in representative_gnm_set:
        op_txt_handle.write('%s\t%s__%s__%s__%s__%s__representative\n' % (each_gnm, each_gnm, gnm_to_cluster_dict[each_gnm], gnm_to_cpl_dict[each_gnm], gnm_to_ctm_dict[each_gnm], gnm_host_genus_dict[gnm_id]))
    else:
        op_txt_handle.write('%s\t%s__%s__%s__%s__%s\n' % (each_gnm, each_gnm, gnm_to_cluster_dict[each_gnm], gnm_to_cpl_dict[each_gnm], gnm_to_ctm_dict[each_gnm], gnm_host_genus_dict[gnm_id]))
op_txt_handle.close()

