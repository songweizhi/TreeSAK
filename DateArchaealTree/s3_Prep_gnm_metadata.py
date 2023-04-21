
# file in
gnm_id_txt      = '/Users/songweizhi/Desktop/DateArTree/01_genome_selection/0_genome_metadata/gnm_id.txt'
gnm_quality_txt = '/Users/songweizhi/Desktop/DateArTree/01_genome_selection/0_genome_metadata/d__Archaea_o_rs_gnms_Betts38_Spang81_Williams39_Wu2_quality_reformatted.txt'
gnm_size_txt    = '/Users/songweizhi/Desktop/DateArTree/01_genome_selection/0_genome_metadata/d__Archaea_o_rs_gnms_Betts38_Spang81_Williams39_Wu2_size.txt'
gnm_taxon_txt   = '/Users/songweizhi/Desktop/DateArTree/01_genome_selection/0_genome_metadata/gtdbtk.ar53.summary.tsv'

# file out
itol_label_out  = '/Users/songweizhi/Desktop/DateArTree/01_genome_selection/0_genome_metadata/iTOL_gnm_label.txt'


# read in gnm_quality_txt
gnm_quality_dict = dict()
for each_gnm in open(gnm_quality_txt):
    each_gnm_split = each_gnm.strip().split('\t')
    gnm_id = each_gnm_split[0]
    gnm_cpl = each_gnm_split[1]
    gnm_ctm = each_gnm_split[2]
    gnm_quality_dict[gnm_id] = '%s__%s' % (gnm_cpl, gnm_ctm)

# read in gnm_taxon_txt
max_tax_str_len = 0
gnm_taxon_dict = dict()
for each_gnm in open(gnm_taxon_txt):
    each_gnm_split = each_gnm.strip().split('\t')
    gnm_id = each_gnm_split[0]
    gnm_taxon = each_gnm_split[1]
    gnm_taxon_split = gnm_taxon.split(';')
    gnm_taxon_pco = ';'.join(gnm_taxon_split[1:4])
    if len(gnm_taxon_pco) > max_tax_str_len:
        max_tax_str_len = len(gnm_taxon_pco)
    gnm_taxon_dict[gnm_id] = gnm_taxon_pco

# read in gnm_size_txt
gnm_size_mbp_dict = dict()
for each_gnm in open(gnm_size_txt):
    each_gnm_split = each_gnm.strip().split('\t')
    gnm_id = each_gnm_split[0]
    gnm_size_mbp = each_gnm_split[3]
    gnm_size_mbp_dict[gnm_id] = gnm_size_mbp

# write out
itol_label_out_handle = open(itol_label_out, 'w')
itol_label_out_handle.write('LABELS\nSEPARATOR TAB\nDATA\n')
for each_gnm in open(gnm_id_txt):
    gnm_id = each_gnm.strip()
    gnm_quality_str = gnm_quality_dict.get(gnm_id, 'NA__NA')
    gnm_quality_str = gnm_quality_str + '_'*(12 - len(gnm_quality_str))
    gnm_size_str = gnm_size_mbp_dict.get(gnm_id, 'NA')
    gnm_size_str = gnm_size_str + '_'*(4 - len(gnm_size_str))

    gnm_taxon_str = gnm_taxon_dict.get(gnm_id, 'NA')
    gnm_taxon_str = gnm_taxon_str + '_'*(max_tax_str_len - len(gnm_taxon_str))
    itol_label_out_handle.write('%s\t%s__%s__%s__%s\n' % (gnm_id, gnm_taxon_str, gnm_quality_str, gnm_size_str, gnm_id))
itol_label_out_handle.close()
