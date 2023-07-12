import os


def get_gnm_to_host_dict(gnm_host_txt, gnm_host_txt_manual):

    gnm_to_host_dict = dict()
    for each_gnm in open(gnm_host_txt):
        each_gnm_split = each_gnm.strip().split('\t')
        gnm_id = each_gnm_split[0]
        gnm_host = each_gnm_split[1]
        gnm_to_host_dict[gnm_id] = gnm_host
    for each_gnm in open(gnm_host_txt_manual):
        each_gnm_split = each_gnm.strip().split('\t')
        gnm_id = each_gnm_split[0]
        gnm_host = each_gnm_split[1]
        gnm_to_host_dict[gnm_id] = gnm_host

    return gnm_to_host_dict


def get_gnm_to_biosample_dict(meta_1, prokaryotes_txt, gtdb_meta_data_txt):

    gnm_set = set()
    gnm_set_no_version = set()
    no_version_to_id_dict = dict()
    for each_gnm in open(meta_1):
        if not each_gnm.startswith('Genome\t'):
            each_gnm_split = each_gnm.strip().split('\t')
            gnm_id = each_gnm_split[0]
            if ('GCA' in gnm_id) or ('GCF' in gnm_id):
                gnm_id_no_version = gnm_id.split('.')[0]
                gnm_set.add(each_gnm.strip())
                gnm_set_no_version.add(gnm_id_no_version)
                no_version_to_id_dict[gnm_id_no_version] = gnm_id

    # read in prokaryotes_txt
    col_index = {}
    gnm_id_to_biosample_dict = dict()
    for genome_record in open(prokaryotes_txt):
        genome_record_split = genome_record.strip().split('\t')
        if genome_record.startswith('#Organism/Name'):
            col_index = {key: i for i, key in enumerate(genome_record_split)}
        else:
            assembly_id = genome_record_split[col_index['Assembly Accession']]
            assembly_id_no_version = assembly_id.split('.')[0]
            sample_id = genome_record_split[col_index['BioSample Accession']]

            assembly_id_no_version_GCF = assembly_id_no_version
            assembly_id_no_version_GCF = assembly_id_no_version_GCF.replace('GCA', 'GCF')

            if (assembly_id_no_version in gnm_set_no_version):
                gnm_id = no_version_to_id_dict[assembly_id_no_version]
                gnm_id_to_biosample_dict[gnm_id] = sample_id
            elif assembly_id_no_version_GCF in gnm_set_no_version:
                gnm_id = no_version_to_id_dict[assembly_id_no_version_GCF]
                gnm_id_to_biosample_dict[gnm_id] = sample_id

    # read in gtdb_meta_data_txt
    col_index = {}
    for each_line in open(gtdb_meta_data_txt):
        each_line_split = each_line.strip().split('\t')

        if each_line.startswith('accession\t'):
            col_index = {key: i for i, key in enumerate(each_line_split)}
        else:
            assembly_id = each_line_split[col_index['accession']]
            if 'GB_GCA_' in assembly_id:
                assembly_id = assembly_id.replace('GB_GCA_', 'GCA_')
            elif 'RS_GCF_' in assembly_id:
                assembly_id = assembly_id.replace('RS_GCF_', 'GCF_')

            biosample_id = each_line_split[col_index['ncbi_biosample']]

            gnm_id_to_biosample_dict[assembly_id] = biosample_id

    return gnm_id_to_biosample_dict


########################################################################################################################

# file in
meta_1                                  = '/Users/songweizhi/Documents/Research/Sponge_Hologenome/0_metadata/0_metadata.txt'
meta_db_genome_bad_ones                 = '/Users/songweizhi/Documents/Research/Sponge_Hologenome/0_metadata/meta_db_genome_bad_ones.txt'
sponge_associated_genome_taxon_r214_txt = '/Users/songweizhi/Documents/Research/Sponge_Hologenome/0_metadata/sponge_associated_genome_GTDB_r214.ar53.summary.tsv'
sponge_associated_genome_quality_txt    = '/Users/songweizhi/Documents/Research/Sponge_Hologenome/0_metadata/sponge_associated_genome_quality.txt'
gnm_size_txt                            = '/Users/songweizhi/Documents/Research/Sponge_Hologenome/0_metadata/gnm_size.txt'
gnm_host_txt                            = '/Users/songweizhi/Documents/Research/Sponge_Hologenome/0_metadata/Sponge_MAGs_1677_host.txt'
gnm_host_txt_manual                     = '/Users/songweizhi/Documents/Research/Sponge_Hologenome/0_metadata/manually_added_genome_to_host.csv'
biosample_metadata_txt                  = '/Users/songweizhi/Documents/Research/Sponge_Hologenome/0_metadata/Biosample_metadata_reformatted.txt'
prokaryotes_txt                         = '/Users/songweizhi/Documents/Research/Sponge_Hologenome/0_metadata/prokaryotes.txt'
sponge_taxonomy_txt                     = '/Users/songweizhi/Documents/Research/Sponge_Hologenome/0_metadata/Sponge_full_lineage_GTDB_format.txt'
meta_db_genome                          = '/Users/songweizhi/Documents/Research/Sponge_Hologenome/1_get_GTDB_MAGs/Ar_50_5_taxonomy.txt'
gtdb_meta_data_txt                      = '/Users/songweizhi/DB/GTDB_r214/ar53_metadata_r214.tsv'

# file out
meta_final                              = '/Users/songweizhi/Documents/Research/Sponge_Hologenome/0_metadata/0_metadata_final.txt'

########################################################################################################################

# read in sponge taxonomy
sponge_genus_to_taxonomy_dict = dict()
for each_line in open(sponge_taxonomy_txt):
    each_line_split = each_line.strip().split(';')
    sponge_g = each_line_split[5]
    for each_r in each_line_split:
        if each_r.startswith('g__'):
            sponge_g = each_r
    tax_str = each_line.split(sponge_g)[0] + sponge_g
    sponge_genus_to_taxonomy_dict[sponge_g] = tax_str

# read in biosample metadata
biosample_metadata_dict = dict()
for each_sample in open(biosample_metadata_txt):
    each_sample_split = each_sample.strip().split('\t')
    biosample_metadata_dict[each_sample_split[0]] = each_sample_split[1]

# get gnm_to_biosample_dict
gnm_to_biosample_dict = get_gnm_to_biosample_dict(meta_1, prokaryotes_txt, gtdb_meta_data_txt)

# read in host information
gnm_to_host_dict = get_gnm_to_host_dict(gnm_host_txt, gnm_host_txt_manual)

# read in gnm_size_txt
gnm_size_dict = dict()
for each_gnm in open(gnm_size_txt):
    if not each_gnm.startswith('Genome\t'):
        each_gnm_split  = each_gnm.strip().split('\t')
        gnm_id          = each_gnm_split[0]
        gnm_size        = each_gnm_split[1]
        gnm_size_dict[gnm_id] = gnm_size

# read in sponge_associated_genome_taxon_r214_txt
sponge_gnm_to_taxon_dict = dict()
for each_gnm in open(sponge_associated_genome_taxon_r214_txt):
    if not each_gnm.startswith('user_genome\t'):
        each_gnm_split  = each_gnm.strip().split('\t')
        gnm_id          = each_gnm_split[0]
        gnm_taxon       = each_gnm_split[1]
        sponge_gnm_to_taxon_dict[gnm_id] = gnm_taxon

# read in aponge_associated_genome_quality_txt
sponge_gnm_to_cpl_dict = dict()
sponge_gnm_to_ctm_dict = dict()
for each_gnm in open(sponge_associated_genome_quality_txt):
    if not each_gnm.startswith('Genome\t'):
        each_gnm_split  = each_gnm.strip().split('\t')
        gnm_id          = each_gnm_split[0]
        gnm_cpl         = each_gnm_split[1]
        gnm_ctm         = each_gnm_split[2]
        sponge_gnm_to_cpl_dict[gnm_id] = gnm_cpl
        sponge_gnm_to_ctm_dict[gnm_id] = gnm_ctm

# read in db genome metadata
db_gnm_to_cpl_dict = dict()
db_gnm_to_ctm_dict = dict()
db_gnm_to_taxon_dict = dict()
for each_gnm in open(meta_db_genome):
    if not each_gnm.startswith('Genome\t'):
        each_gnm_split  = each_gnm.strip().split('\t')
        gnm_id          = each_gnm_split[0]
        gnm_cpl         = each_gnm_split[1]
        gnm_ctm         = each_gnm_split[2]
        gnm_taxon       = each_gnm_split[4]
        db_gnm_to_cpl_dict[gnm_id]   = gnm_cpl
        db_gnm_to_ctm_dict[gnm_id]   = gnm_ctm
        db_gnm_to_taxon_dict[gnm_id] = gnm_taxon

gnm_with_host = 0
meta_final_handle = open(meta_final, 'w')
for each_gnm in open(meta_1):
    if each_gnm.startswith('Genome\t'):
        meta_final_handle.write('Genome\tSource\tAlias\tCompleteness\tContamination\tSize\tTaxon\tBiosample\tHost\tHost_taxon\n')
    else:
        each_gnm_split = each_gnm.strip().split('\t')
        gnm_id      = each_gnm_split[0]
        gnm_source  = each_gnm_split[1]
        gnm_alias   = each_gnm_split[2]

        # get gnm cpl
        if gnm_id in db_gnm_to_cpl_dict:
            gnm_cpl = db_gnm_to_cpl_dict[gnm_id]
        elif gnm_id in sponge_gnm_to_cpl_dict:
            gnm_cpl = sponge_gnm_to_cpl_dict[gnm_id]
        else:
            gnm_cpl = 'na'

        # get gnm ctm
        if gnm_id in db_gnm_to_ctm_dict:
            gnm_ctm = db_gnm_to_ctm_dict[gnm_id]
        elif gnm_id in sponge_gnm_to_ctm_dict:
            gnm_ctm = sponge_gnm_to_ctm_dict[gnm_id]
        else:
            gnm_ctm = 'na'

        # get gnm_size
        gnm_size = gnm_size_dict.get(gnm_id, 'na')

        # get gnm taxon
        if gnm_id in db_gnm_to_taxon_dict:
            gnm_taxon = db_gnm_to_taxon_dict[gnm_id]
        elif gnm_id in sponge_gnm_to_taxon_dict:
            gnm_taxon = sponge_gnm_to_taxon_dict[gnm_id]
        else:
            gnm_taxon = 'na'

        if gnm_taxon == 'na':
            gnm_genus = 'na'
        elif 'Undefined' in gnm_taxon:
            gnm_genus = 'na'
        else:
            gnm_genus = gnm_taxon.split(';')[5]

        # get sponge host
        if gnm_id in gnm_to_host_dict:
            gnm_host = gnm_to_host_dict[gnm_id]
            gnm_with_host += 1
        elif gnm_alias in gnm_to_host_dict:
            gnm_host = gnm_to_host_dict[gnm_alias]
            gnm_with_host += 1
        else:
            gnm_host = 'na'

        # get biosample_id
        biosample_id = gnm_to_biosample_dict.get(gnm_id, 'na')

        # get biosample_type
        gnm_host_to_write = biosample_metadata_dict.get(biosample_id, 'na')
        if gnm_host != 'na':
            gnm_host_to_write = gnm_host

        # get Host_taxon
        if gnm_host_to_write == 'na':
            host_taxon = 'na'
        elif gnm_host_to_write == 'nonsponge':
            host_taxon = 'nonsponge'
        else:
            gnm_host_g = gnm_host_to_write.split('_')[0]
            if '(coral)' in gnm_host_to_write:
                gnm_host_g += '(coral)'
            host_taxon = sponge_genus_to_taxonomy_dict.get('g__' + gnm_host_g, 'na')

        # write out
        to_write_str = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (gnm_id, gnm_source, gnm_alias, gnm_cpl, gnm_ctm, gnm_size, gnm_taxon, biosample_id, gnm_host_to_write, host_taxon)
        meta_final_handle.write(to_write_str + '\n')
meta_final_handle.close()


########################################################################################################################
########################################################################################################################
########################################################################################################################

# file in
host_group_color = '/Users/songweizhi/Documents/Research/Sponge_Hologenome/0_metadata/group_color.txt'
gnms_to_ignore   = {'CHO2_bin2'}

# file out
itol_file_dir           = '/Users/songweizhi/Documents/Research/Sponge_Hologenome/0_metadata_iTOL'
gnm_to_genus_txt        = '%s/gnm_to_genus.txt'           % itol_file_dir
gnm_to_host_g_txt       = '%s/gnm_to_host_1_g.txt'        % itol_file_dir
gnm_to_host_f_txt       = '%s/gnm_to_host_2_f.txt'        % itol_file_dir
gnm_to_host_o_txt       = '%s/gnm_to_host_3_o.txt'        % itol_file_dir
gnm_to_host_sc_txt      = '%s/gnm_to_host_4_sc.txt'       % itol_file_dir
gnm_to_host_c_txt       = '%s/gnm_to_host_5_c.txt'        % itol_file_dir

gnm_to_genus_txt_itol   = '%s/iTOL_gnm_to_genus.txt'      % itol_file_dir
gnm_to_host_g_txt_itol  = '%s/iTOL_gnm_to_host_1_g.txt'   % itol_file_dir
gnm_to_host_f_txt_itol  = '%s/iTOL_gnm_to_host_2_f.txt'   % itol_file_dir
gnm_to_host_o_txt_itol  = '%s/iTOL_gnm_to_host_3_o.txt'   % itol_file_dir
gnm_to_host_sc_txt_itol = '%s/iTOL_gnm_to_host_4_sc.txt'  % itol_file_dir
gnm_to_host_c_txt_itol  = '%s/iTOL_gnm_to_host_5_c.txt'   % itol_file_dir

gnm_to_genus_txt_handle = open(gnm_to_genus_txt, 'w')
gnm_to_host_g_txt_handle = open(gnm_to_host_g_txt, 'w')
gnm_to_host_f_txt_handle = open(gnm_to_host_f_txt, 'w')
gnm_to_host_o_txt_handle = open(gnm_to_host_o_txt, 'w')
gnm_to_host_sc_txt_handle = open(gnm_to_host_sc_txt, 'w')
gnm_to_host_c_txt_handle = open(gnm_to_host_c_txt, 'w')
col_index = {}
for each_gnm in open(meta_final):
    each_gnm_split = each_gnm.strip().split('\t')
    if each_gnm.startswith('Genome\t'):
        col_index = {key: i for i, key in enumerate(each_gnm_split)}
    else:
        gnm_id               = each_gnm_split[col_index['Genome']]

        if gnm_id not in gnms_to_ignore:
            gnm_taxon_str_split  = each_gnm_split[col_index['Taxon']].split(';')
            host_taxon_str_split = each_gnm_split[col_index['Host_taxon']].split(';')

            # write out to gnm_to_genus_txt
            gnm_genus = ''
            for each_gnm_r in gnm_taxon_str_split:
                if each_gnm_r.startswith('g__'):
                    gnm_genus = each_gnm_r
            gnm_to_genus_txt_handle.write('%s\t%s\n' % (gnm_id, gnm_genus))

            # write out host taxon
            if host_taxon_str_split == ['nonsponge']:
                gnm_to_host_g_txt_handle.write('%s\t%s\n' % (gnm_id, 'nonsponge'))
                gnm_to_host_f_txt_handle.write('%s\t%s\n' % (gnm_id, 'nonsponge'))
                gnm_to_host_o_txt_handle.write('%s\t%s\n' % (gnm_id, 'nonsponge'))
                gnm_to_host_sc_txt_handle.write('%s\t%s\n' % (gnm_id, 'nonsponge'))
                gnm_to_host_c_txt_handle.write('%s\t%s\n' % (gnm_id, 'nonsponge'))
            elif host_taxon_str_split == ['na']:
                gnm_to_host_g_txt_handle.write('%s\t%s\n' % (gnm_id, 'na'))
                gnm_to_host_f_txt_handle.write('%s\t%s\n' % (gnm_id, 'na'))
                gnm_to_host_o_txt_handle.write('%s\t%s\n' % (gnm_id, 'na'))
                gnm_to_host_sc_txt_handle.write('%s\t%s\n' % (gnm_id, 'na'))
                gnm_to_host_c_txt_handle.write('%s\t%s\n' % (gnm_id, 'na'))
            else:
                for each_host_r in host_taxon_str_split:
                    if each_host_r.startswith('g__'):
                        gnm_to_host_g_txt_handle.write('%s\t%s\n' % (gnm_id, each_host_r))
                    elif each_host_r.startswith('f__'):
                        gnm_to_host_f_txt_handle.write('%s\t%s\n' % (gnm_id, each_host_r))
                    elif each_host_r.startswith('o__'):
                        gnm_to_host_o_txt_handle.write('%s\t%s\n' % (gnm_id, each_host_r))
                    elif each_host_r.startswith('sc__'):
                        gnm_to_host_sc_txt_handle.write('%s\t%s\n' % (gnm_id, each_host_r))
                    elif each_host_r.startswith('c__'):
                        gnm_to_host_c_txt_handle.write('%s\t%s\n' % (gnm_id, each_host_r))
gnm_to_genus_txt_handle.close()
gnm_to_host_g_txt_handle.close()
gnm_to_host_f_txt_handle.close()
gnm_to_host_o_txt_handle.close()
gnm_to_host_sc_txt_handle.close()
gnm_to_host_c_txt_handle.close()

biosak_cmd_genus   = 'BioSAK iTOL -ColorRange -lg %s -lt Genus -o %s'          % (gnm_to_genus_txt,   gnm_to_genus_txt_itol)
biosak_cmd_host_g  = 'BioSAK iTOL -ColorStrip -lg %s -gc %s -lt Host_g -o %s'  % (gnm_to_host_g_txt,  host_group_color, gnm_to_host_g_txt_itol)
biosak_cmd_host_f  = 'BioSAK iTOL -ColorStrip -lg %s -gc %s -lt Host_f -o %s'  % (gnm_to_host_f_txt,  host_group_color, gnm_to_host_f_txt_itol)
biosak_cmd_host_o  = 'BioSAK iTOL -ColorStrip -lg %s -gc %s -lt Host_o -o %s'  % (gnm_to_host_o_txt,  host_group_color, gnm_to_host_o_txt_itol)
biosak_cmd_host_sc = 'BioSAK iTOL -ColorStrip -lg %s -gc %s -lt Host_sc -o %s' % (gnm_to_host_sc_txt, host_group_color, gnm_to_host_sc_txt_itol)
biosak_cmd_host_c  = 'BioSAK iTOL -ColorStrip -lg %s -gc %s -lt Host_c -o %s'  % (gnm_to_host_c_txt,  host_group_color, gnm_to_host_c_txt_itol)
os.system(biosak_cmd_genus)
os.system(biosak_cmd_host_g)
os.system(biosak_cmd_host_f)
os.system(biosak_cmd_host_o)
os.system(biosak_cmd_host_sc)
os.system(biosak_cmd_host_c)

########################################################################################################################
########################################################################################################################
########################################################################################################################
