
prokaryotes_txt = '/Users/songweizhi/Desktop/prokaryotes.txt'
gnm_id_txt = '/Users/songweizhi/Desktop/2.txt'
# meta_1     = '/Users/songweizhi/Documents/Research/Sponge_Hologenome/0_metadata/0_metadata.txt'



def get_gnm_to_biosample_dict(meta_1, prokaryotes_txt):
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

    return gnm_id_to_biosample_dict

gnm_id_to_biosample_dict = get_gnm_to_biosample_dict(gnm_id_txt, prokaryotes_txt)
print(gnm_id_to_biosample_dict)
print(len(gnm_id_to_biosample_dict))


