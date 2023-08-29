
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


for each_gnm in open('/Users/songweizhi/Desktop/195_genomes.txt'):
    gnm_id = each_gnm.strip()
    host_taxon = gnm_host_dict[gnm_id]
    if host_taxon not in ['nonsponge', 'na']:
        print(host_taxon)
