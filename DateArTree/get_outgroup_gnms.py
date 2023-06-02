
########################################################################################################################

# file in
gnm_taxon_txt               = '/Users/songweizhi/Desktop/DateArTree/0_HGT_ALE_r214/d__Archaea_o_rs_r214.txt'
outgroup_phyla_txt          = '/Users/songweizhi/Desktop/DateArTree/0_HGT_ALE_r214/outgroup_phyla.txt'

# file out
outgroup_gnm_txt            = '/Users/songweizhi/Desktop/DateArTree/0_HGT_ALE_r214/outgroup_gnm.txt'
outgroup_gnm_with_taxon_txt = '/Users/songweizhi/Desktop/DateArTree/0_HGT_ALE_r214/outgroup_gnm_with_taxon.txt'


'''
cd /Users/songweizhi/Desktop/DateArTree/0_HGT_ALE_r214
BioSAK iTOL -Labels -ll outgroup_gnm_with_taxon.txt -o outgroup_gnm_with_taxon_iTOL.txt
'''

########################################################################################################################

outgroup_phyla_set = set()
for each_p in open(outgroup_phyla_txt):
    outgroup_phyla_set.add(each_p.strip())

gnm_to_phylum_dict = dict()
outgroup_gnm_set = set()
for each_gnm in open(gnm_taxon_txt):
    if not each_gnm.startswith('Genome\t'):
        each_gnm_split = each_gnm.strip().split('\t')
        gnm_id = each_gnm_split[0]
        gnm_taxon = each_gnm_split[5]
        gnm_phylum = gnm_taxon.split(';')[1]
        if gnm_phylum in outgroup_phyla_set:
            outgroup_gnm_set.add(gnm_id)
            gnm_to_phylum_dict[gnm_id] = gnm_phylum

outgroup_gnm_txt_handle = open(outgroup_gnm_txt, 'w')
outgroup_gnm_with_taxon_txt_handle = open(outgroup_gnm_with_taxon_txt, 'w')
for each_gnm in outgroup_gnm_set:
    gnm_phylum = gnm_to_phylum_dict[each_gnm]
    outgroup_gnm_txt_handle.write('%s\n' % each_gnm)
    outgroup_gnm_with_taxon_txt_handle.write('%s\t%s\n' % (each_gnm, gnm_phylum))
outgroup_gnm_txt_handle.close()
outgroup_gnm_with_taxon_txt_handle.close()
