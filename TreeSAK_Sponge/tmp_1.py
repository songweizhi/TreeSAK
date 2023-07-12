
meta_final                              = '/Users/songweizhi/Documents/Research/Sponge_Hologenome/0_metadata/0_metadata_final.txt'

col_index = {}
for each_gnm in open(meta_final):
    each_gnm_split = each_gnm.strip().split('\t')
    if each_gnm.startswith('Genome\t'):
        #print(each_gnm_split)
        col_index = {key: i for i, key in enumerate(each_gnm_split)}
    else:
        gnm_id               = each_gnm_split[col_index['Genome']]
        gnm_host             = each_gnm_split[col_index['Host']]
        gnm_biosample        = each_gnm_split[col_index['Biosample']]
        gnm_taxon_str_split  = each_gnm_split[col_index['Taxon']].split(';')
        host_taxon_str_split = each_gnm_split[col_index['Host_taxon']].split(';')
        if gnm_host not in ['nonsponge', 'na']:
            #print('%s\t%s' % (gnm_id, gnm_biosample))
            pass
            #print(gnm_biosample)


import xml.etree.ElementTree as ET
tree = ET.parse('/Users/songweizhi/Desktop/bb.xml')
root = tree.getroot()

print('####################')

for child in root.iter():

    if len(child.attrib) == 0:
        pass
        #print(child.tag, child.text, sep='\t')
    else:
        print(child.tag, child.attrib, child.text, sep='\t')

        print(child.attrib.items())
        print()
        #print(child.attrib['host'])


print('####################')
#
# for rank in root.iter('SampleData'):
#     print(rank.text)
#     attrib = rank.attrib
#     print(len(attrib))

