
gnm_tax_txt = '/Users/songweizhi/Desktop/0_metadata_final.txt'

gnm_id_list = []
for each_gnm in open(gnm_tax_txt):
    each_gnm_split = each_gnm.strip().split('\t')
    gnm_id_list.append(each_gnm_split[0])


for each_gnm in open('/Users/songweizhi/Desktop/aaa.txt'):
    gnm_id = each_gnm.strip()
    if gnm_id in gnm_id_list:
        print(gnm_id)