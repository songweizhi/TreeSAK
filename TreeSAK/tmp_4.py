
interested_gnm_set = set()
for each_gnm in open('/Users/songweizhi/Desktop/Nitrosopumilaceae_50_5_dRep97_195_id.txt'):
    interested_gnm_set.add(each_gnm.strip())
print(interested_gnm_set)
print(len(interested_gnm_set))


for each in open('/Users/songweizhi/Documents/Research/Sponge/2_metadata/0_metadata_final_renamed.txt'):
    each_split = each.strip().split('\t')
    if each_split[0] in interested_gnm_set:
        print('%s\t%s' % (each_split[0], each_split[6]))



