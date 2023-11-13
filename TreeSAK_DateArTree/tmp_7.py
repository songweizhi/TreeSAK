
interested_ko_txt   = '/Users/songweizhi/Desktop/interested_ko.txt'
locus_to_ko_txt     = '/Users/songweizhi/Desktop/Bradyrhizobium_neotropicale_HKCCT5'


interested_ko_set = set()
for each_ko in open(interested_ko_txt):
    interested_ko_set.add(each_ko.strip())


for each_locus in open(locus_to_ko_txt):
    each_locus_split = each_locus.strip().split('\t')
    ko_list = each_locus_split[1:]

    interested_locus = False
    for each_k in ko_list:
        if each_k in interested_ko_set:
            interested_locus = True

    if interested_locus is True:
        print(each_locus.strip())
print(interested_ko_set)

