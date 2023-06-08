import os


dir = '/Users/songweizhi/Desktop/aaa'
dir = '/home-user/wzsong/DateArTree_GTDB_r214_2/OMA_wd2/Cache/AllAll'

sub_dir_list = next(os.walk(dir))[1]

total= 0
for each_f in sub_dir_list:
    pwd_sub_dir = '%s/%s' % (dir, each_f)
    sub_dir_list_2 = next(os.walk(pwd_sub_dir))[1]
    total += len(sub_dir_list_2)
print(total)

