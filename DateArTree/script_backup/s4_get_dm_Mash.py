import matplotlib.pyplot as plt
import numpy as np


mash_op_txt = '/Users/songweizhi/Desktop/DateArTree/01_genome_selection/all_vs_all_mash_1000000.txt'
dm_out      = '/Users/songweizhi/Desktop/DateArTree/01_genome_selection/all_vs_all_mash_1000000_dm.txt'


gnm_id_set = set()
pairwise_dist_dict = dict()
for each_ani in open(mash_op_txt):
    each_ani_split = each_ani.strip().split('\t')
    gnm1_id        = each_ani_split[0].split('.fna')[0]
    gnm2_id        = each_ani_split[1].split('.fna')[0]
    dist_value      = float(each_ani_split[2])
    key_str = '___'.join(sorted([gnm1_id, gnm2_id]))
    if key_str not in pairwise_dist_dict:
        pairwise_dist_dict[key_str] = dist_value
    else:
        pass
        #if abs(pairwise_dist_dict[key_str] - dist_value) >= 0.01:
        #print('%s\t%s\t%s\t%s' % (abs(pairwise_dist_dict[key_str] - dist_value),pairwise_dist_dict[key_str], dist_value, key_str))
    gnm_id_set.add(gnm1_id)
    gnm_id_set.add(gnm2_id)


gnm_id_list_sorted = sorted([i for i in gnm_id_set])


dm_out_handle = open(dm_out, 'w')
dm_out_handle.write('\t%s\n' % ('\t'.join(gnm_id_list_sorted)))
for gl in gnm_id_list_sorted:
    current_dist_list = [gl]
    for gr in gnm_id_list_sorted:
        gl_gr = '___'.join(sorted([gl, gr]))
        dist_value = pairwise_dist_dict[gl_gr]
        current_dist_list.append(dist_value)
    print(current_dist_list)
    current_dist_list_str = [str(i) for i in current_dist_list]
    dm_out_handle.write('%s\n' % ('\t'.join(current_dist_list_str)))
dm_out_handle.close()

# plt.hist(dist_value_list)
# plt.show()

