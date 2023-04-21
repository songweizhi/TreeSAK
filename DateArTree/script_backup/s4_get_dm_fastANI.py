import matplotlib.pyplot as plt
import numpy as np

'''
fastANI --ql gnm_id_with_ext.txt --rl gnm_id_with_ext.txt -o all_vs_all_fastani.txt -t 12 --minFrag 1
fastANI --ql gnm_id_with_ext.txt --rl gnm_id_with_ext.txt -o all_vs_all_fastani.txt -t 12 --minFrag 1 --fragLen 1000
fastANI -q GCA_015521515.1.gtdb.fna -r GCA_015522845.1.gtdb.fna -o 000.txt -t 12 --minFrag 1
'''


fastani_op_txt = '/Users/songweizhi/Desktop/DateArTree/01_genome_selection/d__Archaea_o_rs_gnms_Betts38_Spang81_Williams39_Wu2_all_vs_all_fastani_1.txt'


gnm_id_set = set()
pairwise_ani_dict = dict()
for each_ani in open(fastani_op_txt):
    each_ani_split = each_ani.strip().split('\t')
    gnm1_id        = each_ani_split[0].split('.fna')[0]
    gnm2_id        = each_ani_split[1].split('.fna')[0]
    ani_value      = float(each_ani_split[2])
    key_str = '___'.join(sorted([gnm1_id, gnm2_id]))
    if key_str not in pairwise_ani_dict:
        pairwise_ani_dict[key_str] = ani_value
    else:
        if abs(pairwise_ani_dict[key_str] - ani_value) >= 5:
            #print('%s\t%s\t%s\t%s' % (abs(pairwise_ani_dict[key_str] - ani_value),pairwise_ani_dict[key_str], ani_value, key_str))
            pass
    gnm_id_set.add(gnm1_id)
    gnm_id_set.add(gnm2_id)

ani_value_list = []
for gl in gnm_id_set:
    for gr in gnm_id_set:
        gl_gr = '___'.join(sorted([gl, gr]))
        #ani = pairwise_ani_dict[gl_gr]
        ani = pairwise_ani_dict.get(gl_gr, 0)
        if ani >= 80:
            print(ani)
        ani_value_list.append(ani)



plt.hist(ani_value_list)
plt.show()