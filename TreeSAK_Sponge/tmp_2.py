
og_to_des_txt           = '/Users/songweizhi/Documents/Research/Sponge_Hologenome/8_ALE_wd_all_OGs/ALE1_arcog_stats_copy_pct.txt'
Transfer_propensity_txt = '/Users/songweizhi/Documents/Research/Sponge_Hologenome/8_ALE_wd_all_OGs/ALE3_op_dir_c0.75/Transfer_propensity.txt'
op_txt                  = '/Users/songweizhi/Documents/Research/Sponge_Hologenome/8_ALE_wd_all_OGs/OG_transfer_propensity_with_description.txt'

og_to_des_dict = dict()
for each in open(og_to_des_txt):
    each_split = each.strip().split('\t')
    og_id = each_split[0]
    og_des = each_split[3]
    og_to_des_dict[og_id] = og_des

op_txt_handle = open(op_txt, 'w')
for each in open(Transfer_propensity_txt):
    if not each.startswith('OG\tTransfer_propensity'):
        each_split = each.strip().split('\t')
        og_id = each_split[0]
        propensity = float(each_split[1])
        od_desc = og_to_des_dict.get(og_id, 'NA')
        op_txt_handle.write('%s\t%s\t%s\n' % (og_id, propensity, od_desc))
op_txt_handle.close()
