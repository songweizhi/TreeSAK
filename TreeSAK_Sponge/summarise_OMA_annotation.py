import os
import glob


def sep_path_basename_ext(file_in):

    f_path, file_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'
    f_base, f_ext = os.path.splitext(file_name)

    return f_path, f_base, f_ext


########################################################################################################################

file_dir = '/Users/songweizhi/Desktop/ALE1_arcog_stats_copy_pct'
file_ext = 'txt'

summary_txt = '/Users/songweizhi/Desktop/ALE1_arcog_stats_copy_pct.txt'

########################################################################################################################

processed_og_num = 0
summary_txt_handle = open(summary_txt, 'w')
file_re = '%s/*.%s' % (file_dir, file_ext)
file_list = glob.glob(file_re)
for each_file in file_list:
    f_path, f_base, f_ext = sep_path_basename_ext(each_file)
    og_id = f_base.replace('_arcog_stats_copy_pct', '')
    cog_id_set = set()
    matched_cog_list = []
    matched_des_set = set()
    for each_line in open(each_file):
        if not each_line.startswith('arCOG\t'):
            each_line_split = each_line.strip().split('\t')
            cog_id = each_line_split[0]
            gene_pct = float(each_line_split[1])
            cog_des = each_line_split[2]
            if cog_des != 'Uncharacterized protein':
                matched_cog_list.append(each_line_split)
                cog_id_set.add(cog_id)
                matched_des_set.add(cog_des)

    if len(matched_des_set) == 0:
        processed_og_num += 1
        pass
    else:
        if len(matched_cog_list) == 1:
            summary_txt_handle.write('%s\t%s\n' % (og_id, '\t'.join(matched_cog_list[0])))
            processed_og_num += 1
        else:
            if len(matched_des_set) == 1:
                cog_id_list = []
                total_gene_pct = 0
                cog_des = ''
                for each in matched_cog_list:
                    cog_id = each[0]
                    gene_pct = float(each[1])
                    cog_des = each[2]
                    cog_id_list.append(cog_id)
                    total_gene_pct += gene_pct

                if total_gene_pct > 100:
                    total_gene_pct = 100.0

                summary_txt_handle.write('%s\t%s\t%s\t%s\n' % (og_id, ','.join(cog_id_list), total_gene_pct, cog_des))
                processed_og_num += 1
            else:
                max_cog_id = ''
                max_gene_pct = 0
                max_cog_des = ''
                for each in matched_cog_list:
                    cog_id = each[0]
                    gene_pct = float(each[1])
                    cog_des = each[2]
                    if gene_pct > max_gene_pct:
                        max_gene_pct = gene_pct
                        max_cog_id = cog_id
                        max_cog_des = cog_des

                if max_gene_pct >= 30:
                    summary_txt_handle.write('%s\t%s\t%s\t%s\n' % (og_id, max_cog_id, max_gene_pct, max_cog_des))
                    processed_og_num += 1
                else:
                    print(og_id)
                    for each in matched_cog_list:
                        cog_id = each[0]
                        gene_pct = float(each[1])
                        cog_des = each[2]
                        print(each)

summary_txt_handle.close()
print(len(file_list))
print(processed_og_num)
