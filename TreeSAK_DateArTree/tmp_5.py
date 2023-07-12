import os
import glob


def sep_path_basename_ext(file_in):
    f_path, file_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'
    f_base, f_ext = os.path.splitext(file_name)
    return f_path, f_base, f_ext


file_dir = '/Users/songweizhi/Desktop/000'
file_ext = 'txt'
op_txt   = '/Users/songweizhi/Desktop/000.txt'


file_re = '%s/*.%s' % (file_dir, file_ext)
file_list = glob.glob(file_re)


op_txt_handle = open(op_txt, 'w')
cog_to_desc_dict = dict()
for each_file in sorted(file_list):
    f_path, f_base, f_ext = sep_path_basename_ext(each_file)
    og_id = f_base.split('_arcog_stats_copy')[0]

    cog_with_max_gene = ''
    max_gene_num = 0
    for each_line in open(each_file):
        if not each_line.startswith('arCOG\tcopy\tDescription'):
            each_line_split = each_line.strip().split('\t')
            if len(each_line_split) == 3:
                cog_id   = each_line_split[0]
                gene_num = int(each_line_split[1])
                cog_desc = each_line_split[2]
                if gene_num > max_gene_num:
                    max_gene_num = gene_num
                    cog_with_max_gene = cog_id
                    cog_to_desc_dict[cog_id] = cog_desc

    if cog_with_max_gene != '':
        print('%s\t%s\t%s' % (og_id, cog_with_max_gene, cog_to_desc_dict.get(cog_with_max_gene, '')))
        op_txt_handle.write('%s\t%s\t%s\n' % (og_id, cog_with_max_gene, cog_to_desc_dict.get(cog_with_max_gene, '')))
op_txt_handle.close()
