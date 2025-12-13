
para_txt = '/Users/songweizhi/Desktop/555/batch_itol_para.txt'

para_dict = dict()
for each_line in open(para_txt):
    if not each_line.startswith('#'):
        if len(each_line.strip()) > 0:
            para_without_comment = each_line.strip().split('#')[0].strip()
            para_without_comment_split = para_without_comment.split('\t')
            para_dict[para_without_comment_split[0]] = para_without_comment_split[1]

print(para_dict)