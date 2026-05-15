from Bio import SeqIO


# for each_line in open('3167.txt'):
#     og_id = each_line.strip()
#     faa_raw = '/project/spongeholobiont/Sponge_r232/AOA_2279_plus_r232_2383_OMA_230_OMA_wd/Output/the_rest_files/%s.faa' % og_id
#     faa_new = '/project/spongeholobiont/Sponge_r232/AOA_2279_plus_r232_2383_OMA_230_OMA_wd/Output/faa_3167/%s.faa'       % og_id
#     faa_new_handle = open(faa_new, 'w')
#     for each_seq in SeqIO.parse(faa_raw, 'fasta'):
#         seq_id = each_seq.id
#         gnm_id = '_'.join(seq_id.split('_')[:-1])
#         if gnm_id not in ['GCA029775995_1', 'GCA021162925_1']:
#             faa_new_handle.write('>%s\n' % seq_id)
#             faa_new_handle.write('%s\n'  % str(each_seq.seq))
#         else:
#             print(gnm_id)
#     faa_new_handle.close()


for each_line in open('/Users/songweizhi/Desktop/640.txt'):
    og_id = each_line.strip()
    cmd = 'mafft-einsi --thread 3 --quiet %s.faa > %s.aln; iqtree2 -m LG+G+I -bb 1000 --wbtl -nt 3 -s %s.aln -pre %s' % (og_id,og_id,og_id,og_id)
    print(cmd)


# import os
# import glob
#
# file_dir = '/Users/songweizhi/Desktop/log_dir'
# file_ext = 'log'
# file_re = '%s/*.%s' % (file_dir, file_ext)
# file_list = glob.glob(file_re)
#
# for each in file_list:
#     with open(each) as file:
#         last_line = ""
#         for line in file:
#             last_line = line
#         if 'Date and Time:' in last_line.strip():
#             print(each)
