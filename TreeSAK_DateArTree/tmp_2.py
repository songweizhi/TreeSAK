import os
import glob


og_file_dir             = '/Users/songweizhi/Desktop/DateArTree/01_genome_selection/0done_Tom_Williams_2020-39_archaeal_genomes/williams2019/scos'
concatenated_aln_fasta  = '/Users/songweizhi/Desktop/DateArTree/01_genome_selection/0done_Tom_Williams_2020-39_archaeal_genomes/williams2019/concatenates/43genes_concat.fasta'
blastp_op_dir           = '/Users/songweizhi/Desktop/DateArTree/01_genome_selection/0done_Tom_Williams_2020-39_archaeal_genomes/williams2019/scos_blastp'
blastp_cmd_txt          = '/Users/songweizhi/Desktop/DateArTree/01_genome_selection/0done_Tom_Williams_2020-39_archaeal_genomes/williams2019/scos_blastp_cmds.txt'



og_file_re = '%s/*.fa' % og_file_dir
og_file_list = [os.path.basename(file_name) for file_name in glob.glob(og_file_re)]

blastp_cmd_txt_handle = open(blastp_cmd_txt, 'w')
for og_file in og_file_list:
    pwd_og_file   = '%s/%s' % (og_file_dir, og_file)
    pwd_blastp_op = '%s/%s.txt' % (blastp_op_dir, og_file)
    blastp_cmd    = 'blastp -query %s -subject %s -out %s -outfmt 6' % (pwd_og_file, concatenated_aln_fasta, pwd_blastp_op)
    blastp_cmd_txt_handle.write(blastp_cmd + '\n')
blastp_cmd_txt_handle.close()

# cd /Users/songweizhi/Desktop/DateArTree/01_genome_selection/0done_Tom_Williams_2020-39_archaeal_genomes/williams2019
# BioSAK exe_cmds -c scos_blastp_cmds.txt -t 10


og_num = 0
for og_file in og_file_list:
    pwd_blastp_op = '%s/%s.txt' % (blastp_op_dir, og_file)
    found_in = 0
    for each_line in open(pwd_blastp_op):
        each_line_split = each_line.strip().split('\t')
        iden = float(each_line_split[2])
        aln_len = int(each_line_split[3])
        evalue = float(each_line_split[10])
        #if (evalue <= 0.001) and (iden >= 90):
        if (iden >= 85) and (aln_len >= 50):
            found_in += 1
    if found_in >= 41:
        print(og_file)
        og_num += 1
print(og_num)


bac_num = 10
total   = 92
ar      = 0
eu      = 10
ar_eu   = 82
