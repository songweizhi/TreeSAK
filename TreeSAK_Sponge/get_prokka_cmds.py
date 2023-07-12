
selected_genome_txt = '/Users/songweizhi/Desktop/aaa.txt'
gnm_dir             = '/home-user/wzsong/DateArTree_GTDB_r214_2/2_combined_genomes'
gnm_ext             = 'fna'
prokka_cmds_txt     = '/Users/songweizhi/Desktop/aaa_prokka.txt'

prokka_cmds_txt_handle = open(prokka_cmds_txt, 'w')
for each_gnm in open(selected_genome_txt):
    gnm_id = each_gnm.strip().split('.fna')[0]
    pwd_gnm_file = '%s/%s.%s' % (gnm_dir, gnm_id, gnm_ext)
    prokka_cmd     = 'prokka --force --compliant --cpus 1 --kingdom Archaea --prefix %s --locustag %s --strain %s --outdir %s_prokka_wd %s' % (gnm_id, gnm_id, gnm_id, gnm_id, pwd_gnm_file)
    prokka_cmds_txt_handle.write(prokka_cmd + '\n')
prokka_cmds_txt_handle.close()

