
selected_genome_txt = '/Users/songweizhi/Desktop/d__Archaea_o_rs.txt'
gnm_dir             = '/home-user/wzsong/DateArTree_GTDB_r214/01_genome_selection/d__Archaea_o_rs'
gnm_ext             = 'fna'
prokka_cmds_txt     = '/Users/songweizhi/Desktop/prokka_cmds.txt'

prokka_cmds_txt_handle = open(prokka_cmds_txt, 'w')
for each_gnm in open(selected_genome_txt):
    if not each_gnm.startswith('Genome\tCompleteness\tContamination'):

        # get gnm_id
        if '\t' in each_gnm:
            gnm_id = each_gnm.strip().split('\t')[0]
        elif ' ' in each_gnm:
            gnm_id = each_gnm.strip().split(' ')[0]
        else:
            gnm_id = each_gnm.strip()

        pwd_gnm_file = '%s/%s.%s' % (gnm_dir, gnm_id, gnm_ext)
        prokka_cmd     = 'prokka --force --compliant --cpus 1 --kingdom Archaea --prefix %s --locustag %s --strain %s --outdir %s_prokka_wd %s' % (gnm_id, gnm_id, gnm_id, gnm_id, pwd_gnm_file)
        print(prokka_cmd)
        prokka_cmds_txt_handle.write(prokka_cmd + '\n')
prokka_cmds_txt_handle.close()
