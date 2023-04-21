
prokka_prefix_txt = '/Users/songweizhi/Desktop/DateArTree/01_genome_selection/prokka_prefix.txt'

for each_gnm in open(prokka_prefix_txt):
    each_gnm_split = each_gnm.strip().split('\t')
    gnm_file       = each_gnm_split[1]
    prefix         = each_gnm_split[0]
    prokka_cmd     = 'prokka --force --compliant --cpus 1 --kingdom Archaea --prefix %s --locustag %s --strain %s --outdir %s_prokka_wd %s' % (prefix, prefix, prefix, prefix, gnm_file)
    print(prokka_cmd)
