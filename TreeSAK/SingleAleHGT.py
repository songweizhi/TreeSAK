import os
import argparse


def sep_path_basename_ext(file_in):

    # separate path and file name
    f_path, file_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'

    # separate file basename and extension
    f_base, f_ext = os.path.splitext(file_name)

    return f_path, f_base, f_ext


SingleAleHGT_usage = '''
============================================ SingleAleHGT example commands ============================================

TreeSAK SingleAleHGT -i concatenated.fasta -s genome.treefile -fc 0.3 -c genome_taxon.txt -color phylum_color.txt -api S1kZZuDHc0d5M7J5vLnUNQ -t 9 -f -o demo_SingleAleHGT_wd

=======================================================================================================================
'''

def SingleAleHGT(args):

    faa_in                      = args['faa']
    msa_in                      = args['msa']
    op_dir                      = args['o']
    genome_tree_file_rooted     = args['s']
    API_key                     = args['api']
    hgt_freq_cutoff             = args['fc']
    ar_phylum_color_code_txt    = args['color']
    genome_taxon_txt            = args['c']
    force_overwrite             = args['f']
    trim_msa                    = args['trim']
    docker_image                = args['docker']
    num_threads                 = args['t']

    ######################################## check input files #######################################

    # if docker_image is True, check if docker is activated
    if (faa_in is not None) and (msa_in is None):
        f_path, f_base, f_ext = sep_path_basename_ext(faa_in)
    elif (faa_in is None) and (msa_in is not None):
        f_path, f_base, f_ext = sep_path_basename_ext(msa_in)
    else:
        print('Please specify either -faa or -msa, program exited!')
        exit()

    ######################################## define file name ########################################

    ale1_op_dir = '%s/ALE1_op_dir'      % op_dir
    ale2_op_dir = '%s/ALE2_op_dir'      % op_dir
    ale4_op_dir = '%s/ALE4_op_dir'      % op_dir
    log_txt     = '%s/log.txt'          % op_dir
    msa_file    = '%s/%s.aln'           % (ale1_op_dir, f_base)
    msa_trimmed = '%s/%s_trimmed.aln'   % (ale1_op_dir, f_base)
    tree_prefix = '%s/%s'               % (ale1_op_dir, f_base)

    ###################################### create output folder ######################################

    if os.path.isdir(op_dir) is True:
        if force_overwrite is True:
            os.system('rm -r %s' % op_dir)
        else:
            print('%s exist, program exited!' % op_dir)
            exit()
    os.mkdir(op_dir)
    os.mkdir(ale1_op_dir)

    ##################################################################################################

    # run mafft-einsi
    if (faa_in is not None) and (msa_in is None):
        mafft_cmd = 'mafft-einsi --thread %s --quiet %s > %s' % (num_threads, faa_in, msa_file)

        with open(log_txt, 'a') as log_txt_handle:
            log_txt_handle.write(mafft_cmd + '\n')
        os.system(mafft_cmd)
        msa_file_for_next_step = msa_file
    else:
        msa_file_for_next_step = msa_in

    # run trimal
    if trim_msa is True:
        trimal_cmd = 'trimal -in %s -out %s -automated1' % (msa_file_for_next_step, msa_trimmed)
        with open(log_txt, 'a') as log_txt_handle:
            log_txt_handle.write(trimal_cmd + '\n')
        os.system(trimal_cmd)
        iqtree2_cmd = 'iqtree2 -m LG+G+I -bb 1000 --wbtl -nt %s -s %s -pre %s' % (num_threads, msa_trimmed, tree_prefix)
        with open(log_txt, 'a') as log_txt_handle:
            log_txt_handle.write(iqtree2_cmd + '\n')
        os.system(iqtree2_cmd)
    else:
        iqtree2_cmd = 'iqtree2 -m LG+G+I -bb 1000 --wbtl -nt %s -s %s -pre %s' % (num_threads, msa_file_for_next_step, tree_prefix)
        with open(log_txt, 'a') as log_txt_handle:
            log_txt_handle.write(iqtree2_cmd + '\n')
        os.system(iqtree2_cmd)

    # run ALE2
    ale2_cmd = 'TreeSAK ALE2 -i %s -s %s -t %s -f -runALE -docker %s -o %s' % (ale1_op_dir, genome_tree_file_rooted, num_threads, docker_image, ale2_op_dir)
    with open(log_txt, 'a') as log_txt_handle:
        log_txt_handle.write(ale2_cmd + '\n')
    os.system(ale2_cmd)

    # run ALE4
    ale4_cmd = 'TreeSAK ALE4 -i1 %s -i2 %s -c %s -color %s -o %s -fc %s -f -api %s' % (ale1_op_dir, ale2_op_dir, genome_taxon_txt, ar_phylum_color_code_txt, ale4_op_dir, hgt_freq_cutoff, API_key)
    with open(log_txt, 'a') as log_txt_handle:
        log_txt_handle.write(ale4_cmd + '\n')
    os.system(ale4_cmd)


if __name__ == '__main__':

    SingleAleHGT_parser = argparse.ArgumentParser()
    SingleAleHGT_parser.add_argument('-faa',    required=False, default=None,               help='input aa file, e.g., OMA0001.faa')
    SingleAleHGT_parser.add_argument('-msa',    required=False, default=None,               help='input MSA file, e.g., OMA0001.aln')
    SingleAleHGT_parser.add_argument('-o',      required=True,                              help='output dir, e.g., SingleAleHGT_wd')
    SingleAleHGT_parser.add_argument('-s',      required=True,                              help='rooted species tree')
    SingleAleHGT_parser.add_argument('-c',      required=True,                              help='genome_taxon, GTDB format')
    SingleAleHGT_parser.add_argument('-color',  required=True,                              help='phylum color code')
    SingleAleHGT_parser.add_argument('-fc',     required=False, type=float, default=0.5,    help='hgt_freq_cutoff, default: 0.5')
    SingleAleHGT_parser.add_argument('-mld',    required=False, type=int, default=5,        help='donor_node_min_leaf_num, default: 5')
    SingleAleHGT_parser.add_argument('-mlr',    required=False, type=int, default=5,        help='recipient_node_min_leaf_num, default: 5')
    SingleAleHGT_parser.add_argument('-trim',   required=False, action="store_true",        help='trim MSA')
    SingleAleHGT_parser.add_argument('-docker', required=False, default=None,               help='Docker image, if ALE was installed with Docker, e.g., gregmich/alesuite_new')
    SingleAleHGT_parser.add_argument('-itol',   required=False, default='batch_access_tmp', help='iTOL project_name, default: batch_access_tmp')
    SingleAleHGT_parser.add_argument('-api',    required=True,                              help='iTOL API key')
    SingleAleHGT_parser.add_argument('-t',      required=False, type=int, default=6,        help='number of threads, default: 6')
    SingleAleHGT_parser.add_argument('-f',      required=False, action="store_true",        help='force overwrite')
    args = vars(SingleAleHGT_parser.parse_args())
    SingleAleHGT(args)


'''

cd /Users/songweizhi/Desktop/DateArTree/01_HGT_ALE_with_OMA/ALE1_op_dir_OMA05484_OMA07484_trimmed
trimal -in ../ALE1_op_dir_OMA05484_OMA07484/concatenated.fasta -out concatenated.fasta -automated1
iqtree2 -m LG+G+I -bb 1000 --wbtl -nt 10 -s concatenated.fasta -pre OMA05484_OMA07484
cd /Users/songweizhi/Desktop/DateArTree/01_HGT_ALE_with_OMA
TreeSAK ALE2 -i ALE1_op_dir_OMA05484_OMA07484_trimmed -s genome_tree.newick -t 10 -f -runALE -docker gregmich/alesuite_new -o ALE2_op_dir_OMA05484_OMA07484_trimmed
TreeSAK ALE4 -i1 ALE1_op_dir_OMA05484_OMA07484_trimmed -i2 ALE2_op_dir_OMA05484_OMA07484_trimmed -c genome_taxon.txt -color phylum_color.txt -o ALE4_op_dir_OMA05484_OMA07484_trimmed_0.01 -fc 0.01 -f -api S1kZZuDHc0d5M7J5vLnUNQ

cd /Users/songweizhi/Desktop/DateArTree/01_HGT_ALE_with_OMA
/usr/local/bin/python3.7 /Users/songweizhi/PycharmProjects/TreeSAK/TreeSAK/SingleAleHGT.py -msa ALE1_op_dir_OMA05484_OMA07484_trimmed/concatenated.fasta -s genome_tree_rooted_noEU.treefile -fc 0.3 -c genome_taxon.txt -color phylum_color.txt -api S1kZZuDHc0d5M7J5vLnUNQ -t 9 -f -o demo_SingleAleHGT_wd -trim

cd /Users/songweizhi/Desktop/DateArTree/01_HGT_ALE_with_OMA/demo_SingleAleHGT_wd
TreeSAK ALE2 -i ALE1_op_dir -s ../genome_tree.newick -t 10 -f -runALE -docker gregmich/alesuite_new -o ALE2_op_dir
TreeSAK ALE4 -i1 ALE1_op_dir_OMA05484_OMA07484_trimmed -i2 ALE2_op_dir_OMA05484_OMA07484_trimmed -c genome_taxon.txt -color phylum_color.txt -o ALE4_op_dir_OMA05484_OMA07484_trimmed_0.01 -fc 0.01 -f -api S1kZZuDHc0d5M7J5vLnUNQ

/usr/local/bin/python3.7 /Users/songweizhi/PycharmProjects/TreeSAK/TreeSAK/SingleAleHGT.py -o demo_SingleAleHGT_wd -msa ALE1_op_dir/OMA15312.aln -s genome_tree_rooted_noEU.treefile -fc 0.3 -c genome_taxon.txt -color phylum_color.txt -api S1kZZuDHc0d5M7J5vLnUNQ -t 10 -f -trim -docker gregmich/alesuite_new
/usr/local/bin/python3.7 /Users/songweizhi/PycharmProjects/TreeSAK/TreeSAK/SingleAleHGT.py -o OMA01402_ALE_HGT_wd -msa ALE1_op_dir/OMA01402.aln -s genome_tree_rooted_noEU.treefile -fc 0.3 -c genome_taxon.txt -color phylum_color.txt -api S1kZZuDHc0d5M7J5vLnUNQ -t 10 -f -trim -docker gregmich/alesuite_new
/usr/local/bin/python3.7 /Users/songweizhi/PycharmProjects/TreeSAK/TreeSAK/SingleAleHGT.py -o OMA01402_ALE_HGT_wd_no_trim -msa ALE1_op_dir/OMA01402.aln -s genome_tree_rooted_noEU.treefile -fc 0.3 -c genome_taxon.txt -color phylum_color.txt -api S1kZZuDHc0d5M7J5vLnUNQ -t 10 -f -docker gregmich/alesuite_new

'''
