import os
import glob
from Bio import SeqIO


def sep_path_basename_ext(file_in):

    # separate path and file name
    f_path, file_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'

    # separate file basename and extension
    f_base, f_ext = os.path.splitext(file_name)

    return f_path, f_base, f_ext


faa_dir         = '/Users/songweizhi/Documents/Research/Sponge_Hologenome/3_combined_genomes_faa'
GeneContent_dir = '/Users/songweizhi/Desktop/GeneContent'



gene_content_file_re = '%s/*.%s' % (GeneContent_dir, 'txt')
gene_content_file_list = glob.glob(gene_content_file_re)

for each_file in gene_content_file_list:

    _, f_base, _ = sep_path_basename_ext(each_file)
    if '(' in f_base:
        f_base = f_base.split('(')[0]

    gene_num_by_ale = 0
    for each_line in open(each_file):
        gene_num_by_ale += 1

    gene_num_by_faa = 0
    pwd_faa = '%s/%s.faa' % (faa_dir, f_base.replace('GCA', 'GCA_').replace('GCF', 'GCF_').replace('bin', '_bin').replace('CYMC67496', 'CYMC_67496'))
    for each_seq in SeqIO.parse(pwd_faa, 'fasta'):
        gene_num_by_faa += 1

    print('%s\t%s\t%s' % (f_base, gene_num_by_faa, gene_num_by_ale))

