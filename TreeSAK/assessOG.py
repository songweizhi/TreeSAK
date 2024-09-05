import os
import glob
from Bio import SeqIO


def sep_path_basename_ext(file_in):

    f_path, f_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'
    f_base, f_ext = os.path.splitext(f_name)

    return f_name, f_path, f_base, f_ext[1:]


def get_gnm_og_cov(og_dir, og_ext, og_cov_txt):

    og_file_re   = '%s/*.%s' % (og_dir, og_ext)
    og_file_list = glob.glob(og_file_re)

    gnm_to_og_dict = dict()
    for og_file in og_file_list:
        _, _, og_id, _ = sep_path_basename_ext(og_file)
        for each_seq in SeqIO.parse(og_file, 'fasta'):
            seq_id = each_seq.id
            gnm_id = '_'.join(seq_id.split('_')[:-1])
            if gnm_id not in gnm_to_og_dict:
                gnm_to_og_dict[gnm_id] = set()
            gnm_to_og_dict[gnm_id].add(og_id)

    og_cov_txt_handle = open(og_cov_txt, 'w')
    for each_gnm in sorted(list(gnm_to_og_dict.keys())):
        gnm_og_set = gnm_to_og_dict[each_gnm]
        og_cov = len(gnm_og_set)*100/len(og_file_list)
        og_cov = float("{0:.2f}".format(og_cov))
        og_cov_txt_handle.write('%s\t%s\n' % (each_gnm, og_cov))
    og_cov_txt_handle.close()


og_dir      = '/Users/songweizhi/Desktop/OrthologousGroupsFasta_cov95'
og_ext      = 'fa'
og_cov_txt  = '/Users/songweizhi/Desktop/gnm_og_cov.txt'

get_gnm_og_cov(og_dir, og_ext, og_cov_txt)

