import os
import glob
from Bio import SeqIO


def sep_path_basename_ext(file_in):

    f_path, f_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'
    f_base, f_ext = os.path.splitext(f_name)
    f_ext = f_ext[1:]

    return f_name, f_path, f_base, f_ext


file_dir = '/Users/songweizhi/Desktop/999/03_AOA_genomes_1369_dRep85_263_GTDB_SCG_best50'
file_ext = 'fa'
gnm_txt  = '/Users/songweizhi/Desktop/999/03_AOA_genomes_1369_dRep85_255.txt'
op_dir   = '/Users/songweizhi/Desktop/999/03_AOA_genomes_1369_dRep85_255_GTDB_SCG_best50'


gnm_set = set()
for each in open(gnm_txt):
    gnm_set.add(each.strip())


file_re = '%s/*.%s' % (file_dir, file_ext)
file_list = glob.glob(file_re)

for each in file_list:
    f_name, f_path, f_base, f_ext = sep_path_basename_ext(each)
    op_fa = '%s/%s' % (op_dir, f_name)
    op_fa_handle = open(op_fa, 'w')
    for each_seq in SeqIO.parse(each, 'fasta'):
        seq_id = each_seq.id
        gnm_id = '_'.join(each_seq.id.split('_')[:-1])
        if gnm_id in gnm_set:
            op_fa_handle.write('>%s\n' % seq_id)
            op_fa_handle.write('%s\n' % str(each_seq.seq))
        else:
            print(seq_id)
    op_fa_handle.close()
