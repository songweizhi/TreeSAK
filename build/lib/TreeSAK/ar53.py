import os
import glob
import argparse
from Bio import SeqIO


ar53_usage = '''
========================= ar53 example commands =========================

TreeSAK ar53 -i seq_dir -x fa -g interested_gnms.txt -o seq_dir_renamed

=========================================================================
'''

def sep_path_basename_ext(file_in):

    f_path, f_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'
    f_base, f_ext = os.path.splitext(f_name)
    f_ext = f_ext[1:]

    return f_name, f_path, f_base, f_ext


def ar53(args):

    marker_seq_dir      = args['i']
    marker_seq_ext      = args['x']
    interested_gnm_txt  = args['g']
    op_dir              = args['o']
    force_overwrite     = args['f']

    # create output folder
    if os.path.isdir(op_dir) is True:
        if force_overwrite is True:
            os.system('rm -r %s' % op_dir)
        else:
            print('%s exist, program exited!' % op_dir)
            exit()
    os.mkdir(op_dir)

    interested_gnm_set = set()
    if os.path.isfile(interested_gnm_txt) is True:
        for each_gnm in open(interested_gnm_txt):
            interested_gnm_set.add(each_gnm.strip())
        if len(interested_gnm_set) == 0:
            print('No genome provided in %s, program exited!' % interested_gnm_txt)
            exit()

    marker_seq_re   = '%s/*.%s' % (marker_seq_dir, marker_seq_ext)
    marker_seq_list = sorted(glob.glob(marker_seq_re))

    if len(marker_seq_list) == 0:
        print('No file found in %s, program exited!' % marker_seq_dir)
        exit()

    for marker_seq in marker_seq_list:

        f_name, _, _, _ = sep_path_basename_ext(marker_seq)
        pwd_op_file = '%s/%s' % (op_dir, f_name)

        pwd_op_file_handle = open(pwd_op_file, 'w')
        for each_seq in SeqIO.parse(marker_seq, 'fasta'):
            gnm_id = each_seq.id

            if os.path.isfile(interested_gnm_txt) is False:
                pwd_op_file_handle.write('>%s_XXX\n' % gnm_id)
                pwd_op_file_handle.write('%s\n' % str(each_seq.seq))
            else:
                if gnm_id in interested_gnm_set:
                    pwd_op_file_handle.write('>%s_XXX\n' % gnm_id)
                    pwd_op_file_handle.write('%s\n' % str(each_seq.seq))
        pwd_op_file_handle.close()

    print('Done!')


if __name__ == '__main__':

    ar53_parser = argparse.ArgumentParser()
    ar53_parser.add_argument('-i',   required=True,                          help='sequence folder')
    ar53_parser.add_argument('-x',   required=True,                          help='file extension')
    ar53_parser.add_argument('-g',   required=False, default=None,           help='interested genomes, no file extension, one id per line')
    ar53_parser.add_argument('-o',   required=True,                          help='output folder')
    ar53_parser.add_argument('-f',   required=False, action="store_true",    help='force overwrite')
    args = vars(ar53_parser.parse_args())
    ar53(args)

