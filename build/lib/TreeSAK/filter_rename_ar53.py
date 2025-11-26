import os
import glob
import argparse
from Bio import SeqIO


filter_rename_ar53_usage = '''
===================================== filter_rename_ar53 example commands =====================================

TreeSAK filter_rename_ar53 -i seq_dir -x fa -g interested_gnms.txt -m interested_marker.txt -o seq_dir_renamed

===============================================================================================================
'''

def sep_path_basename_ext(file_in):

    f_path, f_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'
    f_base, f_ext = os.path.splitext(f_name)
    f_ext = f_ext[1:]

    return f_name, f_path, f_base, f_ext


def get_shared_uniq_elements(list_1, list_2):
    shared_set = set(list_1).intersection(list_2)
    list_1_uniq = []
    for e1 in list_1:
        if e1 not in shared_set:
            list_1_uniq.append(e1)
    list_2_uniq = []
    for e2 in list_2:
        if e2 not in shared_set:
            list_2_uniq.append(e2)
    return shared_set, list_1_uniq, list_2_uniq


def filter_rename_ar53(args):

    marker_seq_dir          = args['i']
    marker_seq_ext          = args['x']
    interested_gnm_txt      = args['g']
    interested_marker_txt   = args['m']
    op_dir                  = args['o']
    force_overwrite         = args['f']

    # create output folder
    if os.path.isdir(op_dir) is True:
        if force_overwrite is True:
            os.system('rm -r %s' % op_dir)
        else:
            print('%s exist, program exited!' % op_dir)
            exit()
    os.mkdir(op_dir)

    interested_marker_set = set()
    if os.path.isfile(interested_marker_txt) is True:
        for each_marker in open(interested_marker_txt):
            interested_marker_set.add(each_marker.strip())
        if len(interested_marker_set) == 0:
            print('No marker provided in %s, program exited!' % interested_gnm_txt)
            exit()

    interested_gnm_set = set()
    if os.path.isfile(interested_gnm_txt) is True:
        for each_gnm in open(interested_gnm_txt):
            interested_gnm_set.add(each_gnm.strip())
        if len(interested_gnm_set) == 0:
            print('No genome provided in %s, program exited!' % interested_gnm_txt)
            exit()

    marker_seq_re   = '%s/*.%s' % (marker_seq_dir, marker_seq_ext)
    marker_seq_list = [os.path.basename(i) for i in glob.glob(marker_seq_re)]
    marker_seq_id_list = [('.'.join(i.split('.')[:-1])) for i in marker_seq_list]
    if len(marker_seq_list) == 0:
        print('No file found in %s, program exited!' % marker_seq_dir)
        exit()

    marker_to_process = set()
    if len(interested_marker_set) == 0:
        marker_to_process = marker_seq_id_list
    else:
        shared_set, marker_seq_uniq, interested_marker_uniq = get_shared_uniq_elements(marker_seq_id_list, interested_marker_set)
        if len(interested_marker_uniq) > 0:
            print('Sequences for the following interested markers were not found:')
            print(','.join(interested_marker_uniq))
        marker_to_process = shared_set

    for marker_id in marker_to_process:
        pwd_file_in  = '%s/%s.%s' % (marker_seq_dir, marker_id, marker_seq_ext)
        pwd_file_out = '%s/%s.%s' % (op_dir, marker_id, marker_seq_ext)
        pwd_op_file_handle = open(pwd_file_out, 'w')
        for each_seq in SeqIO.parse(pwd_file_in, 'fasta'):
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

    filter_rename_ar53_parser = argparse.ArgumentParser()
    filter_rename_ar53_parser.add_argument('-i',   required=True,                          help='sequence folder')
    filter_rename_ar53_parser.add_argument('-x',   required=True,                          help='file extension')
    filter_rename_ar53_parser.add_argument('-g',   required=False, default=None,           help='interested genome, no ext, one id per line')
    filter_rename_ar53_parser.add_argument('-m',   required=False, default=None,           help='interested marker, no ext, one id per line')
    filter_rename_ar53_parser.add_argument('-o',   required=True,                          help='output folder')
    filter_rename_ar53_parser.add_argument('-f',   required=False, action="store_true",    help='force overwrite')
    args = vars(filter_rename_ar53_parser.parse_args())
    filter_rename_ar53(args)
