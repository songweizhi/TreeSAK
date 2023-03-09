import os
import argparse
from Bio import SeqIO


get_arCOG_seq_usage = '''
=========================== get_arCOG_seq example commands ===========================

TreeSAK get_arCOG_seq -id cog_id.txt -db_dir /Users/songweizhi/DB/arCOG18 -o op_dir

# required db files
ar18.ar14.02.csv, arCOG_names_220807.txt and ar18.fa

======================================================================================
'''


def get_arCOG_seq(args):

    cog_id_txt          = args['i']
    db_dir              = args['db_dir']
    op_dir              = args['o']
    force_create_dir    = args['f']

    ar18_ar14_02_csv = '%s/ar18.ar14.02.csv'        % db_dir
    cog_des_txt      = '%s/arCOG_names_220807.txt'  % db_dir
    ar18_fa          = '%s/ar18.fa'                 % db_dir
    cog_metadata_txt = '%s/metadata.txt'            % op_dir


    if os.path.isdir(op_dir) is True:
        if force_create_dir is True:
            os.system('rm -r %s' % op_dir)
        else:
            print('Output folder already exist, program exited!')
            exit()
    os.system('mkdir %s' % op_dir)

    cog_des_dict = dict()
    for each_cog in open(cog_des_txt, encoding="ISO-8859-1"):
        each_cog_split = each_cog.strip().split('\t')
        cog_des_dict[each_cog_split[0]] = each_cog_split[1:]

    cog_id_set = set()
    for each_cog in open(cog_id_txt):
        cog_id_set.add(each_cog.strip().replace('ArCOG', 'arCOG'))

    seq_id_set = set()
    seq_to_arcog_dict = dict()
    arcog_to_seq_dict = dict()
    for each_line in open(ar18_ar14_02_csv):
        each_line_split = each_line.strip().split(',')
        arcog_id = each_line_split[6]
        seq_id = each_line_split[2]

        if arcog_id in cog_id_set:
            seq_id_set.add(seq_id)
            if arcog_id not in arcog_to_seq_dict:
                arcog_to_seq_dict[arcog_id] = {seq_id}
            else:
                arcog_to_seq_dict[arcog_id].add(seq_id)

            if seq_id not in seq_to_arcog_dict:
                seq_to_arcog_dict[seq_id] = {arcog_id}
            else:
                seq_to_arcog_dict[seq_id].add(arcog_id)

    # write out sequence by arCOG
    for each_seq in SeqIO.parse(ar18_fa, 'fasta'):
        seq_id = each_seq.id
        if seq_id in seq_id_set:
            seq_cog_set = seq_to_arcog_dict.get(seq_id, [])
            seq_cog_list = [i for i in seq_cog_set]
            if len(seq_cog_list) == 1:
                pwd_fa = '%s/%s.fa' % (op_dir, seq_cog_list[0])
                with open(pwd_fa, 'a') as pwd_fa_handle:
                    pwd_fa_handle.write('>%s\n' % seq_id)
                    pwd_fa_handle.write('%s\n' % str(each_seq.seq))

    # write out metadata
    cog_metadata_txt_handle = open(cog_metadata_txt, 'w')
    for each_c in sorted([i for i in cog_id_set]):
        each_c_desc = '\t'.join(cog_des_dict[each_c])
        cog_metadata_txt_handle.write('%s\t%s\n' % (each_c, each_c_desc))
    cog_metadata_txt_handle.close()


if __name__ == '__main__':

    # initialize the options parser
    parser = argparse.ArgumentParser()
    parser.add_argument('-i',      required=True,                       help='arCOD id file, one id per line')
    parser.add_argument('-db_dir', required=True,                       help='database folder')
    parser.add_argument('-o',      required=True,                       help='output folder')
    parser.add_argument('-f',      required=False, action="store_true", help='force overwrite existing output folder')
    args = vars(parser.parse_args())
    get_arCOG_seq(args)
