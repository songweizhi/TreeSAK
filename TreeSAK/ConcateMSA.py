import os
import glob
import argparse
from Bio import SeqIO
from Bio import AlignIO


ConcateMSA_usage = '''
================= ConcateMSA example commands =================

TreeSAK ConcateMSA -i aln -x aln -p concatenated -gene2gnm

# output file include:
concatenated.fasta
concatenated.phylip
concatenated.partition.txt

===============================================================
'''


def ConcateMSA(args):

    msa_dir                 = args['i']
    msa_ext                 = args['x']
    op_prefix               = args['p']
    gene2gnm                = args['gene2gnm']

    concatenated_msa_phy    = '%s.phylip'           % op_prefix
    concatenated_msa_fasta  = '%s.fasta'            % op_prefix
    partition_file          = '%s.partition.txt'    % op_prefix

    msa_file_re             = '%s/*.%s'             % (msa_dir, msa_ext)
    msa_file_list           = [os.path.basename(file_name) for file_name in glob.glob(msa_file_re)]
    msa_file_list_sorted    = sorted(msa_file_list)

    complete_gnm_set = set()
    for each_msa_file in msa_file_list:
        pwd_msa = '%s/%s' % (msa_dir, each_msa_file)
        for each_seq in SeqIO.parse(pwd_msa, 'fasta'):
            seq_id = each_seq.id
            if gene2gnm is True:
                seq_id = '_'.join(seq_id.split('_')[:-1])
            complete_gnm_set.add(seq_id)

    complete_gnm_list_sorted = sorted([i for i in complete_gnm_set])

    # initialize concatenated msa dict
    gnm_to_seq_dict = {i: '' for i in complete_gnm_list_sorted}
    msa_len_dict = dict()
    for each_msa_file in msa_file_list_sorted:
        msa_id = each_msa_file.split('.' + msa_ext)[0]

        # read in msa
        current_msa_len = 0
        current_msa_len_set = set()
        pwd_current_msa = '%s/%s' % (msa_dir, each_msa_file)
        current_msa_seq_dict = dict()
        for each_seq in SeqIO.parse(pwd_current_msa, 'fasta'):
            seq_id = each_seq.id
            if gene2gnm is True:
                seq_id = '_'.join(seq_id.split('_')[:-1])
            complete_gnm_set.add(seq_id)
            current_msa_seq_dict[seq_id] = str(each_seq.seq)
            current_msa_len_set.add(len(each_seq.seq))
            current_msa_len = len(each_seq.seq)

        if len(current_msa_len_set) != 1:
            print('Sequences with different length were found in %s, program exited!' % each_msa_file)
            exit()

        msa_len_dict[msa_id] = current_msa_len

        # add sequence to concatenated msa dict
        for each_gnm in complete_gnm_list_sorted:
            msa_seq = current_msa_seq_dict.get(each_gnm, current_msa_len*'-')
            gnm_to_seq_dict[each_gnm] += msa_seq

    # write out concatenated msa
    concatenated_msa_handle = open(concatenated_msa_fasta, 'w')
    for each_gnm in complete_gnm_list_sorted:
        concatenated_msa_handle.write('>%s\n' % each_gnm)
        concatenated_msa_handle.write('%s\n' % gnm_to_seq_dict[each_gnm])
    concatenated_msa_handle.close()

    # write out partition file
    end_pos = 0
    partition_file_handle = open(partition_file, 'w')
    for each_m in msa_file_list_sorted:
        gene_id = each_m.split('.' + msa_ext)[0]
        current_m_len = msa_len_dict[gene_id]
        partition_file_handle.write('%s = %s-%s\n' % (each_m, (end_pos + 1), (end_pos + current_m_len)))
        end_pos += current_m_len
    partition_file_handle.close()

    # convert msa in fasta to phy
    AlignIO.convert(concatenated_msa_fasta, 'fasta', concatenated_msa_phy, 'phylip-relaxed')



if __name__ == '__main__':

    # initialize the options parser
    parser = argparse.ArgumentParser()
    parser.add_argument('-i',           required=True,                       help='input MSA folder')
    parser.add_argument('-x',           required=False, default='aln',       help='input file extension')
    parser.add_argument('-p',           required=True,                       help='output prefix')
    parser.add_argument('-gene2gnm',    required=False, action="store_true", help='gene id to gnm id, split sequence id before the last _')

    args = vars(parser.parse_args())
    ConcateMSA(args)
