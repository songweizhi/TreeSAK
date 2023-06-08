import os
import argparse
from Bio import AlignIO


SingleLinePhy_usage = '''
======== SingleLinePhy example commands ========

TreeSAK SingleLinePhy -i in.phy -o out.phy

================================================
'''


def SingleLinePhy(args):

    phy_in  = args['i']
    phy_out = args['o']

    # check input file
    if os.path.isfile(phy_in) is False:
        print('input file not found, program exited!')
        exit()

    alignment = AlignIO.read(phy_in, 'phylip-relaxed')

    max_seq_id_len = 0
    for each_seq in alignment:
        seq_id_len = len(each_seq.id)
        if seq_id_len > max_seq_id_len:
            max_seq_id_len = seq_id_len

    with open(phy_out, 'w') as msa_out_handle:
        msa_out_handle.write('%s %s\n' % (len(alignment), alignment.get_alignment_length()))
        for each_seq in alignment:
            seq_id = each_seq.id
            seq_id_with_space = '%s%s' % (seq_id, ' ' * (max_seq_id_len + 2 - len(seq_id)))
            msa_out_handle.write('%s%s\n' % (seq_id_with_space, str(each_seq.seq)))

    print('Done!')


if __name__ == '__main__':

    # initialize the options parser
    parser = argparse.ArgumentParser()
    parser.add_argument('-i',   required=True,                          help='input file')
    parser.add_argument('-o',   required=True,                          help='output file')
    args = vars(parser.parse_args())
    SingleLinePhy(args)
