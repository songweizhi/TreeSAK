import argparse
from Bio import SeqIO


gap_stats_usage = '''
======= gap_stats example commands =======

TreeSAK gap_stats -i msa.fasta

==========================================
'''


def gap_stats(args):

    msa_in_fa = args['i']

    gap_pct_dict = dict()
    for each_seq in SeqIO.parse(msa_in_fa, 'fasta'):
        seq_id = each_seq.id
        seq_str = str(each_seq.seq)
        gap_pct = seq_str.count('-')*100/len(seq_str)
        gap_pct = float("{0:.2f}".format(gap_pct))
        gap_pct_dict[seq_id] = gap_pct

    gap_pct_sorted = sorted(gap_pct_dict.items(), key=lambda x: x[1])

    print('Sequence\tGap')
    for each_seq in gap_pct_sorted:
        print('%s\t%s' % (each_seq[0], each_seq[1]))


if __name__ == '__main__':

    gap_stats_parser = argparse.ArgumentParser()
    gap_stats_parser.add_argument('-i', required=True, help='MSA in fasta format')
    args = vars(gap_stats_parser.parse_args())
    gap_stats(args)
