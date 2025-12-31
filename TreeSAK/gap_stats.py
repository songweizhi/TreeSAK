import argparse
from Bio import SeqIO


gap_stats_usage = '''
================================== gap_stats example commands ==================================

TreeSAK gap_stats -i msa.fasta -c 40 -o1 stats.txt -o2 msa_maxgap40.fasta -o3 MSA_gap_iTOL.txt

================================================================================================
'''


def gap_stats(args):

    msa_in_fa    = args['i']
    gap_cutoff   = args['c']
    op_stats_txt = args['o1']
    op_msa_file  = args['o2']
    op_itol_file = args['o3']

    if gap_cutoff != 100:
        op_msa_file_handle = open(op_msa_file, 'w')

    qualified_seq_set = set()
    gap_pct_dict = dict()
    for each_seq in SeqIO.parse(msa_in_fa, 'fasta'):
        seq_id = each_seq.id
        seq_str = str(each_seq.seq)
        gap_pct = seq_str.count('-')*100/len(seq_str)
        gap_pct = float("{0:.2f}".format(gap_pct))
        gap_pct_dict[seq_id] = gap_pct

        if (gap_cutoff != 100) and (gap_pct <= gap_cutoff):
            op_msa_file_handle.write('>%s\n' % (seq_id))
            op_msa_file_handle.write('%s\n' % (seq_str))
            qualified_seq_set.add(seq_id)
    if gap_cutoff != 100:
        op_msa_file_handle.close()

    gap_pct_sorted = sorted(gap_pct_dict.items(), key=lambda x: x[1])

    # write stats
    op_stats_txt_handle = open(op_stats_txt, 'w')
    op_stats_txt_handle.write('Sequence\tGap\n')
    for each_seq in gap_pct_sorted:
        op_stats_txt_handle.write('%s\t%s\n' % (each_seq[0], each_seq[1]))
    op_stats_txt_handle.close()

    # prepare itol file
    if op_itol_file != None:
        op_itol_file_handle = open(op_itol_file, 'w')
        op_itol_file_handle.write('DATASET_BINARY\n\nSEPARATOR TAB\nDATASET_LABEL\tMSA_Gap_%s\nCOLOR\tred\n' % gap_cutoff)
        op_itol_file_handle.write('SHOW_LABELS\t1\nLABEL_ROTATION\t45\nLABEL_SHIFT\t5\n')
        op_itol_file_handle.write('FIELD_LABELS\tMSA_Gap_%s\n' % gap_cutoff)
        op_itol_file_handle.write('FIELD_COLORS\tred\n')
        op_itol_file_handle.write('FIELD_SHAPES\t2\n')
        op_itol_file_handle.write('MARGIN\t10\n')
        op_itol_file_handle.write('HORIZONTAL_GRID\t0\n')
        op_itol_file_handle.write('VERTICAL_GRID\t0\n')
        op_itol_file_handle.write('\nDATA\n')
        for each in sorted(list(qualified_seq_set)):
            op_itol_file_handle.write('%s\t1\n' % each)
        op_itol_file_handle.close()


if __name__ == '__main__':

    gap_stats_parser = argparse.ArgumentParser()
    gap_stats_parser.add_argument('-i',  required=True,                             help='MSA in fasta format')
    gap_stats_parser.add_argument('-c',  required=False, default=100, type=float,   help='maximum gap allowed in MSAs, default is 100, i.e., no filtering')
    gap_stats_parser.add_argument('-o1', required=True,                             help='output stats file')
    gap_stats_parser.add_argument('-o2', required=False, default=None,              help='filtered MSA')
    gap_stats_parser.add_argument('-o3', required=False, default=None,              help='iTOL file')
    args = vars(gap_stats_parser.parse_args())
    gap_stats(args)
