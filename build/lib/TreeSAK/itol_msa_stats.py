import os
import math
import argparse
from Bio import SeqIO


iTOL_msa_stats_usage = '''
========= iTOL_msa_stats example command =========

TreeSAK iTOL_msa_stats -i concatenated.phy.fasta

==================================================
'''


def sep_path_basename_ext(file_in):
    f_path, f_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'
    f_base, f_ext = os.path.splitext(f_name)
    f_ext = f_ext[1:]
    return f_name, f_path, f_base, f_ext


def iTOL_msa_stats(args):

    aln_file = args['i']

    _, aln_path, aln_base, _ = sep_path_basename_ext(aln_file)

    stats_txt      = '%s/%s_gap_pct.txt'         % (aln_path, aln_base)
    stats_txt_itol = '%s/%s_gap_pct_iTOL.txt'    % (aln_path, aln_base)

    max_gap_pct = 0
    stats_txt_handle = open(stats_txt, 'w')
    for each_seq in SeqIO.parse(aln_file, 'fasta'):
        seq_id = each_seq.id
        seq_seq = str(each_seq.seq)
        gap_pct = seq_seq.count('-')*100/len(seq_seq)
        gap_pct = float("{0:.2f}".format(gap_pct))
        if gap_pct > max_gap_pct:
            max_gap_pct = gap_pct
        stats_txt_handle.write('%s\t%s\n' % (seq_id, gap_pct))
    stats_txt_handle.close()

    max_scale_value = math.ceil(max_gap_pct/5) * 5
    gap_pct_itol_cmd = 'TreeSAK iTOL -SimpleBar -lv %s -scale 0-25-50-75-100 -lt Gap_Pecentage -o %s' % (stats_txt, stats_txt_itol)
    os.system(gap_pct_itol_cmd)


if __name__ == '__main__':

    iTOL_msa_stats_parser = argparse.ArgumentParser(usage=iTOL_msa_stats_usage)
    iTOL_msa_stats_parser.add_argument('-i', required=True, help='MSA file')
    args = vars(iTOL_msa_stats_parser.parse_args())
    iTOL_msa_stats(args)
