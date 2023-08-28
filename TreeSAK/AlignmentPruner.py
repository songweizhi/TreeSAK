import os
import argparse


AlignmentPruner_usage = '''
================== AlignmentPruner example commands ==================

TreeSAK AlignmentPruner -i in.aln -c -0.2 -o out_RmHeteSite_20.aln
TreeSAK AlignmentPruner -i in.aln -c -0.4 -o out_RmHeteSite_40.aln
TreeSAK AlignmentPruner -i in.aln -c -0.6 -o out_RmHeteSite_60.aln
TreeSAK AlignmentPruner -i in.aln -c -0.8 -o out_RmHeteSite_80.aln

======================================================================
'''


def AlignmentPruner(args):

    msa_in              = args['i']
    heter_cutoff        = args['c']
    msa_out             = args['o']

    pwd_current_file        = os.path.realpath(__file__)
    current_file_path       = '/'.join(pwd_current_file.split('/')[:-1])
    pwd_alignment_pruner_pl = '%s/alignment_pruner.pl' % current_file_path
    alignment_pruner_cmd    = 'perl %s --file %s --conserved_threshold %s > %s' % (pwd_alignment_pruner_pl, msa_in, heter_cutoff, msa_out)
    os.system(alignment_pruner_cmd)


if __name__ == '__main__':

    AlignmentPruner_parser = argparse.ArgumentParser()
    AlignmentPruner_parser.add_argument('-i',   required=True,  help='input MSA')
    AlignmentPruner_parser.add_argument('-c',   required=True,  help='heterogeneity cutoff')
    AlignmentPruner_parser.add_argument('-o',   required=True,  help='output MSA')
    args = vars(AlignmentPruner_parser.parse_args())
    AlignmentPruner(args)
