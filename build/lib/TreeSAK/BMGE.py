import os
import argparse


BMGE_usage = '''
======================= BMGE example commands =======================

# require: java

TreeSAK BMGE -p demo -i input.aln -m BLOSUM30 -esc 0.55

# Settings for calculating split score (Nina Dombrowski): 
# -t AA -m BLOSUM30 -h 0.55

=====================================================================
'''


def BMGE(args):

    op_prefix               = args['p']
    msa_in                  = args['i']
    trim_model              = args['m']
    entropy_score_cutoff    = args['esc']

    # define file name
    msa_out_phylip = '%s.BMGE.phylip' % op_prefix
    msa_out_fasta  = '%s.BMGE.fasta'  % op_prefix
    msa_out_nexus  = '%s.BMGE.nexus'  % op_prefix
    msa_out_html   = '%s.BMGE.html'   % op_prefix

    # specify path to BMGE.jar
    current_file_path   = '/'.join(os.path.realpath(__file__).split('/')[:-1])
    pwd_bmge_jar        = '%s/BMGE.jar' % current_file_path

    # run BMGE
    bmge_cmd = 'java -jar %s -i %s -m %s -t AA -h %s -op %s -of %s -on %s -oh %s' % (pwd_bmge_jar, msa_in, trim_model, entropy_score_cutoff, msa_out_phylip, msa_out_fasta, msa_out_nexus, msa_out_html)
    print('Running %s' % bmge_cmd)
    os.system(bmge_cmd)

    print('Done!')


if __name__ == '__main__':

    BMGE_parser = argparse.ArgumentParser()
    BMGE_parser.add_argument('-p',   required=True,                         help='output prefix')
    BMGE_parser.add_argument('-i',   required=True,                         help='input MSA')
    BMGE_parser.add_argument('-m',   required=False, default='BLOSUM30',    help='trim model, default: BLOSUM30')
    BMGE_parser.add_argument('-esc', required=False, default='0.55',        help='entropy score cutoff, default: 0.55')
    args = vars(BMGE_parser.parse_args())
    BMGE(args)
