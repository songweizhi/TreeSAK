import argparse
from Bio import SeqIO
from Bio import AlignIO


phy2fa_usage = '''
======= phy2fa example commands =======

TreeSAK phy2fa -i msa.phy -o msa.fa

=======================================
'''


def phy2fa(args):

    phylip_in   = args['i']
    fasta_out   = args['o']

    for aln in AlignIO.parse(phylip_in, "phylip"):
        print(aln)

    # alignments = list(AlignIO.parse(phylip_in, "phylip"))
    # print(alignments)
    # records = SeqIO.parse(phylip_in, "phylip")
    # count = SeqIO.write(records, fasta_out, "fasta")
    # print("Converted %i records" % count)


if __name__ == '__main__':

    # initialize the options parser
    phy2fa_parser = argparse.ArgumentParser()
    phy2fa_parser.add_argument('-i',      required=True,   help='input MSA in phylip format')
    phy2fa_parser.add_argument('-o',      required=True,   help='output MSA in fasta format')
    args = vars(phy2fa_parser.parse_args())
    phy2fa(args)
