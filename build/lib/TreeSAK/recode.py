import argparse
from Bio import SeqIO


recode_usage = '''
============================ recode example commands ============================

TreeSAK recode -i msa.fa -m s4 -o recoded_msa_SR4.fa
TreeSAK recode -i msa.fa -m d4 -o recoded_msa_Dayhoff4.fa
TreeSAK recode -i msa.fa -m d6 -o recoded_msa_Dayhoff6.fa

Note:
This script is modified based on the Recode_aa.py from Anja Spang.
It was used to recode AA to SR4 (s4), Dayhoff4 (d4) or Dayhoff6 (d6) categories.  
Please refer to https://doi.org/10.1038/s41467-020-17408-w for details.

Recoding schemes
1. Selenocysteine will be recoded to a gap.
2. s4: A,G,N,P,S,T = A; C,H,W,Y = C; D,E,K,Q,R = G; F,I,L,M,V = T
3. d4: A,G,P,S,T = A; D,E,N,Q = D; H,K,R = H; F,Y,W,I,L,M,V = F; C = ?
4. d6: A,G,P,S,T = A; D,E,N,Q = D; H,K,R = H; F,Y,W = F; I,L,M,V = I; C = C

=================================================================================
'''


def recode(args):

    msa_in   = args['i']
    msa_out  = args['o']
    category = args['m']

    DH6 = {'-': '-', 'A': 'A', 'G': 'A', 'P': 'A', 'S': 'A', 'T': 'A', 'D': 'D', 'E': 'D', 'N': 'D', 'Q': 'D', 'H': 'H', 'K': 'H', 'R': 'H', 'F': 'F', 'Y': 'F', 'W': 'F', 'I': 'I', 'L': 'I', 'M': 'I', 'V': 'I', 'C': 'C'}
    DH4 = {'-': '-', 'A': 'A', 'G': 'A', 'P': 'A', 'S': 'A', 'T': 'A', 'D': 'D', 'E': 'D', 'N': 'D', 'Q': 'D', 'H': 'H', 'K': 'H', 'R': 'H', 'F': 'F', 'Y': 'F', 'W': 'F', 'I': 'F', 'L': 'F', 'M': 'F', 'V': 'F', 'C': '-'}
    SR4 = {'-': '-', 'A': 'A', 'G': 'A', 'N': 'A', 'P': 'A', 'S': 'A', 'T': 'A', 'C': 'C', 'H': 'C', 'W': 'C', 'Y': 'C', 'D': 'G', 'E': 'G', 'K': 'G', 'Q': 'G', 'R': 'G', 'F': 'T', 'I': 'T', 'L': 'T', 'M': 'T', 'V': 'T'}

    msa_out_handle = open(msa_out, 'w')
    for seq_record in SeqIO.parse(msa_in, "fasta"):
        header = str(seq_record.description).strip()
        seq = str(seq_record.seq).strip()
        new_seq = ''
        for item in seq:
            if category in ['D6', 'd6']:
                if item in DH6:
                    new_seq = new_seq + str(DH6.get(item))
                else:
                    new_seq = new_seq + str('-')
            elif category in ['D4', 'd4']:
                if item in DH4:
                    new_seq = new_seq + str(DH4.get(item))
                else:
                    new_seq = new_seq + str('-')
            elif category in ['S4', 's4']:
                if item in SR4:
                    new_seq = new_seq + str(SR4.get(item))
                else:
                    new_seq = new_seq + str('-')
            else:
                print('Please choose recoding scheme from d4, d6 and s4, program exited!')
                exit()

        msa_out_handle.write(">%s\n%s\n" % (header, new_seq))
    msa_out_handle.close()


if __name__ == '__main__':

    recode_parser = argparse.ArgumentParser(description='Recode amino acids to Dayoff 4, Dayoff 6 or SR4 categories')
    recode_parser.add_argument('-i', required=True,  help='input file')
    recode_parser.add_argument('-m', required=True,  help='recoding scheme, choose from d4, d6 or s4')
    recode_parser.add_argument('-o', required=True,  help='output file')
    args = vars(recode_parser.parse_args())
    recode(args)
