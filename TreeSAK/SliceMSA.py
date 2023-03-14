import os
import argparse
from Bio import AlignIO


SliceMSA_usage = '''
================== SliceMSA example commands ==================

BioSAK SliceMSA -i 16S.aln -s 200-300 -o SliceMSA_op
BioSAK SliceMSA -i 16S.aln -s sections.txt -o SliceMSA_op

# example
200-300     export columns 200-300
-100        export columns 1-300
500-        export columns from 500 to the end

===============================================================
'''


def SliceMSA(args):

    msa_in_file         = args['i']
    col_to_select_txt   = args['s']
    op_dir              = args['o']
    force_overwriting   = args['f']

    # check output folder
    if os.path.isdir(op_dir) is True:
        if force_overwriting is True:
            os.system('rm -r %s' % op_dir)
        else:
            print('Output folder already exist, program exited!')
            exit()
    os.system('mkdir %s' % op_dir)

    # read in msa
    msa_in = AlignIO.read(msa_in_file, "fasta")

    # parse provided sections
    section_to_select_list = []
    if os.path.isfile(col_to_select_txt) is False:
        col_to_select_txt_split = col_to_select_txt.strip().split('-')
        if col_to_select_txt == '-':
            section_to_select_list.append(['1', str(msa_in.get_alignment_length())])
        elif col_to_select_txt.startswith('-'):
            section_to_select_list.append(['1', col_to_select_txt_split[1]])
        elif col_to_select_txt.endswith('-'):
            section_to_select_list.append([col_to_select_txt_split[0], str(msa_in.get_alignment_length())])
        else:
            section_to_select_list.append(col_to_select_txt_split)
    else:
        for each_section in open(col_to_select_txt):
            each_section_split = each_section.strip().split('-')
            if each_section == '-':
                section_to_select_list.append(['1', str(msa_in.get_alignment_length())])
            elif each_section.startswith('-'):
                section_to_select_list.append(['1', each_section_split[1]])
            elif each_section.endswith('-'):
                section_to_select_list.append([each_section_split[0], str(msa_in.get_alignment_length())])
            else:
                section_to_select_list.append(each_section_split)

    # write out sections
    for each_section in section_to_select_list:
        pwd_op_file = '%s/%s.fa' % (op_dir, '-'.join(each_section))
        current_section = msa_in[:, (int(each_section[0])-1):(int(each_section[1]))]
        AlignIO.write(current_section, pwd_op_file, 'fasta')

    print('Done!')


if __name__ == '__main__':
    # arguments for rename_seq_parser
    SliceMSA_parser = argparse.ArgumentParser()
    SliceMSA_parser.add_argument('-i',   required=True,                         help='input MSA in fasta format')
    SliceMSA_parser.add_argument('-s',   required=True,                         help='columns to export, e.g. 200-300, -100, 50-')
    SliceMSA_parser.add_argument('-o',   required=True,                         help='output folder')
    SliceMSA_parser.add_argument('-f',   required=False, action="store_true",   help='force overwrite existing output folder')
    args = vars(SliceMSA_parser.parse_args())
    SliceMSA(args)

