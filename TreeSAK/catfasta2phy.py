import os
import glob
from Bio import SeqIO
from Bio import AlignIO


def catfasta2phy(msa_dir, msa_ext, concatenated_msa_phy, partition_file):

    concatenated_msa_fasta = '%s.fasta' % concatenated_msa_phy
    msa_file_re            = '%s/*.%s'  % (msa_dir, msa_ext)
    msa_file_list          = [os.path.basename(file_name) for file_name in glob.glob(msa_file_re)]
    msa_file_list_sorted   = sorted(msa_file_list)

    complete_gnm_set = set()
    for each_msa_file in msa_file_list:
        pwd_msa = '%s/%s' % (msa_dir, each_msa_file)
        for each_seq in SeqIO.parse(pwd_msa, 'fasta'):
            complete_gnm_set.add(each_seq.id)

    complete_gnm_list_sorted = sorted([i for i in complete_gnm_set])

    # initialize concatenated msa dict
    gnm_to_seq_dict = {i: '' for i in complete_gnm_list_sorted}
    msa_len_dict = dict()
    for each_msa_file in msa_file_list_sorted:
        gene_id = each_msa_file.split('.' + msa_ext)[0]

        # read in msa
        current_msa_len = 0
        current_msa_len_set = set()
        pwd_current_msa = '%s/%s' % (msa_dir, each_msa_file)
        current_msa_seq_dict = dict()
        for each_seq in SeqIO.parse(pwd_current_msa, 'fasta'):
            complete_gnm_set.add(each_seq.id)
            current_msa_seq_dict[each_seq.id] = str(each_seq.seq)
            current_msa_len_set.add(len(each_seq.seq))
            current_msa_len = len(each_seq.seq)

        if len(current_msa_len_set) != 1:
            print('Sequences with different length were found in %s, program exited!' % each_msa_file)
            exit()

        msa_len_dict[gene_id] = current_msa_len

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


msa_dir                = '/Users/songweizhi/Desktop/s06_identified_marker_aln_trimmed'
msa_ext                = 'aln'
concatenated_msa_phy   = '/Users/songweizhi/Desktop/s06_identified_marker_aln_trimmed_concatenated.phy'
partition_file         = '/Users/songweizhi/Desktop/s06_identified_marker_aln_trimmed_concatenated_partition.txt'
# catfasta2phy(msa_dir, msa_ext, concatenated_msa_phy, partition_file)



msa_file = '/Users/songweizhi/Desktop/PA_75_DeltaLL_75_concatenated.phy'
msa_file_subset = '/Users/songweizhi/Desktop/PA_75_DeltaLL_75_concatenated_subset.phy'

from Bio import AlignIO

def slice_msa_by_col(msa_in, range_str, msa_out):
    alignment = AlignIO.read(msa_in, 'phylip-relaxed')

    range_l = int(range_str.split('-')[0]) - 1
    range_r = int(range_str.split('-')[1])

    aln_subset = alignment[:, range_l:range_r]
    AlignIO.write(aln_subset, msa_out, 'phylip-relaxed')


def slice_msa_by_col_manual(msa_in, range_str, msa_out):
    alignment = AlignIO.read(msa_in, 'phylip-relaxed')

    range_l = int(range_str.split('-')[0]) - 1
    range_r = int(range_str.split('-')[1])
    aln_subset = alignment[:, range_l:range_r]

    max_seq_id_len = 0
    for each_seq in aln_subset:
        seq_id_len = len(each_seq.id)
        if seq_id_len > max_seq_id_len:
            max_seq_id_len = seq_id_len
    print(max_seq_id_len)

    with open(msa_out, 'w') as msa_out_handle:
        msa_out_handle.write('%s %s\n' % (len(aln_subset), aln_subset.get_alignment_length()))
        for each_seq in aln_subset:
            seq_id = each_seq.id
            seq_id_with_space = '%s%s' % (seq_id, ' '*(max_seq_id_len + 2 - len(seq_id)))
            print(seq_id_with_space)
            msa_out_handle.write('%s%s\n' % (seq_id_with_space, str(each_seq.seq)))


#    AlignIO.write(aln_subset, msa_out, 'phylip-relaxed')


slice_range = ['1-500', '501-1000', '1001-1500', '1501-2000', '2001-2500', '2501-3000', '3001-3500', '3501-4000', '4001-4500', '4501-4879']

for each_range in slice_range:
    pwd_msa_op = '/Users/songweizhi/Desktop/%s.phy' % each_range
    slice_msa_by_col_manual(msa_file, each_range, pwd_msa_op)


def fa2phy(fasta_in, phy_out):
    alignment = AlignIO.read(fasta_in, 'fasta')
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
