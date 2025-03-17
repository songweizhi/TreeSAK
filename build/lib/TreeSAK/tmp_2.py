from Bio import SeqIO


def filter_by_gap(file_in, max_gap_pct, file_out):
    file_out_handle = open(file_out, 'w')
    for each_seq in SeqIO.parse(file_in, 'fasta'):
        seq_str = str(each_seq.seq)
        gap_num = seq_str.count('-')
        gap_pct = gap_num*100 / len(seq_str)
        if gap_pct <= float(max_gap_pct):
            file_out_handle.write('>%s\n%s\n' % (each_seq.id, seq_str))
    file_out_handle.close()


file_in     = '/Users/songweizhi/Desktop/Jianwei_Maggie/mmseqs_cov0.85_iden0.35/All_PeptidaseS9-RiPPs_dna100_realRiPPS9100_addref.mmseqs.iden0.35.cov0.85.representatives.trimal.aln'
max_gap_pct = 40
file_out    = '/Users/songweizhi/Desktop/Jianwei_Maggie/mmseqs_cov0.85_iden0.35/All_PeptidaseS9-RiPPs_dna100_realRiPPS9100_addref.mmseqs.iden0.35.cov0.85.representatives.trimal.maxgap40.aln'

filter_by_gap(file_in, max_gap_pct, file_out)
