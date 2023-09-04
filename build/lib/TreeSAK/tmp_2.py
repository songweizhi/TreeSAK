from Bio import SeqIO

fa = '/Users/songweizhi/Desktop/Phylogenomics/archaeal_marker_gene_sets/Yang_70.fasta'

for each_seq in SeqIO.parse(fa, 'fasta'):
    seq_file_subset = '/Users/songweizhi/Desktop/Phylogenomics/archaeal_marker_gene_sets/%s.faa' % each_seq.id
    seq_file_subset_handle = open(seq_file_subset, 'w')
    seq_file_subset_handle.write('>%s\n' % each_seq.id)
    seq_file_subset_handle.write('%s\n' % str(each_seq.seq))
    seq_file_subset_handle.close()

