
gtdb_r207_ar_meta = '/Users/songweizhi/DB/GTDB_r207/ar53_metadata_r207.tsv'


ar_p_set = set()
ar_c_set = set()
ar_o_set = set()
ar_f_set = set()
ar_g_set = set()
for each_gnm in open(gtdb_r207_ar_meta):
    if not each_gnm.startswith('accession\t'):
        each_gnm_split = each_gnm.strip().split('\t')
        gnm_tax = each_gnm_split[16]
        gnm_tax_split = gnm_tax.split(';')
        for each_rank in gnm_tax_split:
            if each_rank.startswith('p__'):
                ar_p_set.add(each_rank)
            if each_rank.startswith('c__'):
                ar_c_set.add(each_rank)
            if each_rank.startswith('o__'):
                ar_o_set.add(each_rank)
            if each_rank.startswith('f__'):
                ar_f_set.add(each_rank)
            if each_rank.startswith('g__'):
                ar_g_set.add(each_rank)

print('GTDB r207 phyla\t%s' % len(ar_p_set))
print('GTDB r207 classes\t%s' % len(ar_c_set))
print('GTDB r207 orders\t%s' % len(ar_o_set))
print('GTDB r207 families\t%s' % len(ar_f_set))
print('GTDB r207 genera\t%s' % len(ar_g_set))


for each_p in sorted([i for i in ar_p_set]):
    print(each_p)


