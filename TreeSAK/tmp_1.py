
def gene_gain_and_loss(TableEvents_tsv):

    for each_line in open(TableEvents_tsv):
        if not each_line.startswith('Family\tBranchType\tBranch'):
            each_line_split = each_line.strip().split('\t')
            protein_family_id   = each_line_split[0]
            node_id             = each_line_split[2]
            Losses              = float(each_line_split[5])
            Originations        = float(each_line_split[6])
            Extinctinonprob     = float(each_line_split[9])

            print('%s\t%s\t%s\t%s\t%s' % (protein_family_id, node_id, Losses, '', Extinctinonprob))


TableEvents_tsv = '/Users/songweizhi/Desktop/Sponge_r220/10_dRep95_255_ALE_wd/dRep95_255_ALE3_op_dir_80/TableEvents.tsv'

gene_gain_and_loss(TableEvents_tsv)

