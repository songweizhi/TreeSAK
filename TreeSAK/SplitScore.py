
SplitScore_usage = '''
============================================= SplitScore example commands =============================================

# SplitScore modules
TreeSAK SplitScore1     ->  Step 1: Infer gene tree
TreeSAK SplitScore1OMA  ->  Step 1: Infer gene tree (based on OMA outputs)
TreeSAK SplitScore2     ->  Step 2: Calculate split score

# SplitScore1
TreeSAK SplitScore1 -i OrthologousGroups.txt -s OrthologousGroupsFasta -o step1_op_dir -t 6 -f
TreeSAK SplitScore1 -i OrthologousGroups.txt -s OrthologousGroupsFasta -o step1_op_dir -t 6 -f -u interested_gnm.txt

# SplitScore2
# Please ensure that all the commands produced in step one have been executed before proceeding to step two.
TreeSAK SplitScore2 -i step1_op_dir -g combined_374_genomes.clusters.tsv -d ar53_metadata_r214.tsv -k combined_374_genomes.GTDB.r214.ar53.summary.tsv -f -t 10 -o step_2_op_dir

=======================================================================================================================
'''
