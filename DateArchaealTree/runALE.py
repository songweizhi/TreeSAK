

def runALE(species_tree_file, gene_tree_file):

    gene_tree_ale_file  = '%s.ale' % gene_tree_file
    uml_rec_file        = ''
    uTs_file            = ''

    # obtain the .ale file which contains the CCPs (conditional clade probabilities)
    obtain_ale_file_cmd = 'ALEobserve %s' % gene_tree_file

    # run the reconciliation
    reconciliation_cmd = 'ALEml_undated %s %s' % (species_tree_file, gene_tree_ale_file)

    # parse ALE outputs
    run_ale_splitter_cmd = 'python ale_splitter.py -i S_S_COG3397.ufboot.ale.uml_rec -sftr'
    run_ale_parser_cmd   = 'python ale_parser.py -i FolderWithReconciliations -sft'


species_tree_file = ''
gene_tree_file    = ''

runALE(species_tree_file, gene_tree_file)

