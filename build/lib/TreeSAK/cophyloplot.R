suppressMessages(suppressWarnings(library(optparse)))
suppressMessages(suppressWarnings(library(ape)))

################################ argument parser ###############################

option_list = list(
  make_option(c("-i", "--tree1"),   type="character", default=NULL, help="tree file 1"),
  make_option(c("-s", "--tree2"),   type="character", default=NULL, help="tree file 2"),
  make_option(c("-l", "--link"),    type="character", default=NULL, help="connection file"),
  make_option(c("-o", "--output"),  type="character", default=NULL, help="output plot"));
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

tree_file_1     = opt$tree1
tree_file_2     = opt$tree2
link_table_file = opt$link
op_plot         = opt$output

################################################################################

# tree_file_1     = '/Users/songweizhi/Desktop/Sponge_r226/08_Host_specificity/codiv_data_test/RefSeqs_with_AOA_18S_iden99_g_representatives_with_JL_rooted.treefile'
# tree_file_2     = '/Users/songweizhi/Desktop/Sponge_r226/08_Host_specificity/codiv_data_test/AOA_tree_654_with_fixed_topo1_LG_rooted.treefile'
# link_table_file = '/Users/songweizhi/Desktop/Sponge_r226/08_Host_specificity/codiv_data_test/codiv_sponge_to_AOA_op.txt'

################################################################################

tree_1         <- read.tree(tree_file_1)
tree_2         <- read.tree(tree_file_2)
association_df <- read.table(link_table_file, header = FALSE)
association    <- cbind(association_df$V1, association_df$V2)
link_color     <- association_df$V3

pdf(op_plot, width =10,  height =10)

plot(tree_1, type = "c")
plot(tree_2, type = "c")

cophyloplot(tree_1, tree_2, assoc = association,                # data
            show.tip.label = FALSE, rotate = F,                 # layout
            length.line = 1, space = 200, gap = 3,              # layout
            col=link_color, lwd=1, lty = 1)                     # linking lines

invisible(dev.off())
rm(list=ls())
