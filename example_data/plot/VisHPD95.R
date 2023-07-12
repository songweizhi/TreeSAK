library(ggplot2)
library(optparse)


plot_grouped_HPD95 <- function(data_file, plot_width, plot_height, plot_file){
  
  dat <- read.table(data_file, header = T)
  
  ggplot(dat, aes(x = Var, y = Mean, ymin = Low, ymax = High)) +
    geom_pointrange(aes(col = factor(Test), shape=factor(Shape)),
                    position=position_dodge(width=0.6),  # controls distance between groups
                    linewidth = 0.9,   # line width
                    size=0.75) +       # size of shape
    ylab("95% HPD CI") +
    scale_color_discrete(name="Color") +
    scale_shape_discrete(name="Shape") +
    xlab("") +
    theme_bw() +
    theme(panel.grid.major=element_blank(), 
          panel.grid.minor=element_blank()) +
    theme(axis.text.x=element_text(size=12, angle=45, hjust=1),  # rotate x-axis label at an angle of 45
          axis.text.y=element_text(size=12), 
          legend.text=element_text(size=10)) +
    guides(color=guide_legend(override.aes=list(linetype=0))) +  # remove the short line in color legend 
    guides(shape=guide_legend(override.aes=list(linetype=0)))    # remove the short line in color legend 
  
  
  # write to file
  ggsave(plot_file, width=plot_width, height=plot_height, dpi=300)
}


option_list = list(
  make_option(c("-i", "--datain"),  type="character", default=NULL, help="input data matrix"),
  make_option(c("-x", "--width"),   type="double",    default=8,    help="plot width"),
  make_option(c("-y", "--height"),  type="double",    default=5,    help="plot height"),
  make_option(c("-o", "--plotout"), type="character", default=NULL, help="output plot"));

opt_parser      = OptionParser(option_list=option_list);
opt             = parse_args(opt_parser);
data_matrix_txt = opt$datain
plot_width      = opt$width
plot_height     = opt$height
output_plot     = opt$plotout

plot_grouped_HPD95(data_matrix_txt, plot_width, plot_height, output_plot)
