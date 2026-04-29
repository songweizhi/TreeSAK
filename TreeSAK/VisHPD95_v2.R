library(ggplot2)
library(optparse)


plot_grouped_HPD95 <- function(data_file, plot_width, plot_height, plot_file, x_axis_label_order, break_point_str, v_line_str){
 
  v_line_list      = as.numeric(unlist(strsplit(v_line_str, split = ",")))
  break_point_list = as.numeric(unlist(strsplit(break_point_str, split = ",")))
  label_list       = unlist(strsplit(x_axis_label_order, split = ","))
  
  dat <- read.table(data_file, header = T)
  
  ggplot(dat, aes(x=Post_X, y = Mean, ymin = Low, ymax = High)) +
    geom_linerange(aes(color = Color), linewidth = 0.6, alpha = 0.5) +
    geom_point(aes(color = Color, shape = Shape))+
    scale_x_continuous("",breaks=break_point_list, labels=label_list) + # here you define coordinates for A and B 

  # 
  # ggplot(dat, aes(x=Node, # x=factor(Post_X, level=unlist(strsplit(x_axis_label_order, split = ",")))
  #                 y = Mean, ymin = Low, ymax = High)) +
  #   geom_pointrange(aes(color=factor(Color), shape=factor(Shape)),
  #                   position=position_dodge(width=0.6),                         # controls distance between groups
  #                   linewidth = 0.6,                                            # line width
  #                   size=0.5) +                                                # size of shape
    theme_bw() +                                                                # remove background
    theme(panel.grid.major=element_blank(),                                     # remove grid
          panel.grid.minor=element_blank()) +                                  # remove grid
    xlab("") +                                                                  # x-axis label text
    ylab("95% HPD CI") +  
    geom_vline(xintercept = v_line_list, linewidth = 0.1, linetype = "dashed") + # y-axis label text
    theme(axis.text.x=element_text(size=12, color='black', angle=315, hjust=0),  # x-axis label, rotate at an angle of 45
          axis.text.y=element_text(size=12, color='black'),                     # y-axis label
          legend.text=element_text(size=10))                                 # legend label
    # scale_color_discrete(name="Color") +                                        # customize color legend, title
    # guides(color=guide_legend(override.aes=list(linetype=0))) +                 # customize color legend
    # scale_shape_discrete(name="Shape") +                                        # customize color legend, title
    # guides(shape=guide_legend(override.aes=list(linetype=0, color='grey')))     # customize color legend,

  # write to file
  ggsave(plot_file, width=plot_width, height=plot_height, dpi=300)
}


option_list = list(
  make_option(c("-i", "--datain"),      type="character", default=NULL,    help="input data matrix"),
  make_option(c("-x", "--width"),       type="double",    default=10,      help="plot width"),
  make_option(c("-y", "--height"),      type="double",    default=6,       help="plot height"),
  make_option(c("-l", "--labelorder"),  type="character", default=NULL,    help="X-axis label order"),
  make_option(c("-b", "--breakpoint"),  type="character", default=NULL,    help="X-axis break point"),
  make_option(c("-v", "--vlines"),      type="character", default=NULL,    help="vertical lines"),
  make_option(c("-o", "--plotout"),     type="character", default=NULL,    help="output plot"));

opt_parser         = OptionParser(option_list=option_list);
opt                = parse_args(opt_parser);
data_matrix_txt    = opt$datain
plot_width         = opt$width
plot_height        = opt$height
output_plot        = opt$plotout
x_axis_label_order = opt$labelorder
break_point_str    = opt$breakpoint
v_line_str         = opt$vlines

plot_grouped_HPD95(data_matrix_txt, plot_width, plot_height, output_plot, x_axis_label_order, break_point_str, v_line_str)
