from itolapi import Itol


def itol_tree(tree_file, annotation_file_list, project_name, APIkey, display_mode, op_plot):

    # https://github.com/albertyw/itolapi
    # http://itol.embl.de/help.cgi#batch

    op_plot_ext = op_plot.split('.')[-1]

    # upload tree to iTOL
    itol_uploader = Itol()
    itol_uploader.params['projectName'] = project_name  # better to create a project with a unique name.
    itol_uploader.params['APIkey'] = APIkey             # sine we are using the same account, we can use the same APIkey
    itol_uploader.params['treeName'] = tree_file
    itol_uploader.add_file(tree_file)

    # upload annotation files to iTOL
    for annotation_file in annotation_file_list:
        itol_uploader.add_file(annotation_file)

    status = itol_uploader.upload()
    # import pdb;pdb.set_trace()
    assert status != False

    # the following parameters are optional, refer to https://itol.embl.de/help.cgi#batchExp
    if len(annotation_file_list) == 1:
        datasets_visible_str = '0'
    elif len(annotation_file_list) == 2:
        datasets_visible_str = '0,1'
    elif len(annotation_file_list) == 3:
        datasets_visible_str = '0,1,2'
    else:
        datasets_visible_str = ','.join([str(i) for i in list(range(0, len(annotation_file_list)))])

    itol_exporter = itol_uploader.get_itol_export()
    itol_exporter.set_export_param_value('datasets_visible', datasets_visible_str)
    itol_exporter.set_export_param_value('display_mode', display_mode)
    itol_exporter.set_export_param_value('range_mode', '2')
    itol_exporter.set_export_param_value('dashed_lines', '1')
    # itol_exporter.set_export_param_value('current_font_size', '96')
    itol_exporter.set_export_param_value('line_width', '3')
    itol_exporter.set_export_param_value('vertical_shift_factor', '0.9')
    itol_exporter.set_export_param_value('horizontal_scale_factor', '1')
    itol_exporter.set_export_param_value('format', op_plot_ext)
    itol_exporter.export(op_plot)


tree_file            = '/Users/songweizhi/Desktop/iTOL_tree/tree_file.tree'
tree_plot            = '/Users/songweizhi/Desktop/iTOL_tree/tree_file.tree.pdf'
API_key              = 'S1kZZuDHc0d5M7J5vLnUNQ'
project_name         = 'batch_access_tmp'
display_mode         = '1'  # 1=rectangular, 2=circular, 3=unrooted
annotation_file_list = ['/Users/songweizhi/Desktop/iTOL_tree/Taxon.txt',
                        '/Users/songweizhi/Desktop/iTOL_tree/LifeStyle.txt',
                        '/Users/songweizhi/Desktop/iTOL_tree/Abundance.txt',
                        '/Users/songweizhi/Desktop/iTOL_tree/MAG_Size.txt']

itol_tree(tree_file, annotation_file_list, project_name, API_key, display_mode, tree_plot)

