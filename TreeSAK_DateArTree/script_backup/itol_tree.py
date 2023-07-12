from itolapi import Itol


def itol_tree(tree_file, annotation_file_list, project_name, APIkey, display_mode, op_plot):

    # https://github.com/albertyw/itolapi
    # http://itol.embl.de/help.cgi#batch

    op_plot_ext = op_plot.split('.')[-1]

    # upload tree to iTOL
    itol_uploader = Itol()
    itol_uploader.params['projectName'] = project_name  # better to create a project with a unique name.
    itol_uploader.params['APIkey'] = APIkey  # sine we are the same account, we can use the same APIkey
    itol_uploader.params['treeName'] = tree_file
    itol_uploader.add_file(tree_file)

    # upload annotation files to iTOL
    for annotation_file in annotation_file_list:
        itol_uploader.add_file(annotation_file)

    status = itol_uploader.upload()
    #import pdb;pdb.set_trace()
    assert status != False

    # the following parameters are optional, refer to https://itol.embl.de/help.cgi#batchExp
    itol_exporter = itol_uploader.get_itol_export()
    itol_exporter.set_export_param_value('datasets_visible', '0')
    itol_exporter.set_export_param_value('display_mode', display_mode)
    itol_exporter.set_export_param_value('range_mode', '2')
    itol_exporter.set_export_param_value('dashed_lines', '0')
    itol_exporter.set_export_param_value('format', op_plot_ext)
    itol_exporter.export(op_plot)


tree_file           = '/Users/songweizhi/Desktop/DateArTree/0_HGT_ALE/ale_op_dir/OG0000006_genome_tree_for_ALE.treefile'
annotation_file_1   = '/Users/songweizhi/Desktop/DateArTree/0_HGT_ALE/ale_op_dir/OG0000006_iTOL_connection.txt'
tree_plot           = '/Users/songweizhi/Desktop/DateArTree/0_HGT_ALE/ale_op_dir/OG0000006_genome_tree_for_ALE.treefile.png'
project_name        = 'batch_access_tmp'
API_key             = 'S1kZZuDHc0d5M7J5vLnUNQ'
display_mode        = '1'  # # 1=rectangular, 2=circular, 3=unrooted

annotation_file_list = [annotation_file_1]
itol_tree(tree_file, annotation_file_list, project_name, API_key, display_mode, tree_plot)
