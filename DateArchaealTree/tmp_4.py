import os
from itolapi import Itol


def sep_path_basename_ext(file_in):

    # separate path and file name
    file_path, file_name = os.path.split(file_in)
    if file_path == '':
        file_path = '.'

    # separate file basename and extension
    file_basename, file_extension = os.path.splitext(file_name)

    return file_path, file_basename, file_extension


def itol_tree(tree_file, annotation_file_list, project_name, APIkey, display_mode, op_plot):

    # https://github.com/albertyw/itolapi
    # http://itol.embl.de/help.cgi#batch

    op_plot_ext = op_plot.split('.')[-1]
    tree_file_path, tree_file_basename, tree_file_extension = sep_path_basename_ext(tree_file)

    # upload tree to iTOL
    itol_uploader = Itol()
    itol_uploader.params['projectName'] = project_name  # better to create a project with a unique name.
    itol_uploader.params['APIkey'] = APIkey  # sine we are the same account, we can use the same APIkey
    itol_uploader.params['treeName'] = tree_file_basename
    itol_uploader.add_file(tree_file)

    # upload annotation files to iTOL
    for annotation_file in annotation_file_list:
        if os.path.isfile(annotation_file) is True:
            itol_uploader.add_file(annotation_file)
        else:
            print('%s not found!' % annotation_file)

    status = itol_uploader.upload()
    # import pdb;pdb.set_trace()
    assert status != False

    # the following parameters are optional, refer to https://itol.embl.de/help.cgi#batchExp
    itol_exporter = itol_uploader.get_itol_export()
    itol_exporter.set_export_param_value('datasets_visible', '0')
    itol_exporter.set_export_param_value('display_mode', display_mode)
    itol_exporter.set_export_param_value('range_mode', '2')
    itol_exporter.set_export_param_value('dashed_lines', '0')
    itol_exporter.set_export_param_value('format', op_plot_ext)
    itol_exporter.export(op_plot)


wd                              = '/Users/songweizhi/Desktop/DateArTree/0_HGT_ALE/ale_op_dir'
og_id                           = 'OG0000007'

tree_file                       = '%s/%s_for_ALE.ufboot_ALE_renamed_genome_tree.tree'           % (wd, og_id)
plot_file                       = '%s/%s_for_ALE.ufboot_ALE_renamed_genome_tree.tree.png'       % (wd, og_id)
connection_file                 = '%s/%s_iTOL_connection.txt'                                   % (wd, og_id)
project_name                    = 'batch_access_tmp'
API_key                         = 'S1kZZuDHc0d5M7J5vLnUNQ'
display_mode                    = '1'  # # 1=rectangular, 2=circular, 3=unrooted


tree_file                       = '/Users/songweizhi/Desktop/DateArTree/0_HGT_ALE/ale_op_dir/OG0000007_with_both_renamed.tree'
plot_file                       = '/Users/songweizhi/Desktop/DateArTree/0_HGT_ALE/ale_op_dir/OG0000007_with_both_renamed.tree.png'
connection_file                 = '/Users/songweizhi/Desktop/DateArTree/0_HGT_ALE/ale_op_dir/OG0000007_iTOL_connection_renamed.txt'
leaf_label_file                 = '/Users/songweizhi/Desktop/DateArTree/0_HGT_ALE/ale_op_dir/OG0000007_iTOL_genome_pco.txt'

annotation_file_list = [connection_file, leaf_label_file]
itol_tree(tree_file, annotation_file_list, project_name, API_key, display_mode, plot_file)

