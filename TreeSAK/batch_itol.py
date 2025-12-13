import os
import glob
import argparse
from pathlib import Path
from itolapi import Itol


batch_itol_usage = '''
======================================= batch_itol example commands =======================================

TreeSAK batch_itol -f -api API_key -ip batch_access_tmp -a annotation_files.txt -i input.tree -o out.pdf
TreeSAK batch_itol -f -api API_key -ip batch_access_tmp -a annotation_files.txt -i tree_dir -x tree -o out_pdf

Manual
https://github.com/albertyw/itolapi
http://itol.embl.de/help.cgi#batch

# An example of the parameter file is available here
# to be added

===========================================================================================================
'''

def sep_path_basename_ext(file_in):

    f_path, f_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'
    f_base, f_ext = os.path.splitext(f_name)
    f_ext = f_ext[1:]

    return f_name, f_path, f_base, f_ext


def itol_single_tree(tree_file, annotation_files_txt, project_name, APIkey, parameter_dict, op_plot):

    annotation_file_list = []
    for each_file in open(annotation_files_txt):
        annotation_file_list.append(each_file.strip())

    op_plot_ext = op_plot.split('.')[-1]

    # upload tree to iTOL
    itol_uploader = Itol()
    itol_uploader.params['projectName'] = project_name  # better to create a project with a unique name.
    itol_uploader.params['APIkey']      = APIkey  # sine we are the same account, we can use the same APIkey
    itol_uploader.params['treeName']    = tree_file
    itol_uploader.add_file(Path(tree_file))

    # upload annotation files to iTOL
    for annotation_file in annotation_file_list:
        itol_uploader.add_file(Path(annotation_file))

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

    parameter_dict.get('', 'to be added')
    parameter_dict.get('', '')


    # for a full list of options, go to https://itol.embl.de/help.cgi#batchExp
    itol_exporter = itol_uploader.get_itol_export()
    itol_exporter.set_export_param_value('internal_scale', parameter_dict.get('internal_scale', '0'))
    itol_exporter.set_export_param_value('datasets_visible', datasets_visible_str)
    itol_exporter.set_export_param_value('display_mode', parameter_dict.get('display_mode', '1'))
    itol_exporter.set_export_param_value('vertical_shift_factor', parameter_dict.get('vertical_shift_factor', '1'))
    itol_exporter.set_export_param_value('horizontal_scale_factor', parameter_dict.get('horizontal_scale_factor', '0.9'))

    # range
    itol_exporter.set_export_param_value('range_mode', parameter_dict.get('range_mode', '2'))                                           # Possible values: 0,1 or 2 (0=off, 1=cover labels only, 2=cover full clades)
    itol_exporter.set_export_param_value('include_ranges_legend', parameter_dict.get('include_ranges_legend', '0'))

    # label
    # itol_exporter.set_export_param_value('current_font_size', '12') # the default looks good
    itol_exporter.set_export_param_value('current_font_name', parameter_dict.get('current_font_name', 'Courier'))
    itol_exporter.set_export_param_value('default_label_color', parameter_dict.get('default_label_color', '#000000'))

    # branch
    itol_exporter.set_export_param_value('line_width', parameter_dict.get('line_width', '2'))
    itol_exporter.set_export_param_value('dashed_lines', parameter_dict.get('dashed_lines', '1'))
    itol_exporter.set_export_param_value('default_branch_color', parameter_dict.get('default_branch_color', '#000000'))

    # bootstrap
    itol_exporter.set_export_param_value('metadata_source', parameter_dict.get('metadata_source', 'bootstrap'))                         # Which metadata source to use for bootstrap display options
    itol_exporter.set_export_param_value('bootstrap_display', parameter_dict.get('bootstrap_display', '1'))                             # possible values: 0 or 1
    itol_exporter.set_export_param_value('bootstrap_type', parameter_dict.get('bootstrap_type', '2'))                                   # Possible values: 1, 2, 3 or 4 (1=Symbol, 2=Text label, 3=Branch color and 4=Branch width)
    itol_exporter.set_export_param_value('bootstrap_label_size', parameter_dict.get('bootstrap_label_size', '15'))                      # in pixels, integer >= 9
    itol_exporter.set_export_param_value('bootstrap_label_percent_factor', parameter_dict.get('bootstrap_label_percent_factor', '10'))  # in pixels, integer >= 9

    # write out
    itol_exporter.set_export_param_value('format', op_plot_ext)
    itol_exporter.export(op_plot)


def batch_itol(args):

    tree_file_dir           = args['i']
    tree_file_ext           = args['x']
    annotation_files_txt    = args['a']
    op_file_dir             = args['o']
    force_overwrite         = args['f']
    API_key                 = args['api']
    project_name            = args['ip']
    para_txt                = args['para']

    para_dict = dict()
    if para_txt is not None:
        if os.path.isfile(para_txt) is False:
            print('The specified parameter file does not exist, program exited!')
            exit()
        else:
            for each_line in open(para_txt):
                if not each_line.startswith('#'):
                    if len(each_line.strip()) > 0:
                        para_without_comment = each_line.strip().split('#')[0].strip()
                        para_without_comment_split = para_without_comment.split('\t')
                        para_dict[para_without_comment_split[0]] = para_without_comment_split[1]

    if os.path.isfile(tree_file_dir) is True:
        itol_single_tree(tree_file_dir, annotation_files_txt, project_name, API_key, para_dict, op_file_dir)
    elif os.path.isdir(tree_file_dir) is True:
        file_re = '%s/*.%s' % (tree_file_dir, tree_file_ext)
        file_list = glob.glob(file_re)

        if len(file_list) == 0:
            print('no file found in %s, please check file extension, program exited!' % tree_file_dir)
            exit()

        # create output folder
        if os.path.isdir(op_file_dir) is True:
            if force_overwrite is True:
                os.system('rm -r %s' % op_file_dir)
            else:
                print('Output folder detected, program exited!')
                exit()
        os.system('mkdir %s' % op_file_dir)

        for each_file in sorted(file_list):
            f_name, f_path, f_base, f_ext = sep_path_basename_ext(each_file)
            op_pdf = '%s/%s.pdf' % (op_file_dir, f_base)

            itol_single_tree(each_file, annotation_files_txt, project_name, API_key, para_dict, op_pdf)
    else:
        print('please provide input file with -i, program exited!')
        exit()


if __name__ == '__main__':

    batch_itol_parser = argparse.ArgumentParser(usage=batch_itol_usage)
    batch_itol_parser.add_argument('-i',        required=True,                              help='input tree file or folder')
    batch_itol_parser.add_argument('-x',        required=False, default=None,               help='file extension')
    batch_itol_parser.add_argument('-o',        required=True,                              help='output file or folder')
    batch_itol_parser.add_argument('-a',        required=False, default=None,               help='a txt file contain absolute to all annotation files')
    batch_itol_parser.add_argument('-para',     required=False, default=None,               help='parameter file')
    batch_itol_parser.add_argument('-api',      required=True,                              help='iTOL API key')
    batch_itol_parser.add_argument('-ip',       required=False, default='batch_access_tmp', help='iTOL project name, default: batch_access_tmp')
    batch_itol_parser.add_argument('-f',        required=False, action="store_true",        help='force overwrite')
    args = vars(batch_itol_parser.parse_args())
    batch_itol(args)
