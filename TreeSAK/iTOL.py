import math
import os.path
import random
import argparse
import seaborn as sns


'''
TreeSAK iTOL -ColorLabel -lg Mag_genus.txt -gc genus_color.txt-o ColorLabel_genus.txt
'''

iTOL_usage = '''
==================================== iTOL example commands ====================================

# Example commands
TreeSAK iTOL -Labels -ll new_Mag_name.txt -o Mag_name_iTOL_Labels.txt
TreeSAK iTOL -ColoredLabel -lg leaf_group.txt -gc group_color.txt -o iTOL_ColoredLabel.txt
TreeSAK iTOL -ColoredLabel -lg leaf_group.txt -gc group_color.txt -o iTOL_ColoredLabel.txt -ll new_Mag_name.txt
TreeSAK iTOL -Binary -lm Binary_matrix.txt -lt Enzyme -o Presence_Absence_iTOL.txt -cc "#F1BB83"
TreeSAK iTOL -Binary -lm Binary_matrix.txt -lt Enzyme -o Presence_Absence_iTOL.txt -cc col_color.txt
TreeSAK iTOL -BinaryID -id MagID.txt -lt dRep95 -o dRep95_representatives_iTOL.txt
TreeSAK iTOL -Heatmap -lm MagAbundance.txt -lt Abundance -o Heatmap_abundance.txt
TreeSAK iTOL -SimpleBar -lv MagSize.txt -scale 0-3-6-9 -lt Size -o SimpleBar_size.txt
TreeSAK iTOL -ColorStrip -lg MagTaxon.txt -lt Phylum -gc phylum_color.txt -o ColorStrip_taxon.txt
TreeSAK iTOL -ColorRange -lg MagTaxon.txt -lt Phylum -gc phylum_color.txt -o ColorRange_taxon.txt
TreeSAK iTOL -ColorRange -taxon Taxonomy.txt -rank f -lt Family -o ColorRange_taxon.txt
TreeSAK iTOL -ColorClade -lg Mag_genus.txt -gc genus_color.txt-o ColorClade_genus.txt
TreeSAK iTOL -ExternalShape -lm identity_matrix.txt -lt Identity -scale 25-50-75-100 -o ExternalShape_identity.txt
TreeSAK iTOL -PieChart -lv MagCompleteness.txt -lt Completeness -o PieChart_completeness.txt
TreeSAK iTOL -Collapse -lg MagTaxon.txt -o Collapse_by_taxon.txt

# Leaf-to-Group file format (-lg, tab separated, no header)
genome_1	Bacteria
genome_2	Archaea

# taxonomy file format (-taxon, tab separated, GTDB-format taxononomy string)
genome_1	d__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Dongiales;f__Dongiaceae;g__Dongia;s__Dongia mobilis
genome_2	d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Arenicellales;f__LS-SOB;g__VYGS01;s__

# Group-to-Color and Column-to-Color file format (-gc and -cc, tab separated, no header)
bac	#CCCC00
ar	#9999FF
# Please note you can only specify one color for Binary data, provide with -gc lightblue or -gc "#85C1E9"

# Leaf-to-Value file format (-lv, tab separated, no header)
genome_1	6.15
genome_2	6.63

# Leaf-to-Matrix file format (-lm, tab separated, header required!!!)
Genome_id SampleA   SampleB   SampleC
genome_1	6.15    2.23    1.56
genome_2	6.63    1.72    2.55

===============================================================================================
'''


def get_color_list(color_num):

    if color_num <= 8:
        color_list_combined = ['#3787c0', '#39399f', '#ffb939', '#399f39', '#9f399f', '#fb694a', '#9f9f39', '#959595']

    elif 8 < color_num <= 16:
        color_list_combined = ['#2b7bba', '#89bedc', '#2e2e99', '#8a8acc', '#ffa500', '#ffc55c', '#2e992e', '#8acc8a', '#992e99', '#cc8acc', '#d52221', '#fc8161', '#99992e', '#cccc8a', '#5c5c5c', '#adadad']

    else:
        color_num_each = math.ceil(color_num/8) + 2
        color_list_1 = sns.color_palette('Blues', n_colors=color_num_each).as_hex()
        color_list_2 = sns.light_palette('navy',   n_colors=color_num_each).as_hex()
        color_list_3 = sns.light_palette('orange', n_colors=color_num_each).as_hex()
        color_list_4 = sns.light_palette('green',  n_colors=color_num_each).as_hex()
        color_list_5 = sns.light_palette('purple', n_colors=color_num_each).as_hex()
        color_list_6 = sns.color_palette('Reds',  n_colors=color_num_each).as_hex()
        color_list_7 = sns.light_palette('olive',  n_colors=color_num_each).as_hex()
        color_list_8 = sns.color_palette('Greys', n_colors=color_num_each).as_hex()

        color_list_combined = []
        for color_list in [color_list_1, color_list_2, color_list_3, color_list_4, color_list_5, color_list_6, color_list_7, color_list_8]:
            for color in color_list[2:][::-1]:
                color_list_combined.append(color)

    color_list_to_return = random.sample(color_list_combined, color_num)

    color_list_to_return_sorted = []
    for color_to_return in color_list_combined:
        if color_to_return in color_list_to_return:
            color_list_to_return_sorted.append(color_to_return)

    random.shuffle(color_list_to_return_sorted)

    return color_list_to_return_sorted


def scale_str_to_size_list(scale_str):

    scale_list = scale_str.split('-')
    scale_list = [float(i) for i in scale_list]

    shape_size_list = []
    if scale_list[0] == 0:
        shape_size_list = [0]
        for each_value in scale_list[1:-1]:
            current_size = each_value/scale_list[-1]
            shape_size_list.append(current_size)
        shape_size_list.append(1)

    if scale_list[0] != 0:
        shape_size_list = [0.1]
        interval_num = len(scale_list) - 1
        interval_value = (1 - 0.1)/interval_num
        n = 1
        for each_value in scale_list[1:-1]:
            shape_size_list.append(interval_value * n + 0.1)
            n += 1
        shape_size_list.append(1)

    return shape_size_list


def iTOL(args):

    # read in arguments
    Labels                  = args['Labels']
    ColoredLabel            = args['ColoredLabel']
    ColorStrip              = args['ColorStrip']
    ColorRange              = args['ColorRange']
    ColorClade              = args['ColorClade']
    ColorLabel              = args['ColorLabel']
    SimpleBar               = args['SimpleBar']
    Heatmap                 = args['Heatmap']
    ExternalShape           = args['ExternalShape']
    Binary                  = args['Binary']
    BinaryID                = args['BinaryID']
    Connection              = args['Connection']
    PieChart                = args['PieChart']
    Collapse                = args['Collapse']
    MultiStyleLabel         = args['MultiStyleLabel']
    leaf_id_txt             = args['id']
    LeafGroup               = args['lg']
    GroupColor              = args['gc']
    ColumnColor_txt         = args['cc']
    LeafValue               = args['lv']
    LeafLabel               = args['ll']
    LeafMatrix              = args['lm']
    d2r                     = args['dr']
    scale_str               = args['scale']
    show_color_strip_legend = args['legend']
    LegendTitle             = args['lt']
    show_strip_labels       = args['show_strip_labels']
    FileOut                 = args['o']
    BinaryShape             = args['BinaryShape']
    BinaryColor             = args['BinaryColor']

    # General
    STRIP_WIDTH                 = 100
    MARGIN                      = 20
    branch_width                = 2

    # SimpleBar
    SimpleBar_COLOR             = 'grey'
    SimpleBar_WIDTH             = 300
    SimpleBar_HEIGHT_FACTOR     = 0.8
    SimpleBar_BORDER_WIDTH      = 0
    SimpleBar_SCALE_COLOR       = '#696969'
    SimpleBar_SCALE_WIDTH       = 1
    SimpleBar_SCALE_DASHED      = 1
    SimpleBar_SCALE_FontSize    = 2

    # Heatmap
    Heatmap_STRIP_WIDTH         = 30

    # check the number of specified file type
    True_num = 0
    for file_type in [Labels, ColorStrip, ColorRange, SimpleBar, Heatmap, ExternalShape, Binary, BinaryID, PieChart, Collapse, ColorClade, ColorLabel, ColoredLabel]:
        if file_type is True:
            True_num += 1

    if True_num == 0:
        print('Please specify one file type, choose from -ColorStrip, -ColorRange, -SimpleBar, -Heatmap, -ExternalShape, -Binary, -BinaryID or -Collapse')
        exit()
    if True_num > 1:
        print('Please specify one file type ONLY, choose from -ColorStrip, -ColorRange, -SimpleBar, -Heatmap, -ExternalShape, -Binary, -BinaryID or -Collapse')
        exit()

    ####################################################################################################################

    if ColoredLabel is True:

        leaf_set = set()
        leaf_label_dict = dict()
        if LeafLabel is not None:
            for each in open(LeafLabel):
                each_split = each.strip().split('\t')
                leaf_label_dict[each_split[0]] = each_split[1]
                leaf_set.add(each_split[0])

        leaf_group_dict = dict()
        for each in open(LeafGroup):
            each_split = each.strip().split('\t')
            leaf_group_dict[each_split[0]] = each_split[1]
            leaf_set.add(each_split[0])

        group_color_dict = dict()
        for each in open(GroupColor):
            each_split = each.strip().split('\t')
            group_color_dict[each_split[0]] = each_split[1]

        FileOut_handle = open(FileOut, 'w')
        FileOut_handle.write('DATASET_TEXT\n')
        FileOut_handle.write('SEPARATOR TAB\n')
        FileOut_handle.write('DATASET_LABEL\t%s\n' % 'ColoredLabel')
        FileOut_handle.write('\nDATA\n')
        for leaf in leaf_set:
            leaf_label = leaf_label_dict.get(leaf, leaf)
            leaf_group = leaf_group_dict.get(leaf, 'na')
            leaf_color = group_color_dict.get(leaf_group, '#000000')
            FileOut_handle.write('%s\t%s\t-1\t%s\tnormal\t1\t0\n' % (leaf, leaf_label, leaf_color))
        FileOut_handle.close()

    ####################################################################################################################

    # Prepare ColorStrip and ColorRange file
    if (ColorStrip is True) or (ColorRange is True):

        Leaf_to_Group_dict = {}
        Group_list = []
        for each_leaf in open(LeafGroup):
            each_leaf_split = each_leaf.strip().split('\t')
            Leaf_to_Group_dict[each_leaf_split[0]] = each_leaf_split[1]
            if each_leaf_split[1] not in Group_list:
                Group_list.append(each_leaf_split[1])

        Group_to_Color_dict = {}
        if GroupColor is None:
            color_list = get_color_list(len(Group_list))
            Group_to_Color_dict = dict(zip(Group_list, color_list))
        else:
            # get groups with provided color
            Group_to_provided_Color_dict = dict()
            group_with_provided_color_list = []
            for each_group in open(GroupColor):
                each_group_split = each_group.strip().split('\t')
                group_id = each_group_split[0]
                color_code = each_group_split[1]
                if group_id in Group_list:
                    Group_to_provided_Color_dict[group_id] = color_code
                    group_with_provided_color_list.append(group_id)

            # assign colors to the rest groups
            group_without_color_list = []
            for each_group in Group_list:
                if each_group not in group_with_provided_color_list:
                    group_without_color_list.append(each_group)
            if len(group_without_color_list) > 0:
                color_list_unprovided = get_color_list(len(group_without_color_list))
                Group_to_Color_dict_unprovided = dict(zip(group_without_color_list, color_list_unprovided))
                for each_group in Group_to_Color_dict_unprovided:
                    Group_to_Color_dict[each_group] = Group_to_Color_dict_unprovided[each_group]

            # combine two dict
            for each_group in Group_to_provided_Color_dict:
                Group_to_Color_dict[each_group] = Group_to_provided_Color_dict[each_group]

        group_list = [i for i in Group_to_Color_dict]
        color_list = [Group_to_Color_dict[i] for i in group_list]

        FileOut_handle = open(FileOut, 'w')

        # write out header
        if ColorStrip is True:
            FileOut_handle.write('DATASET_COLORSTRIP\n')
            FileOut_handle.write('SEPARATOR TAB\n')
            FileOut_handle.write('DATASET_LABEL\t%s\n' % LegendTitle)
            FileOut_handle.write('SHOW_LABELS\t1\n')
            FileOut_handle.write('LABEL_SHIFT\t30\n')
            FileOut_handle.write('LABEL_ROTATION\t45\n')
            FileOut_handle.write('SIZE_FACTOR\t2\n')
            FileOut_handle.write('COLOR_BRANCHES\t0\n')

        if ColorRange is True:
            FileOut_handle.write('TREE_COLORS\n')
            FileOut_handle.write('SEPARATOR TAB\n')
            FileOut_handle.write('DATASET_LABEL\t%s_ColorRange\n' % LegendTitle)

        # write out strip attributes
        if ColorStrip is True:
            FileOut_handle.write('\n# customize strip attributes here\n')
            FileOut_handle.write('STRIP_WIDTH\t%s\n' % STRIP_WIDTH)
            FileOut_handle.write('MARGIN\t%s\n'      % MARGIN)

            if show_color_strip_legend is False:
                FileOut_handle.write('\n# customize labels on the trip\n')
                if show_strip_labels is True:
                    FileOut_handle.write('SHOW_STRIP_LABELS\t1\n')
                else:
                    FileOut_handle.write('SHOW_STRIP_LABELS\t0\n')
                FileOut_handle.write('STRIP_LABEL_POSITION\tcenter\n')
                FileOut_handle.write('STRIP_LABEL_ROTATION\t90\n')
            else:
                FileOut_handle.write('\n# uncomment to show labels on the trip\n')
                if show_strip_labels is True:
                    FileOut_handle.write('# SHOW_STRIP_LABELS\t1\n')
                else:
                    FileOut_handle.write('# SHOW_STRIP_LABELS\t0\n')
                FileOut_handle.write('# STRIP_LABEL_POSITION\tcenter\n')
                FileOut_handle.write('# STRIP_LABEL_ROTATION\t90\n')

        # write out legend info
        if show_color_strip_legend is True:
            FileOut_handle.write('\n# customize legend\n')
            FileOut_handle.write('LEGEND_TITLE\t%s\n' % LegendTitle)
            FileOut_handle.write('LEGEND_SHAPES\t%s\n' % '\t'.join(['1' for i in Group_to_Color_dict]))
            FileOut_handle.write('LEGEND_COLORS\t%s\n' % '\t'.join(color_list))
            FileOut_handle.write('LEGEND_LABELS\t%s\n' % '\t'.join(group_list))
        else:
            FileOut_handle.write('\n# uncomment if you want to show the legend\n')
            FileOut_handle.write('# LEGEND_TITLE\t%s\n' % LegendTitle)
            FileOut_handle.write('# LEGEND_SHAPES\t%s\n' % '\t'.join(['1' for i in Group_to_Color_dict]))
            FileOut_handle.write('# LEGEND_COLORS\t%s\n' % '\t'.join(color_list))
            FileOut_handle.write('# LEGEND_LABELS\t%s\n' % '\t'.join(group_list))

        FileOut_handle.write('\n# provide data here\nDATA\n')
        for leaf in Leaf_to_Group_dict:
            node_group = Leaf_to_Group_dict[leaf]
            leaf_color = Group_to_Color_dict[node_group]

            if ColorStrip is True:
                FileOut_handle.write('%s\t%s\t%s\n' % (leaf, leaf_color, node_group))
            if ColorRange is True:
                FileOut_handle.write('%s\trange\t%s\t%s\n' % (leaf, leaf_color, node_group))
        FileOut_handle.close()

    ####################################################################################################################

    # Prepare SimpleBar file
    if SimpleBar is True:

        if scale_str is None:
            print('Please provide scale values for barchart, in format 0-3-6-9')
            exit()

        # read in leaf value into dict
        leaf_value_dict = {}
        max_value = None
        for leaf_value in open(LeafValue):
            leaf_value_split = leaf_value.strip().split('\t')
            leaf_value_dict[leaf_value_split[0]] = float(leaf_value_split[1])

            # get max value
            if max_value == None:
                max_value = float(leaf_value_split[1])
            else:
                if float(leaf_value_split[1]) > max_value:
                    max_value = float(leaf_value_split[1])

        SimpleBar_FileOut_handle = open(FileOut, 'w')
        SimpleBar_FileOut_handle.write('DATASET_SIMPLEBAR\n')
        SimpleBar_FileOut_handle.write('# Reference: https://itol.embl.de/help/dataset_simplebar_template.txt\n')
        SimpleBar_FileOut_handle.write('\nSEPARATOR TAB\n')
        SimpleBar_FileOut_handle.write('\n# customize barchart attributes here\n')
        SimpleBar_FileOut_handle.write('DATASET_LABEL\t%s\n'    % LegendTitle)
        SimpleBar_FileOut_handle.write('COLOR\t%s\n'            % SimpleBar_COLOR)
        SimpleBar_FileOut_handle.write('WIDTH\t%s\n'            % SimpleBar_WIDTH)
        SimpleBar_FileOut_handle.write('MARGIN\t%s\n'           % MARGIN)
        SimpleBar_FileOut_handle.write('HEIGHT_FACTOR\t%s\n'    % SimpleBar_HEIGHT_FACTOR)
        SimpleBar_FileOut_handle.write('BORDER_WIDTH\t%s\n'     % SimpleBar_BORDER_WIDTH)

        # write out scale attributes
        scale_attributes_list = []
        for scale_value in scale_str.split('-'):
            scale_attributes = '%s-%s-%s-%s-%s-%s' % (scale_value, scale_value, SimpleBar_SCALE_COLOR, SimpleBar_SCALE_WIDTH, SimpleBar_SCALE_DASHED, SimpleBar_SCALE_FontSize)
            scale_attributes_list.append(scale_attributes)

        SimpleBar_FileOut_handle.write('\n# customize scale attributes here\n')
        SimpleBar_FileOut_handle.write('# format: VALUE-LABEL-COLOR-WIDTH-DASHED-LABEL_SCALE_FACTOR, LABEL_SCALE_FACTOR controls font size\n')
        SimpleBar_FileOut_handle.write('DATASET_SCALE\t%s\n' % '\t'.join(scale_attributes_list))
        SimpleBar_FileOut_handle.write('\n# provide data here\n')
        SimpleBar_FileOut_handle.write('DATA\n')

        for leaf in leaf_value_dict:
            SimpleBar_FileOut_handle.write('%s\t%s\n' % (leaf, leaf_value_dict[leaf]))

        SimpleBar_FileOut_handle.close()

    ####################################################################################################################

    # Prepare Binary file
    if Binary is True:

        Column_to_Color_dict = {}
        ColumnColor = 'red'
        if ColumnColor_txt is not None:
            if os.path.isfile(ColumnColor_txt) is True:
                for each_col in open(ColumnColor_txt):
                    each_col_split = each_col.strip().split('\t')
                    Column_to_Color_dict[each_col_split[0]] = each_col_split[1]
            else:
                ColumnColor = ColumnColor_txt

        # get color
        if GroupColor is not None:
            if os.path.isfile(GroupColor) is True:
                print('Only one color can be specified for Binary data, here are two examples:\n-gc lightblue\n-gc "#85C1E9"')
                print('Program exited!')
                exit()
            else:
                Binary_color = GroupColor
        else:
            Binary_color = 'red'

        Binary_FileOut_handle = open(FileOut, 'w')
        line_index = 0
        for each_line in open(LeafMatrix):
            each_line_split = each_line.strip().split('\t')
            if line_index == 0:

                # get header list
                if each_line.startswith('\t'):
                    col_name_list = each_line_split
                else:
                    col_name_list = each_line_split[1:]

                col_color_list = [Column_to_Color_dict.get(i, ColumnColor) for i in col_name_list]

                Binary_FileOut_handle.write('DATASET_BINARY\n\nSEPARATOR TAB\nDATASET_LABEL\t%s\nCOLOR\t%s\n' % (LegendTitle, Binary_color))
                Binary_FileOut_handle.write('SHOW_LABELS\t1\nLABEL_ROTATION\t45\nLABEL_SHIFT\t5\n')
                Binary_FileOut_handle.write('FIELD_LABELS\t%s\n' % '\t'.join(col_name_list))
                Binary_FileOut_handle.write('FIELD_COLORS\t%s\n' % '\t'.join(col_color_list))
                Binary_FileOut_handle.write('MARGIN\t10\n')
                Binary_FileOut_handle.write('HORIZONTAL_GRID\t0\n')
                Binary_FileOut_handle.write('VERTICAL_GRID\t0\n')
                Binary_FileOut_handle.write('FIELD_SHAPES\t%s\n' % '\t'.join(['2'] * len(col_name_list)))
                Binary_FileOut_handle.write('#1: rectangle\n')
                Binary_FileOut_handle.write('#2: circle\n')
                Binary_FileOut_handle.write('#3: star\n')
                Binary_FileOut_handle.write('#4: right pointing triangle\n')
                Binary_FileOut_handle.write('#5: left pointing triangle\n')
                Binary_FileOut_handle.write('#6: check mark\n')
                Binary_FileOut_handle.write('\nDATA\n')
            else:
                Binary_FileOut_handle.write(each_line)
            line_index += 1
        Binary_FileOut_handle.close()

    ####################################################################################################################

    # Prepare BinaryID file
    if BinaryID is True:
        BinaryID_FileOut_handle = open(FileOut, 'w')
        BinaryID_FileOut_handle.write('DATASET_BINARY\n\nSEPARATOR TAB\nDATASET_LABEL\t%s\nCOLOR\tred\n' % LegendTitle)
        BinaryID_FileOut_handle.write('SHOW_LABELS\t1\nLABEL_ROTATION\t45\nLABEL_SHIFT\t5\n')
        BinaryID_FileOut_handle.write('FIELD_LABELS\t%s\n' % LegendTitle)
        BinaryID_FileOut_handle.write('FIELD_COLORS\t%s\n' % BinaryColor)
        BinaryID_FileOut_handle.write('# FIELD_SHAPES: 1: rectangle; 2: circle; 3: star; 4: right pointing triangle; 5: left pointing triangle; 6: check mark\n')
        BinaryID_FileOut_handle.write('FIELD_SHAPES\t%s\n' % BinaryShape)
        BinaryID_FileOut_handle.write('MARGIN\t10\n')
        BinaryID_FileOut_handle.write('HORIZONTAL_GRID\t0\n')
        BinaryID_FileOut_handle.write('VERTICAL_GRID\t0\n')
        BinaryID_FileOut_handle.write('\nDATA\n')
        for each_line in open(leaf_id_txt):
            leaf_id = each_line.strip().split()[0]
            BinaryID_FileOut_handle.write('%s\t1\n' % leaf_id)
        BinaryID_FileOut_handle.close()

    ####################################################################################################################

    # Prepare Heatmap file
    if Heatmap is True:

        n = 0
        col_name_list = []
        leaf_matrix_dict = {}
        for leaf_matrix in open(LeafMatrix):
            leaf_matrix_split = leaf_matrix.strip().split('\t')
            if n == 0:
                col_name_list = leaf_matrix_split[1:]
            else:
                leaf_matrix_dict[leaf_matrix_split[0]] = leaf_matrix_split[1:]
            n += 1

        Heatmap_FileOut_handle = open(FileOut, 'w')
        Heatmap_FileOut_handle.write('DATASET_HEATMAP\n')
        Heatmap_FileOut_handle.write('# Reference https://itol.embl.de/help/dataset_heatmap_template.txt\n')
        Heatmap_FileOut_handle.write('\nSEPARATOR TAB\n')
        Heatmap_FileOut_handle.write('\n# customize heatmap attributes here\n')
        Heatmap_FileOut_handle.write('MARGIN\t%s\n'          % MARGIN)
        Heatmap_FileOut_handle.write('STRIP_WIDTH\t%s\n'     % Heatmap_STRIP_WIDTH)
        Heatmap_FileOut_handle.write('\n# customize legend here\n')
        Heatmap_FileOut_handle.write('AUTO_LEGEND\t1\n')
        Heatmap_FileOut_handle.write('DATASET_LABEL\t%s\n' % LegendTitle)
        Heatmap_FileOut_handle.write('USE_MID_COLOR\t1\n')
        Heatmap_FileOut_handle.write('COLOR_MIN\t#2980B9\n')
        Heatmap_FileOut_handle.write('COLOR_MID\t#ECF0F1\n')
        Heatmap_FileOut_handle.write('COLOR_MAX\t#E74C3C\n')
        Heatmap_FileOut_handle.write('\n# customize value range here. By default, color gradients will be calculated based on dataset values\n')
        Heatmap_FileOut_handle.write('# USER_MIN_VALUE	0\n')
        Heatmap_FileOut_handle.write('# USER_MID_VALUE	5\n')
        Heatmap_FileOut_handle.write('# USER_MAX_VALUE	10\n')
        Heatmap_FileOut_handle.write('\n# customize column name here\n')
        Heatmap_FileOut_handle.write('FIELD_LABELS\t%s\n' % '\t'.join(col_name_list))
        Heatmap_FileOut_handle.write('\n# Provide data here\n')
        Heatmap_FileOut_handle.write('DATA\n')
        for leaf in leaf_matrix_dict:
            Heatmap_FileOut_handle.write('%s\t%s\n' % (leaf, '\t'.join(leaf_matrix_dict[leaf])))
        Heatmap_FileOut_handle.close()

    ####################################################################################################################

    if Labels is True:
        if os.path.isfile(LeafLabel) is False:
            print('leaf to label file not found, please provide with -ll, program exited!')
            exit()

        # write out header
        Labels_FileOut_handle = open(FileOut, 'w')
        Labels_FileOut_handle.write('LABELS\nSEPARATOR TAB\nDATA\n')
        for each_ll in open(LeafLabel):
            Labels_FileOut_handle.write(each_ll)
        Labels_FileOut_handle.close()

    ####################################################################################################################

    if Connection is True:

        if os.path.isfile(d2r) is False:
            print('donor to recipient file not found, please provide with -dr, program exited!')
            exit()

        Connection_FileOut_handle = open(FileOut, 'w')
        Connection_FileOut_handle.write('DATASET_CONNECTION\nSEPARATOR TAB\nDATASET_LABEL\tdemo_connections\n')
        Connection_FileOut_handle.write('COLOR\t#ff0ff0\nDRAW_ARROWS\t1\nARROW_SIZE\t60\nLOOP_SIZE\t100\n')
        Connection_FileOut_handle.write('MAXIMUM_LINE_WIDTH\t10\nCURVE_ANGLE\t45\nCENTER_CURVES\t1\nALIGN_TO_LABELS\t0\nDATA\n')
        for each_connection in open(d2r):
            Connection_FileOut_handle.write(each_connection)
        Connection_FileOut_handle.close()

    ####################################################################################################################

    # Prepare PieChart file
    if PieChart is True:
        PieChart_FileOut_handle = open(FileOut, 'w')
        PieChart_FileOut_handle.write('DATASET_PIECHART\n')
        PieChart_FileOut_handle.write('SEPARATOR TAB\n')
        PieChart_FileOut_handle.write('DATASET_LABEL\t%s\n' % LegendTitle)
        PieChart_FileOut_handle.write('COLOR\t#ff0000\n')
        PieChart_FileOut_handle.write('FIELD_COLORS\t#5DADE2\t#FFFFFF\n')
        PieChart_FileOut_handle.write('FIELD_LABELS\tf1\tf2\n')
        PieChart_FileOut_handle.write('MARGIN\t10\n\n')
        PieChart_FileOut_handle.write('DATA\n')
        for each in open(LeafValue):
            each_split = each.strip().split('\t')
            node_id = each_split[0]
            node_value = float(each_split[1])
            PieChart_FileOut_handle.write('%s\t-1\t1\t%s\t%s\n' % (node_id, node_value, (100 - node_value)))
        PieChart_FileOut_handle.close()

    ####################################################################################################################

    if (ColorClade is True) or (ColorLabel is True):

        # get group_to_leaf_dict
        group_set = set()
        group_to_leaf_dict = dict()
        node_to_group_dict = dict()
        for each_line in open(LeafGroup):
            each_line_split = each_line.strip().split('\t')
            node_id    = each_line_split[0]
            node_group = each_line_split[1]
            group_set.add(node_group)
            node_to_group_dict[node_id] = node_group
            if node_group not in group_to_leaf_dict:
                group_to_leaf_dict[node_group] = [node_id]
            else:
                group_to_leaf_dict[node_group].append(node_id)

        group_list = sorted(list(group_set))

        Group_to_Color_dict = dict()

        # get groups with provided color
        Group_to_provided_Color_dict = dict()
        group_with_provided_color_list = []
        for each_group in open(GroupColor):
            each_group_split = each_group.strip().split('\t')
            group_id = each_group_split[0]
            color_code = each_group_split[1]
            if group_id in group_list:
                Group_to_provided_Color_dict[group_id] = color_code
                group_with_provided_color_list.append(group_id)

        # assign colors to the rest groups
        group_without_color_list = []
        for each_group in group_list:
            if each_group not in group_with_provided_color_list:
                group_without_color_list.append(each_group)
        if len(group_without_color_list) > 0:
            color_list_unprovided = get_color_list(len(group_without_color_list))
            Group_to_Color_dict_unprovided = dict(zip(group_without_color_list, color_list_unprovided))
            for each_group in Group_to_Color_dict_unprovided:
                Group_to_Color_dict[each_group] = Group_to_Color_dict_unprovided[each_group]

        # combine two dict
        for each_group in Group_to_provided_Color_dict:
            Group_to_Color_dict[each_group] = Group_to_provided_Color_dict[each_group]

        ColorClade_Label_FileOut_handle = open(FileOut, 'w')
        ColorClade_Label_FileOut_handle.write('TREE_COLORS\n')
        ColorClade_Label_FileOut_handle.write('SEPARATOR TAB\n')
        ColorClade_Label_FileOut_handle.write('\nDATA\n')
        for grp in group_to_leaf_dict:
            group_member = group_to_leaf_dict[grp]
            group_color = Group_to_Color_dict[grp]
            concate_str = '|'.join(group_member)
            if ColorClade is True:
                ColorClade_Label_FileOut_handle.write('%s\tclade\t%s\tnormal\t%s\n' % (concate_str, group_color, branch_width))
            if ColorLabel is True:
                ColorClade_Label_FileOut_handle.write('%s\tlabel\t%s\tnormal\t%s\n' % (concate_str, group_color, branch_width))
        ColorClade_Label_FileOut_handle.close()

    ####################################################################################################################

    # Prepare Collapse file
    if Collapse is True:

        corresponding_label_file = '%s.label.txt' % FileOut

        # read in grouping file
        group_to_leaf_dict = dict()
        for each_line in open(LeafGroup):
            each_line_split = each_line.strip().split('\t')
            leaf_id         = each_line_split[0]
            node_group      = each_line_split[1]
            if node_group not in group_to_leaf_dict:
                group_to_leaf_dict[node_group] = [leaf_id]
            else:
                group_to_leaf_dict[node_group].append(leaf_id)

        Collapse_FileOut_handle = open(FileOut, 'w')
        Collapse_FileOut_handle.write('COLLAPSE\n')
        Collapse_FileOut_handle.write('DATA\n')
        label_FileOut_handle = open(corresponding_label_file, 'w')
        label_FileOut_handle.write('LABELS\n')
        label_FileOut_handle.write('SEPARATOR TAB\n')
        label_FileOut_handle.write('DATA\n')
        for each_group in group_to_leaf_dict:
            group_member = group_to_leaf_dict[each_group]
            concate_str = '|'.join(group_member)
            Collapse_FileOut_handle.write(concate_str + '\n')
            label_FileOut_handle.write('%s\t%s\n' % (concate_str, each_group))
        Collapse_FileOut_handle.close()
        label_FileOut_handle.close()

        print('iTOL files exported to:\n%s\n%s' % (FileOut, corresponding_label_file))

    ####################################################################################################################

    # Prepare ExternalShape file
    if  MultiStyleLabel is True:
        MultiStyleLabel_FileOut_handle = open(FileOut, 'w')
        MultiStyleLabel_FileOut_handle.write('MULTI_STYLE\n')
        MultiStyleLabel_FileOut_handle.write('SEPARATOR TAB\n')
        MultiStyleLabel_FileOut_handle.write('#each line in the DATA section defines the style for one label part (specified in the first field)\n')
        MultiStyleLabel_FileOut_handle.write('\n')
        MultiStyleLabel_FileOut_handle.write('\n')
        MultiStyleLabel_FileOut_handle.write('DATA\n')
        MultiStyleLabel_FileOut_handle.close()

    ####################################################################################################################

    # Prepare ExternalShape file
    if ExternalShape is True:

        # read in leaf matrix into dict
        n = 0
        col_name_list = []
        leaf_matrix_dict = {}
        for leaf_matrix in open(LeafMatrix):
            leaf_matrix_split = leaf_matrix.strip().split('\t')
            if n == 0:
                col_name_list = leaf_matrix_split[1:]
            else:
                leaf_matrix_dict[leaf_matrix_split[0]] = leaf_matrix_split[1:]
            n += 1

        Column_to_Color_dict = {}
        if ColumnColor_txt is not None:
            for each_col in open(ColumnColor_txt):
                each_col_split = each_col.strip().split('\t')
                Column_to_Color_dict[each_col_split[0]] = each_col_split[1]

        # check if all columns in data matrix are in Column_to_Color_dict
        if ColumnColor_txt is not None:
            unfound_cols = []
            for each_col_header in col_name_list:
                if each_col_header not in Column_to_Color_dict:
                    unfound_cols.append(each_col_header)
            if len(unfound_cols) > 0:
                print('Color code for the following columns are not provided, program exited!')
                print(','.join(unfound_cols))
                exit()

        ExternalShape_FileOut_handle = open(FileOut, 'w')
        ExternalShape_FileOut_handle.write('DATASET_EXTERNALSHAPE\n')
        ExternalShape_FileOut_handle.write('# Reference https://itol.embl.de/help/dataset_external_shapes_template.txt\n')
        ExternalShape_FileOut_handle.write('\nSEPARATOR TAB\n')

        # define scale here
        if scale_str is None:
            print('Please provide scale values for ExternalShapes, e.g., 25-50-75-100, 2-4-6-8-10')
            exit()

        scale_list = scale_str.split('-')
        LEGEND_SHAPES_list = ['2'] * len(scale_list)
        LEGEND_COLORS_list = ['grey'] * len(scale_list)
        SHAPE_SCALES_list  = scale_str_to_size_list(scale_str)

        ExternalShape_FileOut_handle.write('LEGEND_TITLE\t%s\n' % LegendTitle)
        ExternalShape_FileOut_handle.write('LEGEND_SHAPES\t%s\n' % '\t'.join(LEGEND_SHAPES_list))
        ExternalShape_FileOut_handle.write('LEGEND_COLORS\t%s\n' % '\t'.join(LEGEND_COLORS_list))
        ExternalShape_FileOut_handle.write('LEGEND_LABELS\t%s\n' % '\t'.join(scale_list))
        ExternalShape_FileOut_handle.write('LEGEND_SHAPE_SCALES\t%s\n' % '\t'.join(str(i) for i in SHAPE_SCALES_list))
        ExternalShape_FileOut_handle.write('\n# customize attributes here\n')
        ExternalShape_FileOut_handle.write('VERTICAL_GRID\t0\n')
        ExternalShape_FileOut_handle.write('HORIZONTAL_GRID\t0\n')
        ExternalShape_FileOut_handle.write('SHAPE_TYPE\t2\n')
        ExternalShape_FileOut_handle.write('COLOR_FILL\t1\n')
        ExternalShape_FileOut_handle.write('SHAPE_SPACING\t1\n')
        ExternalShape_FileOut_handle.write('SHOW_LABELS\t0\n')
        ExternalShape_FileOut_handle.write('SHOW_VALUES\t0\n')
        ExternalShape_FileOut_handle.write('DASHED_LINES\t0\n')

        if LegendTitle is None:
            ExternalShape_FileOut_handle.write('DATASET_LABEL\tExternalShape\n')
        else:
            ExternalShape_FileOut_handle.write('DATASET_LABEL\t%s\n' % LegendTitle)

        # write column header/color information
        ExternalShape_FileOut_handle.write('\n# customize column header/color here\n')
        ExternalShape_FileOut_handle.write('FIELD_LABELS\t%s\n' % '\t'.join(col_name_list))

        color_list = get_color_list(len(col_name_list))
        random.shuffle(color_list)
        if ColumnColor_txt is not None:
            color_list = [Column_to_Color_dict[i] for i in col_name_list]
        ExternalShape_FileOut_handle.write('FIELD_COLORS\t%s\n' % '\t'.join(color_list))
        ExternalShape_FileOut_handle.write('\n# Provide data here\n')
        ExternalShape_FileOut_handle.write('DATA\n')
        for leaf in leaf_matrix_dict:
            ExternalShape_FileOut_handle.write('%s\t%s\n' % (leaf, '\t'.join(leaf_matrix_dict[leaf])))

        ExternalShape_FileOut_handle.close()

    ####################################################################################################################


if __name__ == '__main__':

    iTOL_parser = argparse.ArgumentParser(usage=iTOL_usage)
    iTOL_parser.add_argument('-Labels',             required=False, action='store_true',    help='Labels')
    iTOL_parser.add_argument('-ColoredLabel',       required=False, action='store_true',    help='ColoredLabel')
    iTOL_parser.add_argument('-MultiStyleLabel',    required=False, action='store_true',    help='MultiStyleLabel')
    iTOL_parser.add_argument('-ColorStrip',         required=False, action='store_true',    help='ColorStrip')
    iTOL_parser.add_argument('-ColorRange',         required=False, action='store_true',    help='ColorRange')
    iTOL_parser.add_argument('-ColorClade',         required=False, action='store_true',    help='ColorClade')
    iTOL_parser.add_argument('-ColorLabel',         required=False, action='store_true',    help='ColorLabel')
    iTOL_parser.add_argument('-ColorLeafLabel',     required=False, action='store_true',    help='ColorLabel')
    iTOL_parser.add_argument('-SimpleBar',          required=False, action='store_true',    help='SimpleBar')
    iTOL_parser.add_argument('-Heatmap',            required=False, action='store_true',    help='Heatmap')
    iTOL_parser.add_argument('-ExternalShape',      required=False, action='store_true',    help='ExternalShape')
    iTOL_parser.add_argument('-Binary',             required=False, action='store_true',    help='Binary')
    iTOL_parser.add_argument('-BinaryID',           required=False, action='store_true',    help='Binary specified IDs as 1')
    iTOL_parser.add_argument('-BinaryShape',        required=False, default='2',            help='Binary Shape, choose from 1(rectangle), 2(circle), 3(star), 4, 5 and 6, default is 2')
    iTOL_parser.add_argument('-BinaryColor',        required=False, default='red',          help='Binary Color, default is red')
    iTOL_parser.add_argument('-Connection',         required=False, action='store_true',    help='Connection')
    iTOL_parser.add_argument('-PieChart',           required=False, action='store_true',    help='PieChart')
    iTOL_parser.add_argument('-Collapse',           required=False, action='store_true',    help='Collapse')
    iTOL_parser.add_argument('-id',                 required=False, default=None,           help='File contains leaf id')
    iTOL_parser.add_argument('-ll',                 required=False, default=None,           help='Leaf Label')
    iTOL_parser.add_argument('-lg',                 required=False, default=None,           help='Leaf Group')
    iTOL_parser.add_argument('-gc',                 required=False, default=None,           help='Specify Group/column Color (optional)')
    iTOL_parser.add_argument('-cc',                 required=False, default=None,           help='Specify Column Color (for ExternalShape format) (optional)')
    iTOL_parser.add_argument('-lv',                 required=False, default=None,           help='Leaf Value')
    iTOL_parser.add_argument('-lm',                 required=False, default=None,           help='Leaf Matrix')
    iTOL_parser.add_argument('-dr',                 required=False, default=None,           help='Donor to Recipient')
    iTOL_parser.add_argument('-scale',              required=False, default=None,           help='Scale Values, in format 0-3-6-9')
    iTOL_parser.add_argument('-lt',                 required=False, default=None,           help='Legend Title')
    iTOL_parser.add_argument('-legend',             required=False, action='store_true',    help='show legend for ColorStrip')
    iTOL_parser.add_argument('-show_strip_labels',  required=False, action='store_true',    help='SHOW_STRIP_LABELS')
    iTOL_parser.add_argument('-o',                  required=True,                          help='Output filename')
    args = vars(iTOL_parser.parse_args())
    iTOL(args)
