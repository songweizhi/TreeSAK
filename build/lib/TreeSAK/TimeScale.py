
def get_timescale_file(args):

    # https://github.com/i-c-stratigraphy/chart/blob/main/chart.ttl

    chart_ttl       = args['i']
    time_range      = args['rg']
    time_unit       = args['u']
    interested_rank = args['rk']
    op_dir          = '.'



    chart_ttl       = '/Users/songweizhi/Desktop/chart.ttl'
    time_range      = '3.5-0'
    time_unit       = 'Ga'              # choose from Ma and Ga
    interested_rank = 'Era,Eon'         # choose from Super-Eon, Eon, Era, Period, Sub-Period, Epoch and Age
    op_dir          = '/Users/songweizhi/Desktop'



    time_range_split = time_range.split('-')
    specified_time_range_l = float(time_range_split[0])
    specified_time_range_r = float(time_range_split[1])

    dod = dict()
    current_ns2_id = ''
    current_ns1_rank = ''
    current_range_start = ''
    current_range_end = ''
    current_color = ''
    for each_line in open(chart_ttl):
        if each_line.startswith('ns2:'):
            current_ns2_id = each_line.strip()[len('ns2:'):]
        if each_line.startswith('    ns1:rank '):
            current_ns1_rank = each_line.strip().split('/rank/')[1].split('>')[0]
        if each_line.startswith('    skos:definition "A time period from '):
            if 'million years ago to the present"' in each_line:
                current_range_start = each_line.split('A time period from ')[-1].split()[0]
                current_range_end = '0'
            else:
                current_range_start = each_line.split(' to ')[0].split()[-1]
                current_range_end   = each_line.split(' to ')[1].split()[0]
        if each_line.startswith('    schema:color "'):
            current_color = each_line.strip().split('"')[1]
            if current_ns1_rank not in dod:
                dod[current_ns1_rank] = dict()
            if current_ns2_id not in dod[current_ns1_rank]:
                dod[current_ns1_rank][current_ns2_id] = dict()
            dod[current_ns1_rank][current_ns2_id] = [float(current_range_start), float(current_range_end), current_color]

    interested_rank_list = interested_rank.split(',')
    for each_rank in dod:

        if each_rank in interested_rank_list:
            itol_file = '%s/iTOL_TimeScale_%s.txt' % (op_dir, each_rank)
            itol_file_handle = open(itol_file, 'w')
            itol_file_handle.write('DATASET_TIMESCALE\nSEPARATOR TAB\n')
            itol_file_handle.write('\n')

            itol_file_handle.write('# position related\n')
            itol_file_handle.write('POSITION_ABOVE\t0\n')
            vertical_shift_index = interested_rank_list.index(each_rank) + 1
            vertical_shift_value = vertical_shift_index * 50
            itol_file_handle.write('VERTICAL_SHIFT\t%s\n' % vertical_shift_value)
            itol_file_handle.write('AUTO_SCALE\t0\n')
            itol_file_handle.write('SCALING_FACTOR\t1\n')
            itol_file_handle.write('COVER_TREE\t0\n')
            itol_file_handle.write('\n')

            itol_file_handle.write('# dataset label related\n')
            itol_file_handle.write('DATASET_LABEL\t%s\n' % each_rank)
            itol_file_handle.write('SHOW_LABELS\t1\n')
            itol_file_handle.write('LABEL_POSITION\tend\n')
            itol_file_handle.write('LABEL_SHIFT_X\t30\n')
            itol_file_handle.write('LABEL_SHIFT_Y\t0\n')
            itol_file_handle.write('SIZE_FACTOR\t2\n')
            itol_file_handle.write('LABEL_ROTATION\t0\n')
            itol_file_handle.write('\n')

            itol_file_handle.write('# range label related\n')
            itol_file_handle.write('SHOW_RANGE_LABELS\t1\n')
            itol_file_handle.write('RANGE_LABEL_SIZE\t10\n')
            itol_file_handle.write('RANGE_LABEL_ROTATION\t0\n')
            itol_file_handle.write('RANGE_LABEL_POSITION\tcenter\n')
            itol_file_handle.write('RANGE_LABEL_SHIFT_X\t0\n')
            itol_file_handle.write('RANGE_LABEL_SHIFT_Y\t0\n')
            itol_file_handle.write('RANGE_LABEL_COLOR\t#000000\n')
            itol_file_handle.write('# RANGE_LABEL_AUTO_COLOR\t1\n')
            itol_file_handle.write('\n')

            itol_file_handle.write('# border related\n')
            itol_file_handle.write('BORDER_WIDTH\t0\n')
            itol_file_handle.write('BORDER_COLOR\t#000000\n')
            itol_file_handle.write('BORDER_DASHED\t0\n')
            itol_file_handle.write('\n')

            itol_file_handle.write('# value label related\n')
            itol_file_handle.write('SHOW_VALUE_LABELS\t0\n')
            itol_file_handle.write('VALUE_LABEL_SIZE_FACTOR\t0.5\n')
            itol_file_handle.write('VALUE_LABEL_ROTATION\t0\n')
            itol_file_handle.write('VALUE_LABEL_SHIFT\t0\n')
            itol_file_handle.write('\n')

            itol_file_handle.write('\nDATA\n')
            for each_range in dod[each_rank]:
                current_range_color = dod[each_rank][each_range][2]
                current_range_l = float(dod[each_rank][each_range][0])
                current_range_r = float(dod[each_rank][each_range][1])
                if current_range_l < current_range_r:
                    current_range_l = float(dod[each_rank][each_range][1])
                    current_range_r = float(dod[each_rank][each_range][0])

                if time_unit == 'Ga':
                    current_range_l = current_range_l/1000
                    current_range_r = current_range_r/1000

                range_l_to_write = ''
                range_r_to_write = ''
                if (current_range_l <= specified_time_range_l) and (current_range_r >= specified_time_range_r):
                    range_l_to_write = current_range_l
                    range_r_to_write = current_range_r
                elif (specified_time_range_l >= current_range_l >= specified_time_range_r) and (current_range_r <= specified_time_range_r):
                    range_l_to_write = current_range_l
                    range_r_to_write = specified_time_range_r
                elif (current_range_l >= specified_time_range_l) and (specified_time_range_l >= current_range_r >= specified_time_range_r):
                    range_l_to_write = specified_time_range_l
                    range_r_to_write = current_range_r
                elif (current_range_l <= specified_time_range_r):
                    pass
                elif (current_range_r >= specified_time_range_l):
                    pass
                elif (current_range_l >= specified_time_range_l) and (current_range_r<= specified_time_range_r):
                    range_l_to_write = specified_time_range_l
                    range_r_to_write = specified_time_range_r

                if (range_l_to_write != '') and (range_r_to_write != ''):
                    itol_data_str = '%s\t%s\t\t\t%s\t%s\t%s\tnormal' % (range_l_to_write, range_r_to_write, current_range_color, current_range_color, each_range)
                    itol_file_handle.write(itol_data_str + '\n')
            itol_file_handle.close()

