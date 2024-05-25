import argparse


deltall_usage = '''
========================= deltall example commands =========================

TreeSAK deltall -i nohup.out -o DeltaLL_op_summary.txt

# This script was wrote to parse the stdout of deltaLL.rb from Sishuo Wang

============================================================================
'''


def deltall(args):

    deltall_stdout_txt = args['i']
    summary_txt        = args['o']

    deltall_op_dict = dict()
    for each_line in open(deltall_stdout_txt):
        if not ((each_line.startswith('WARNING:')) or (each_line.startswith('awk:'))):
            each_line_split = each_line.strip().split('\t')
            marker_id = each_line_split[0]
            value = float(each_line_split[1])
            if marker_id not in deltall_op_dict:
                deltall_op_dict[marker_id] = [value]
            else:
                deltall_op_dict[marker_id].append(value)

    metric_1_dict = dict()
    metric_2_dict = dict()
    for each_marker in deltall_op_dict:
        metric_1_value = float("{0:.2f}".format(deltall_op_dict[each_marker][0]))
        metric_2_value = float("{0:.2f}".format(deltall_op_dict[each_marker][1]))
        metric_1_dict[each_marker] = metric_1_value
        metric_2_dict[each_marker] = metric_2_value

    metric_1_dict_sorted = {k: v for k, v in sorted(metric_1_dict.items(), key=lambda item: item[1])[::-1]}
    metric_2_dict_sorted = {k: v for k, v in sorted(metric_2_dict.items(), key=lambda item: item[1])}

    metric_1_score_dict = dict()
    metric_1_score = 1
    for each_marker_1 in metric_1_dict_sorted:
        metric_1_score_dict[each_marker_1] = metric_1_score
        metric_1_score += 1

    metric_2_score_dict = dict()
    metric_2_score = 1
    for each_marker_2 in metric_2_dict_sorted:
        metric_2_score_dict[each_marker_2] = metric_2_score
        metric_2_score += 1

    overall_score_dict = dict()
    for each_marker in deltall_op_dict:
        metric_score_1 = metric_1_score_dict[each_marker]
        metric_score_2 = metric_2_score_dict[each_marker]
        metric_score_overall = metric_score_1 + metric_score_2
        overall_score_dict[each_marker] = metric_score_overall

    overall_score_dict_sorted = {k: v for k, v in sorted(overall_score_dict.items(), key=lambda item: item[1])}

    summary_txt_handle = open(summary_txt, 'w')
    summary_txt_handle.write('Marker\tmetric1\tmetric1_score\tmetric2\tmetric2_score\toverall_score\n')
    for each_marker in overall_score_dict_sorted:
        metric_value_1 = metric_1_dict[each_marker]
        metric_value_2 = metric_2_dict[each_marker]
        metric_score_1 = metric_1_score_dict[each_marker]
        metric_score_2 = metric_2_score_dict[each_marker]
        metric_score_overall = overall_score_dict_sorted[each_marker]
        summary_txt_handle.write('%s\t%s\t%s\t%s\t%s\t%s\n' % (each_marker, metric_value_1, metric_score_1, metric_value_2, metric_score_2, metric_score_overall))
    summary_txt_handle.close()


if __name__ == '__main__':

    # initialize the options parser
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', required=True, help='input file (e.g., nohup.out)')
    parser.add_argument('-o', required=True, help='output summary')
    args = vars(parser.parse_args())
    deltall(args)
