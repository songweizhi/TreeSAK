
def read_in_posterior_mean(mcmctree_out):

    # read in Posterior mean
    node_to_mean_hpd95_dict = dict()
    current_line = 1
    posterior_mean_header_line = 0
    for each_line in open(mcmctree_out):
        if 'Posterior mean (95% Equal-tail CI) (95% HPD CI) HPD-CI-width' in each_line:
            posterior_mean_header_line = current_line

        if (posterior_mean_header_line != 0) and (current_line > posterior_mean_header_line):
            each_line_split = each_line.strip().split(' ')

            each_line_split_no_empty = []
            for each_element in each_line_split:
                if each_element not in ['', '(']:
                    each_element_value = each_element.replace('(', '').replace(')', '').replace(',', '')
                    each_line_split_no_empty.append(each_element_value)
            print(each_line_split_no_empty)
            print(len(each_line_split_no_empty))
            if len(each_line_split_no_empty) == 9:
                node_id           = each_line_split_no_empty[0]
                value_mean        = each_line_split_no_empty[1]
                value_hpd95_small = each_line_split_no_empty[4]
                value_hpd95_big   = each_line_split_no_empty[5]
                node_to_mean_hpd95_dict[node_id] = [value_mean, value_hpd95_small, value_hpd95_big]
        current_line += 1

    return node_to_mean_hpd95_dict


node_to_mean_95_hpd_dict = read_in_posterior_mean('/Users/songweizhi/Desktop/999/dating_results/topo1_C_clock3_nsample200000_PhyloHessian_LG_G_run1_out.txt')
print(node_to_mean_95_hpd_dict)

