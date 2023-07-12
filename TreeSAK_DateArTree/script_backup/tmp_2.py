
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
            if len(each_line_split_no_empty) == 9:
                node_id           = each_line_split_no_empty[0]
                value_mean        = each_line_split_no_empty[1]
                value_hpd95_small = each_line_split_no_empty[4]
                value_hpd95_big   = each_line_split_no_empty[5]
                node_to_mean_hpd95_dict[node_id] = [value_mean, value_hpd95_small, value_hpd95_big]
        current_line += 1

    return node_to_mean_hpd95_dict


mcmctree_out_1 = '/Users/songweizhi/Desktop/777/M1_mcmc_txt/M1_PA_75_DeltaLL_100_clock3_nsample500000_out.txt'
mcmctree_out_2 = '/Users/songweizhi/Desktop/777/M1_mcmc_txt/M1_PA_75_DeltaLL_100_clock2_nsample500000_out.txt'

node_to_mean_95_hpd_dict_1 = read_in_posterior_mean(mcmctree_out_1)
node_to_mean_95_hpd_dict_2 = read_in_posterior_mean(mcmctree_out_2)

print('node_to_mean_95_hpd_dict_1\t%s' % node_to_mean_95_hpd_dict_1)
print('node_to_mean_95_hpd_dict_2\t%s' % node_to_mean_95_hpd_dict_2)
