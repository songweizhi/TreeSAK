import operator


def keep_highest_rrtc(rrtc_in, rrtc_out):

    rrtc_highest_prob_dict = dict()
    for each_rrtc in open(rrtc_in):
        each_rrtc_split = each_rrtc.strip().split(':')[0].split('\t')
        rrtc_r = each_rrtc.strip().split(':')[0].split('\t')[0]
        rrtc_d = each_rrtc.strip().split(':')[0].split('\t')[1]
        rrtc_v = float(each_rrtc.strip().split(':')[1])
        rrtc_key = '%s___%s' % (rrtc_r, rrtc_d)
        if rrtc_key not in rrtc_highest_prob_dict:
            rrtc_highest_prob_dict[rrtc_key] = rrtc_v
        else:
            if rrtc_v > rrtc_highest_prob_dict[rrtc_key]:
                rrtc_highest_prob_dict[rrtc_key] = rrtc_v

    with open(rrtc_out, 'w') as rrtc_out_handle:
        for each_rrtc in sorted(rrtc_highest_prob_dict.items(), key=operator.itemgetter(1))[::-1]:
            print(each_rrtc)
            rrtc_r = each_rrtc[0].split('___')[0]
            rrtc_d = each_rrtc[0].split('___')[1]
            rrtc_v = each_rrtc[1]
            rrtc_out_handle.write('%s\t%s:%s\n' % (rrtc_r, rrtc_d, rrtc_v))


rrtc_in  = '/Users/songweizhi/Desktop/rrtc.txt'
rrtc_out = '/Users/songweizhi/Desktop/rrtc_out.txt'
keep_highest_rrtc(rrtc_in, rrtc_out)


demo_dict = { 'a': 6, 'b': 2, 'c': 2 }
for each in sorted(demo_dict.items(), key=operator.itemgetter(1))[::-1]:
    print(each[0])
    print(each[1])

