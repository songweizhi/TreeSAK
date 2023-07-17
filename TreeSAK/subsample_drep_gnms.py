import argparse


sample_drep_gnms_usage = '''
============================= sample_drep_gnms example commands =============================

BioSAK sample_drep_gnms -c Cdb.csv -r rep_gnms.txt -k sponge_gnms.txt -o sampled_gnms.txt

=============================================================================================
'''


def cdb_2_gnm_cluster_file(Cdb_file):
    cluster_to_bin_dict = {}
    for each_bin in open(Cdb_file):
        if not each_bin.startswith('genome,secondary_cluster'):
            each_bin_split = each_bin.strip().split(',')
            bin_id = each_bin_split[0]
            secondary_cluster = each_bin_split[1]
            if secondary_cluster not in cluster_to_bin_dict:
                cluster_to_bin_dict[secondary_cluster] = [bin_id]
            else:
                cluster_to_bin_dict[secondary_cluster].append(bin_id)

    return cluster_to_bin_dict


def sample_drep_gnms(args):

    drep_cdb                     = args['c']
    drep_representative_gnm_txt  = args['r']
    gnm_to_keep_txt              = args['k']
    sampled_gnm_txt              = args['o']

    cluster_to_gnm_dict          = cdb_2_gnm_cluster_file(drep_cdb)
    drep_representative_gnm_list = [i.strip() for i in open(drep_representative_gnm_txt)]
    to_keep_gnm_set              = [i.strip() for i in open(gnm_to_keep_txt)]

    subsampled_gnm_set = set()
    for each_c in cluster_to_gnm_dict:
        cluster_gnms = cluster_to_gnm_dict[each_c]
        if len(cluster_gnms) == 1:
            subsampled_gnm_set.add(cluster_gnms[0])
        else:
            rep_g = ''
            g_in_to_keep_list = []
            for each_g in cluster_gnms:
                if each_g in drep_representative_gnm_list:
                    rep_g = each_g
                if each_g in to_keep_gnm_set:
                    g_in_to_keep_list.append(each_g)

            if len(g_in_to_keep_list) == 0:
                subsampled_gnm_set.add(rep_g)
            else:
                for each_to_keep_g in g_in_to_keep_list:
                    subsampled_gnm_set.add(each_to_keep_g)

    sampled_gnm_list_sorted = sorted([i for i in subsampled_gnm_set])
    sampled_gnm_list_sorted_no_ext = ['.'.join(i.split('.')[:-1]) for i in sampled_gnm_list_sorted]

    with open(sampled_gnm_txt, 'w') as sampled_gnm_txt_handle:
        sampled_gnm_txt_handle.write('\n'.join(sampled_gnm_list_sorted_no_ext))


if __name__ == '__main__':

    SliceMSA_parser = argparse.ArgumentParser()
    SliceMSA_parser.add_argument('-c',      required=True,  help='Cdb.csv')
    SliceMSA_parser.add_argument('-r',      required=True,  help='Id of drep representative genomes, with file extension')
    SliceMSA_parser.add_argument('-k',      required=True,  help='ID of genomes to keep, with file extension')
    SliceMSA_parser.add_argument('-o',      required=True,  help='ID of subsampled genomes')
    args = vars(SliceMSA_parser.parse_args())
    sample_drep_gnms(args)
