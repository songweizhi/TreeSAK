import os
import glob
import argparse
from ete3 import Tree
import multiprocessing as mp


ALE2_usage = '''
========================= ALE2 example commands =========================

TreeSAK ALE2 -1 ALE1_op_dir -s genome.treefile -t 10 -f -runALE -docker gregmich/alesuite_new -o ALE2_op_dir

Note: 
Genome names should NOT contain "_", the program will tackle this automatically.

# You can try to add this while building the docker images
--platform linux/arm64/v8

=========================================================================
'''


def sep_path_basename_ext(file_in):

    # separate path and file name
    f_path, file_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'

    # separate file basename and extension
    f_base, f_ext = os.path.splitext(file_name)

    return f_path, f_base, f_ext


def subset_tree(tree_file_in, leaves_to_keep_list, tree_file_out):

    input_tree = Tree(tree_file_in)
    subset_tree = input_tree.copy()
    subset_tree.prune(leaves_to_keep_list, preserve_branch_length=True)
    if tree_file_out is None:
        return subset_tree.write()
    else:
        subset_tree.write(outfile=tree_file_out)


def prepare_ale_ip_worker(arg_list):

    ufboot_in   = arg_list[0]
    ufboot_out  = arg_list[1]

    ufboot_out_handle = open(ufboot_out, 'w')
    for each_gene_tree in open(ufboot_in):
        gene_tree_str = each_gene_tree.strip()
        gene_tree_in = Tree(gene_tree_str, format=1)
        for leaf in gene_tree_in:
            leaf_name_split = leaf.name.split('_')
            gnm_id = '_'.join(leaf_name_split[:-1])
            gene_index = leaf_name_split[-1]
            gnm_id_renamed = gnm_id.replace('_', '')
            gene_id_renamed = '%s_%s' % (gnm_id_renamed, gene_index)
            leaf.name = gene_id_renamed
        gene_tree_str_renamed = gene_tree_in.write()
        ufboot_out_handle.write(gene_tree_str_renamed + '\n')
    ufboot_out_handle.close()


def ALE2(args):

    gene_tree_dir           = args['1']
    genome_tree_file_rooted = args['s']
    force_create_ale_wd     = args['f']
    num_threads             = args['t']
    ale_wd                  = args['o']
    run_ale                 = args['runALE']
    docker_image            = args['docker']
    run_ale_cmds_txt        = '%s_cmds.txt'     % ale_wd

    ufboot_file_re   = '%s/*.ufboot' % gene_tree_dir
    ufboot_file_list = glob.glob(ufboot_file_re)
    og_to_process_list = []
    for each_ufboot in ufboot_file_list:
        _, ufboot_base, _ = sep_path_basename_ext(each_ufboot)
        og_to_process_list.append(ufboot_base)

    # define file name
    gnm_tree_no_underscore          = 'genome_tree.newick'
    gnm_tree_leaf_rename_txt        = 'genome_tree_leaf_rename.txt'
    gnm_tree_no_underscore_in_wd    = '%s/%s'                       % (ale_wd, gnm_tree_no_underscore)

    # create ale_wd
    if force_create_ale_wd is True:
        if os.path.isdir(ale_wd) is True:
            os.system('rm -r %s' % ale_wd)
    os.system('mkdir %s' % ale_wd)

    # prepare genome tree for running ALE
    gnm_tree_leaf_rename_txt_handle = open(gnm_tree_leaf_rename_txt, 'w')
    gnm_tree_in = Tree(genome_tree_file_rooted, format=1)
    rename_dict = dict()
    for leaf in gnm_tree_in:
        leaf_name = leaf.name
        leaf_name_new = leaf_name.replace('_', '')
        gnm_tree_leaf_rename_txt_handle.write('%s\t%s\n' % (leaf_name_new, leaf.name))
        leaf.name = leaf_name_new
        rename_dict[leaf_name] = leaf_name_new
    gnm_tree_leaf_rename_txt_handle.close()

    gnm_tree_in.write(outfile=gnm_tree_no_underscore_in_wd)

    # prepare gene tree for running ALE
    run_ale_cmds_txt_handle = open(run_ale_cmds_txt, 'w')
    prepare_ale_ip_worker_arg_lol = []
    ale_cmd_list = []
    for qualified_og in og_to_process_list:
        pwd_gene_tree_ufboot = '%s/%s.ufboot' % (gene_tree_dir, qualified_og)
        if os.path.isfile(pwd_gene_tree_ufboot) is False:
            print('%s not found, please build gene tree first!' % pwd_gene_tree_ufboot)
        else:
            gene_tree_ufboot            = '%s.ufboot'                           % qualified_og
            ALEobserve_stdout           = '%s.ALEobserve.log'                   % qualified_og
            ALEml_undated_stdout        = '%s.ALEml_undated.log'                % qualified_og
            pwd_gene_tree_ufboot_in     = '%s/%s'                               % (gene_tree_dir, gene_tree_ufboot)
            pwd_gene_tree_ufboot_out    = '%s/%s'                               % (ale_wd, gene_tree_ufboot)
            obtain_ale_file_cmd         = 'ALEobserve %s > %s'                  % (gene_tree_ufboot, ALEobserve_stdout)
            reconciliation_cmd          = 'ALEml_undated %s %s.ufboot.ale > %s' % (gnm_tree_no_underscore, qualified_og, ALEml_undated_stdout)
            if docker_image is not None:
                obtain_ale_file_cmd = 'docker run -v $PWD:$PWD -w $PWD %s %s'   % (docker_image, obtain_ale_file_cmd)
                reconciliation_cmd  = 'docker run -v $PWD:$PWD -w $PWD %s %s'   % (docker_image, reconciliation_cmd)

            current_arg_list = [pwd_gene_tree_ufboot_in, pwd_gene_tree_ufboot_out]
            run_ale_cmds_txt_handle.write('%s; %s\n' % (obtain_ale_file_cmd, reconciliation_cmd))
            ale_cmd_list.append('%s; %s\n' % (obtain_ale_file_cmd, reconciliation_cmd))
            prepare_ale_ip_worker_arg_lol.append(current_arg_list)
    run_ale_cmds_txt_handle.close()

    # prepare input files and job script for running ALE with multiprocessing
    print('Preparing files for running ALE with %s cores for %s OGs' % (num_threads, len(prepare_ale_ip_worker_arg_lol)))
    pool = mp.Pool(processes=num_threads)
    pool.map(prepare_ale_ip_worker, prepare_ale_ip_worker_arg_lol)
    pool.close()
    pool.join()

    # run ALE
    if run_ale is True:
        print('running ALE with %s cores for %s OGs' % (num_threads, len(prepare_ale_ip_worker_arg_lol)))
        os.chdir(ale_wd)
        pool = mp.Pool(processes=num_threads)
        pool.map(os.system, ale_cmd_list)
        pool.close()
        pool.join()

    print('Done!')


if __name__ == '__main__':

    ALE2_parser = argparse.ArgumentParser()
    ALE2_parser.add_argument('-1',      required=True,                         help='ALE1 output directory')
    ALE2_parser.add_argument('-s',      required=True,                         help='rooted species tree')
    ALE2_parser.add_argument('-o',      required=True,                         help='output dir, i.e., OMA working directory')
    ALE2_parser.add_argument('-runALE', required=False, action="store_true",   help='run ALE')
    ALE2_parser.add_argument('-docker', required=False, default=None,          help='Docker image, if ALE was installed with Docker, e.g., gregmich/alesuite_new')
    ALE2_parser.add_argument('-f',      required=False, action="store_true",   help='force overwrite')
    ALE2_parser.add_argument('-t',      required=False, type=int, default=6,   help='number of threads, default: 6')
    args = vars(ALE2_parser.parse_args())
    ALE2(args)
