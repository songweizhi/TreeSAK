import os
import glob
import shutil
from Bio import SeqIO
from Bio import AlignIO
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


time_format = '[%Y-%m-%d %H:%M:%S] '


def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


def force_create_folder(folder_to_create):
    if os.path.isdir(folder_to_create):
        shutil.rmtree(folder_to_create, ignore_errors=True)
        if os.path.isdir(folder_to_create):
            shutil.rmtree(folder_to_create, ignore_errors=True)
            if os.path.isdir(folder_to_create):
                shutil.rmtree(folder_to_create, ignore_errors=True)
                if os.path.isdir(folder_to_create):
                    shutil.rmtree(folder_to_create, ignore_errors=True)

    os.mkdir(folder_to_create)


def sep_path_basename_ext(file_in):

    # separate path and file name
    file_path, file_name = os.path.split(file_in)
    if file_path == '':
        file_path = '.'

    # separate file basename and extension
    file_basename, file_ext = os.path.splitext(file_name)

    return file_path, file_basename, file_ext


def get_no_hidden_folder_list(wd):
    folder_list = []
    for each_folder in os.listdir(wd):
        if not each_folder.startswith('.'):
            folder_list.append(each_folder)
    return folder_list


def unique_list_elements(list_input):

    list_output = []
    for each_element in list_input:
        if each_element not in list_output:
            list_output.append(each_element)

    return list_output


def ctg_depth_and_gbk_to_gene_depth(ctg_depth_file, gbk_file, skip_depth_file_header, gene_depth_file_folder):

    gbk_file_path, gbk_file_basename, gbk_file_extension = sep_path_basename_ext(gbk_file)
    pwd_depth_file = '%s/%s.depth' % (gene_depth_file_folder, gbk_file_basename)

    # read in depth
    ctg_depth_dict = {}
    line = 0
    for ctg in open(ctg_depth_file):

        ctg_split = ctg.strip().split('\t')

        if skip_depth_file_header is True:
            if line > 0:
                ctg_depth_dict[ctg_split[0]] = float(ctg_split[1])
        else:
            ctg_depth_dict[ctg_split[0]] = float(ctg_split[1])

        line += 1

    # get gene depth
    gene_depth_file_handle = open(pwd_depth_file, 'w')
    gene_depth_file_handle.write('Gene\tDepth\n')
    for seq_record in SeqIO.parse(gbk_file, 'genbank'):

        seq_id = seq_record.id
        seq_depth = ctg_depth_dict[seq_id]

        for feature in seq_record.features:
            if feature.type == 'CDS':
                gene_id = feature.qualifiers['locus_tag'][0]
                for_out = '%s\t%s\n' % (gene_id, seq_depth)
                gene_depth_file_handle.write(for_out)

    gene_depth_file_handle.close()


def barh_plotter(num_list, label_list, query_seq_num, query_ko_NA, fig_width, fig_height, plot_file):

    fig, ax = plt.subplots()
    fig.set_size_inches(fig_width, fig_height)

    y_pos = range(len(num_list))
    ax.barh(y_pos, num_list, height=0.8, align='center', alpha=0.2, linewidth=0)
    ax.set_yticks([])  # not show yticks
    ax.invert_xaxis()  # line up bar on right
    ax.invert_yaxis()  # put first number on top
    ax.axis('tight')   # remove extra spaces at the top and bottom, equal to: ax.margins(0, 0)
    # ax.margins(0, 0.01) # customize space percentage

    ax.set_xlabel('Number of gene')
    ax.set_title('Query genes number: %s, genes without KO: %s' % (query_seq_num, query_ko_NA))

    ax2 = ax.twinx()
    ax2.set_ylim(ax.get_ylim())
    ax2.set_yticks(y_pos)
    ax2.set_yticklabels(label_list)

    plt.tight_layout()
    plt.savefig(plot_file, dpi=300)
    plt.close()
    plt.clf()


def AnnotateNorm(file_in, skip_header, value_column, Divisor_value, file_out, file_out_header):

    file_out_handle = open(file_out, 'w')
    file_out_handle.write(file_out_header)
    line_num = 0
    for each_line in open(file_in):

        each_line_split = each_line.strip().split('\t')
        value_str = each_line_split[value_column - 1]

        if (skip_header is True and line_num > 0) or (skip_header is False):
            value_pct = float(value_str) * 100 / Divisor_value
            each_line_split[value_column - 1] = str(float("{0:.2f}".format(value_pct)))
            file_out_handle.write('%s\n' % '\t'.join(each_line_split))

        line_num += 1

    file_out_handle.close()


def get_gene_list_TotalDepth(gene_list, gene_to_depth_dict):

    total_depth = 0
    for gene in gene_list:
        gene_depth = gene_to_depth_dict[gene]
        total_depth += gene_depth

    return total_depth


def catfasta2phy(msa_dir, msa_ext, concatenated_msa_phy, partition_file):

    concatenated_msa_fasta = '%s.fasta' % concatenated_msa_phy
    msa_file_re            = '%s/*.%s'  % (msa_dir, msa_ext)
    msa_file_list          = [os.path.basename(file_name) for file_name in glob.glob(msa_file_re)]
    msa_file_list_sorted   = sorted(msa_file_list)

    complete_gnm_set = set()
    for each_msa_file in msa_file_list:
        pwd_msa = '%s/%s' % (msa_dir, each_msa_file)
        for each_seq in SeqIO.parse(pwd_msa, 'fasta'):
            complete_gnm_set.add(each_seq.id)

    complete_gnm_list_sorted = sorted([i for i in complete_gnm_set])

    # initialize concatenated msa dict
    gnm_to_seq_dict = {i: '' for i in complete_gnm_list_sorted}
    msa_len_dict = dict()
    for each_msa_file in msa_file_list_sorted:
        gene_id = each_msa_file.split('.' + msa_ext)[0]

        # read in msa
        current_msa_len = 0
        current_msa_len_set = set()
        pwd_current_msa = '%s/%s' % (msa_dir, each_msa_file)
        current_msa_seq_dict = dict()
        for each_seq in SeqIO.parse(pwd_current_msa, 'fasta'):
            complete_gnm_set.add(each_seq.id)
            current_msa_seq_dict[each_seq.id] = str(each_seq.seq)
            current_msa_len_set.add(len(each_seq.seq))
            current_msa_len = len(each_seq.seq)

        if len(current_msa_len_set) != 1:
            print('Sequences with different length were found in %s, program exited!' % each_msa_file)
            exit()

        msa_len_dict[gene_id] = current_msa_len

        # add sequence to concatenated msa dict
        for each_gnm in complete_gnm_list_sorted:
            msa_seq = current_msa_seq_dict.get(each_gnm, current_msa_len*'-')
            gnm_to_seq_dict[each_gnm] += msa_seq

    # write out concatenated msa
    concatenated_msa_handle = open(concatenated_msa_fasta, 'w')
    for each_gnm in complete_gnm_list_sorted:
        concatenated_msa_handle.write('>%s\n' % each_gnm)
        concatenated_msa_handle.write('%s\n' % gnm_to_seq_dict[each_gnm])
    concatenated_msa_handle.close()

    # write out partition file
    end_pos = 0
    partition_file_handle = open(partition_file, 'w')
    for each_m in msa_file_list_sorted:
        gene_id = each_m.split('.' + msa_ext)[0]
        current_m_len = msa_len_dict[gene_id]
        partition_file_handle.write('%s = %s-%s\n' % (each_m, (end_pos + 1), (end_pos + current_m_len)))
        end_pos += current_m_len
    partition_file_handle.close()

    # convert msa in fasta to phy
    AlignIO.convert(concatenated_msa_fasta, 'fasta', concatenated_msa_phy, 'phylip-relaxed')

