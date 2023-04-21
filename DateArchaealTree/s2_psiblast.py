import os
from Bio import SeqIO


def get_best_hit_dict(blast_results_txt):

    best_hit_dict = dict()
    best_hit_score_dict = dict()
    for each_line in open(blast_results_txt):
        each_line_split = each_line.strip().split('\t')
        query_id = each_line_split[0]
        subject_id = each_line_split[1]
        bit_score = float(each_line_split[11])
        current_best_hit_score = best_hit_score_dict.get(subject_id, 0)
        if bit_score > current_best_hit_score:
            best_hit_dict[subject_id] = query_id
            best_hit_score_dict[subject_id] = bit_score

    return best_hit_dict


##################################################### first trial ######################################################

# wd = '/Users/songweizhi/Desktop/DateArTree/02_identify_marker_gene_1st_trial'
# pwd_combined_protein                = '%s/combined.protein'                                 % wd
# ar_gnm_id_txt                       = '%s/38_ar_gnm_id.txt'                                 % wd
# mito_gnm_id_txt                     = '%s/44_Mito_protein_27_gnm_id.txt'                    % wd
# combined_gnm_id_txt                 = '%s/combined_gnm_id.txt'                              % wd
# blast_op_dir                        = '%s/rpsblast_op_3'                                    % wd
# best_hit_id_by_marker_dir           = '%s/best_hit_by_marker_1_id'                          % wd
# best_hit_seq_by_marker_dir          = '%s/best_hit_by_marker_2_seq'                         % wd
# best_hit_seq_by_marker_dir_renamed  = '%s/best_hit_by_marker_3_seq_renamed'                 % wd
# best_hit_aln_by_marker_dir          = '%s/best_hit_by_marker_4_aln'                         % wd
# best_hit_aln_by_marker_dir_trimmed  = '%s/best_hit_by_marker_5_aln_trimmed'                 % wd
# best_hit_aln_by_marker_dir_trimmed_c= '%s/best_hit_by_marker_5_aln_trimmed_concatenated'    % wd
# gnm_rename_txt                      = '%s/gnm_renamed.txt'                                  % wd
first_trial = False

######################################################### run1 #########################################################

#wd = '/Users/songweizhi/Desktop/DateArTree/02_identify_marker_gene'
wd = '/Users/songweizhi/Desktop/DateArTree/02_identify_marker_gene_e50'
#wd = '/Users/songweizhi/Desktop/DateArTree/02_identify_marker_gene_e30'

pwd_combined_protein                = '%s/d__Archaea_o_rs_133_gnms_plus_27_mito.faa'                % wd
ar_gnm_id_txt                       = '%s/133_ar_gnm_id.txt'                                        % wd
mito_gnm_id_txt                     = '%s/44_Mito_protein_27_gnm_id.txt'                            % wd
combined_gnm_id_txt                 = '%s/combined_gnm_id.txt'                                      % wd
blast_op_dir                        = '%s/d__Archaea_o_rs_133_gnms_plus_27_mito_faa_rpsblast_op'    % wd
best_hit_id_by_marker_dir           = '%s/best_hit_by_marker_1_id'                                  % wd
best_hit_seq_by_marker_dir          = '%s/best_hit_by_marker_2_seq'                                 % wd
best_hit_seq_by_marker_dir_renamed  = '%s/best_hit_by_marker_3_seq_renamed'                         % wd
best_hit_aln_by_marker_dir          = '%s/best_hit_by_marker_4_aln'                                 % wd
best_hit_aln_by_marker_dir_trimmed  = '%s/best_hit_by_marker_5_aln_trimmed'                         % wd
best_hit_aln_by_marker_dir_trimmed_c= '%s/best_hit_by_marker_5_aln_trimmed_concatenated'            % wd
gnm_rename_txt                      = '%s/gnm_renamed.txt'                                          % wd

########################################################################################################################

catfasta2phyml_pl = '/Users/songweizhi/PycharmProjects/Sponge_Hologenome/Scripts/catfasta2phyml.pl'

os.system('mkdir %s' % best_hit_id_by_marker_dir)
os.system('mkdir %s' % best_hit_seq_by_marker_dir)
os.system('mkdir %s' % best_hit_seq_by_marker_dir_renamed)
os.system('mkdir %s' % best_hit_aln_by_marker_dir)
os.system('mkdir %s' % best_hit_aln_by_marker_dir_trimmed)
os.system('mkdir %s' % best_hit_aln_by_marker_dir_trimmed_c)

if first_trial is True:
    gnm_rename_old2new_dict = dict()
    for each_gnm in open(gnm_rename_txt):
        each_gnm_split = each_gnm.strip().split('\t')
        gnm_rename_old2new_dict[each_gnm_split[1]] = each_gnm_split[0]

# get rpsblast command
for each_line in open(combined_gnm_id_txt):
    gnm_id = each_line.strip()
    rpsblast_cmd = 'rpsblast -db /home-user/wzsong/DB/sing/sing -evalue 1e-30 -outfmt 6 -num_threads 16 -query d__Archaea_o_rs_133_gnms_plus_27_mito_faa_files/%s.faa -out d__Archaea_o_rs_133_gnms_plus_27_mito_faa_rpsblast_op/%s_rpsblast.txt' % (gnm_id, gnm_id)
    #print(rpsblast_cmd)

# get best_hit_dict_by_marker
best_hit_to_gnm_dict = dict()
best_hit_dict_by_marker = dict()
for each_line in open(combined_gnm_id_txt):
    gnm_id = each_line.strip()
    pwd_blast_op = '%s/%s_rpsblast.txt' % (blast_op_dir, gnm_id)
    best_hit_dict = get_best_hit_dict(pwd_blast_op)

    # get best_hit_to_gnm_dict
    for each_mk in best_hit_dict:
        gene_id = best_hit_dict[each_mk]
        best_hit_to_gnm_dict[gene_id] = gnm_id

    for each_m in best_hit_dict:
        if each_m not in best_hit_dict_by_marker:
            best_hit_dict_by_marker[each_m] = [best_hit_dict[each_m]]
        else:
            best_hit_dict_by_marker[each_m].append(best_hit_dict[each_m])

print('best_hit_dict_by_marker')
print(best_hit_dict_by_marker)
print('best_hit_dict_by_marker')

print('best_hit_to_gnm_dict')
print(best_hit_to_gnm_dict)
print('best_hit_to_gnm_dict')

# write out best hits and extract their sequences
for each_marker in best_hit_dict_by_marker:
    current_m_hit_list = best_hit_dict_by_marker[each_marker]

    # marker identified from # genomes
    #print('%s\t%s' % (each_marker, len(current_m_hit_list)))

    marker_hits_txt         = ('%s/%s.txt'  % (best_hit_id_by_marker_dir,  each_marker)).replace(':', '')
    marker_hits_seq         = ('%s/%s.fa'   % (best_hit_seq_by_marker_dir, each_marker)).replace(':', '')
    marker_hits_seq_renamed = ('%s/%s.fa'   % (best_hit_seq_by_marker_dir_renamed, each_marker)).replace(':', '')
    marker_hits_aln         = ('%s/%s.aln'  % (best_hit_aln_by_marker_dir, each_marker)).replace(':', '')
    marker_hits_aln_trimmed = ('%s/%s.aln'  % (best_hit_aln_by_marker_dir_trimmed, each_marker)).replace(':', '')

    with open(marker_hits_txt, 'w') as marker_hits_txt_handle:
        marker_hits_txt_handle.write('\n'.join(current_m_hit_list))

    # extract sequences
    biosak_select_seq_cmd = 'BioSAK select_seq -seq %s -id %s -out %s -option 1 -oneline' % (pwd_combined_protein, marker_hits_txt, marker_hits_seq)
    #print(biosak_select_seq_cmd)
    #os.system(biosak_select_seq_cmd)

    # rename sequences
    marker_hits_seq_renamed_handle = open(marker_hits_seq_renamed, 'w')
    for each_seq in SeqIO.parse(marker_hits_seq, 'fasta'):
        seq_id = each_seq.id
        seq_gnm = best_hit_to_gnm_dict[seq_id]
        #seq_gnm_renamed = gnm_rename_old2new_dict[seq_gnm]
        marker_hits_seq_renamed_handle.write('>%s\n' % seq_gnm)
        marker_hits_seq_renamed_handle.write('%s\n' % str(each_seq.seq))
    marker_hits_seq_renamed_handle.close()

    # run mafft-einsi
    mafft_cmd = 'mafft-einsi --quiet %s > %s' % (marker_hits_seq_renamed, marker_hits_aln)
    #print(mafft_cmd)
    #os.system(mafft_cmd)

    # trim msa
    trimal_cmd = 'trimal -in %s -out %s -automated1' % (marker_hits_aln, marker_hits_aln_trimmed)
    #print(trimal_cmd)
    #os.system(trimal_cmd)

# concatenate trimmed msa
catfasta2phyml_cmd = 'perl %s --sequential --concatenate %s/*.aln > %s/concatenated.phy' % (catfasta2phyml_pl, best_hit_aln_by_marker_dir_trimmed, best_hit_aln_by_marker_dir_trimmed_c)
#print(catfasta2phyml_cmd)
#os.system(catfasta2phyml_cmd)
