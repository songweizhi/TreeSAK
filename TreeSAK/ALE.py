
ALE_usage = '''
================================================= ALE example commands =================================================

# ALE modules
TreeSAK ALE1         -> Step 1: get gene tree
TreeSAK ALE2         -> Step 2: run ALE
TreeSAK ALE3         -> Step 3: parse ALE outputs (including ancestral genome reconstruction, gene family transfer propensity and verticality)
TreeSAK ALE4         -> Filter ALE identified HGTs
TreeSAK ALE5         -> Get RTC file based on ALE detected HGTs
TreeSAK SingleAleHGT -> Perform HGT analysis using ALE for single protein family
TreeSAK ALE6         -> faa ancestral genomes

# Example commands
TreeSAK ALE1 -i OrthologousGroups.txt -s combined_d__Archaea_o_rs.faa -p oma -m 50 -t 12 -jst 3 -f -o ALE1_op_dir
TreeSAK ALE2 -i ALE1_op_dir -s genome_tree_rooted_noEU.treefile -t 10 -f -o ALE2_op_dir -runALE -docker gregmich/alesuite_new
TreeSAK ALE3 -i ALE2_op_dir -o ALE3_op_dir_c0.75 -f -c 0.75
TreeSAK ALE4 -i1 ALE1_op_dir -i2 ALE2_op_dir -c genome_taxon.txt -color phylum_color.txt -o ALE4_op_dir_0.3 -fc 0.3 -f -api your_own_itol_api
TreeSAK ALE4 -i1 ALE1_op_dir -i2 ALE2_op_dir -c genome_taxon.txt -color phylum_color.txt -o ALE4_op_dir_0.5 -fc 0.5 -f -api your_own_itol_api
TreeSAK ALE4 -i1 ALE1_op_dir -i2 ALE2_op_dir -c genome_taxon.txt -color phylum_color.txt -o ALE4_op_dir_0.8 -fc 0.8 -f -api your_own_itol_api
TreeSAK SingleAleHGT -i OMA00001.aln -s genome.treefile -fc 0.3 -c genome_taxon.txt -color phylum_color.txt -api S1kZZuDHc0d5M7J5vLnUNQ -t 9 -f -o OMA00001_ALE_HGT_wd

Note:
Genome names should NOT contain "_".

========================================================================================================================
'''

'''
cd /Users/songweizhi/Desktop/run_ALE_wd
TreeSAK ALE2 -i ALE1_op_dir -s genome_tree_rooted_noEU.treefile -t 10 -f -o ALE2_op_dir -runALE -docker gregmich/alesuite_new
TreeSAK ALE3 -i ALE2_op_dir -c 0.8 -f -o ALE3_op_dir
TreeSAK ALE4 -i1 ALE1_op_dir -i2 ALE2_op_dir -c genome_taxon.txt -color phylum_color.txt -o ALE4_op_dir_0.3 -fc 0.3 -f -api S1kZZuDHc0d5M7J5vLnUNQ
TreeSAK ALE4 -i1 ALE1_op_dir -i2 ALE2_op_dir -c genome_taxon.txt -color phylum_color.txt -o ALE4_op_dir_0.5 -fc 0.5 -f -api S1kZZuDHc0d5M7J5vLnUNQ
TreeSAK ALE4 -i1 ALE1_op_dir -i2 ALE2_op_dir -c genome_taxon.txt -color phylum_color.txt -o ALE4_op_dir_0.8 -fc 0.8 -f -api S1kZZuDHc0d5M7J5vLnUNQ

python3 /Users/songweizhi/PycharmProjects/TreeSAK/TreeSAK/ALE4.py -i1 ALE1_op_dir -i2 ALE2_op_dir -c genome_taxon.txt -color phylum_color.txt -o ALE4_op_dir_0.3 -fc 0.3 -f -api S1kZZuDHc0d5M7J5vLnUNQ
python3 /Users/songweizhi/PycharmProjects/TreeSAK/TreeSAK/ALE4.py -i1 ALE1_op_dir -i2 ALE2_op_dir -c genome_taxon.txt -color phylum_color.txt -o ALE4_op_dir_0.5 -fc 0.5 -f -api S1kZZuDHc0d5M7J5vLnUNQ
python3 /Users/songweizhi/PycharmProjects/TreeSAK/TreeSAK/ALE4.py -i1 ALE1_op_dir -i2 ALE2_op_dir -c genome_taxon.txt -color phylum_color.txt -o ALE4_op_dir_0.8 -fc 0.8 -f -api S1kZZuDHc0d5M7J5vLnUNQ
python3 /Users/songweizhi/PycharmProjects/TreeSAK/TreeSAK/ALE3.py -i ALE2_op_dir -c 0.8 -f -o ALE3_op_dir
'''

'''
cd /Users/songweizhi/Documents/Research/Sponge_Hologenome/6_ALE_wd
python3 /Users/songweizhi/PycharmProjects/TreeSAK/TreeSAK/ALE1.py -i OMA_op_filtered/OrthologousGroups.txt -s OMA_op_filtered/OrthologousGroups.fasta -p oma -m 3 -t 10 -jt 2 -f -o ALE1_op_dir
TreeSAK ALE2 -i ALE1_op_dir -s genome_tree_rooted_noEU.treefile -t 10 -f -o ALE2_op_dir -runALE -docker gregmich/alesuite_new

cd /Users/songweizhi/Documents/Research/Sponge_Hologenome/8_ALE_wd_all_OGs
TreeSAK ALE2 -i ALE1_op_dir_ufboot -s concatenated_rooted.treefile -t 10 -f -o ALE2_op_dir -runALE -docker gregmich/alesuite_new

cd /home-user/wzsong/tmp
TreeSAK ALE2 -i ALE1_op_dir_ufboot -s concatenated_rooted.treefile -t 32 -f -o ALE2_op_dir -runALE -docker gregmich/alesuite_new

cd /Users/songweizhi/Documents/Research/Sponge_Hologenome/8_ALE_wd_all_OGs
TreeSAK ALE3 -i ALE2_op_dir -o ALE3_op_dir_c0.75 -f -c 0.75

cd /Users/songweizhi/Documents/Research/Sponge_Hologenome/8_ALE_wd_all_OGs
/usr/local/bin/python3.7 /Users/songweizhi/PycharmProjects/TreeSAK/TreeSAK/ALE3.py -i ALE2_op_dir -o ALE3_op_dir_c0.75 -f -c 0.75 -a ALE1_arcog_description.txt
'''
