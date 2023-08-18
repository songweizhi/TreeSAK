
ALE_usage = '''
========================= ALE example commands =========================

Modules for performing ALE analysis including:
TreeSAK ALE1  ->  Step 1: get gene tree
TreeSAK ALE2  ->  Step 2: run ALE
TreeSAK ALE3  ->  Step 3: parse ALE output
TreeSAK ALE4  ->  Infer ancestral genome

# Example commands
TreeSAK ALE1 -i OrthologousGroups.txt -s combined_d__Archaea_o_rs.faa -p oma -c genome_taxon.txt -m 50 -n 2 -t 6 -jt 3 -f -o ALE1_op_dir
TreeSAK ALE2 -i ALE1_op_dir -s genome_tree_rooted_noEU.treefile -t 10 -f -o ALE2_op_dir -runALE -docker gregmich/alesuite_new
TreeSAK ALE3 -i1 ALE1_op_dir -i2 ALE2_op_dir -c genome_taxon.txt -color ar_phylum_color_code.txt -o ALE3_op_dir_0.3 -fc 0.3 -f -api S1kZZuDHc0d5M7J5vLnUNQ
TreeSAK ALE3 -i1 ALE1_op_dir -i2 ALE2_op_dir -c genome_taxon.txt -color ar_phylum_color_code.txt -o ALE3_op_dir_0.5 -fc 0.5 -f -api S1kZZuDHc0d5M7J5vLnUNQ
TreeSAK ALE3 -i1 ALE1_op_dir -i2 ALE2_op_dir -c genome_taxon.txt -color ar_phylum_color_code.txt -o ALE3_op_dir_0.8 -fc 0.8 -f -api S1kZZuDHc0d5M7J5vLnUNQ
TreeSAK ALE4 -i ALE2_op_dir -o ALE4_op_dir -f -c 0.8

Note: 
Genome names should NOT contain "_".

========================================================================
'''

'''
cd /Users/songweizhi/Desktop/run_ALE_wd

TreeSAK ALE2 -i ALE1_op_dir -s genome_tree_rooted_noEU.treefile -t 10 -f -o ALE2_op_dir -runALE -docker gregmich/alesuite_new
TreeSAK ALE3 -i1 ALE1_op_dir -i2 ALE2_op_dir -c genome_taxon.txt -color ar_phylum_color_code.txt -o ALE3_op_dir_0.3 -fc 0.3 -f -api S1kZZuDHc0d5M7J5vLnUNQ
TreeSAK ALE3 -i1 ALE1_op_dir -i2 ALE2_op_dir -c genome_taxon.txt -color ar_phylum_color_code.txt -o ALE3_op_dir_0.5 -fc 0.5 -f -api S1kZZuDHc0d5M7J5vLnUNQ
TreeSAK ALE3 -i1 ALE1_op_dir -i2 ALE2_op_dir -c genome_taxon.txt -color ar_phylum_color_code.txt -o ALE3_op_dir_0.8 -fc 0.8 -f -api S1kZZuDHc0d5M7J5vLnUNQ
TreeSAK ALE4 -i ALE2_op_dir -o ALE4_op_dir -f -c 0.8

python3 /Users/songweizhi/PycharmProjects/TreeSAK/TreeSAK/ALE3.py -i1 ALE1_op_dir -i2 ALE2_op_dir -c genome_taxon.txt -color ar_phylum_color_code.txt -o ALE3_op_dir_0.3 -fc 0.3 -f -api S1kZZuDHc0d5M7J5vLnUNQ
python3 /Users/songweizhi/PycharmProjects/TreeSAK/TreeSAK/ALE3.py -i1 ALE1_op_dir -i2 ALE2_op_dir -c genome_taxon.txt -color ar_phylum_color_code.txt -o ALE3_op_dir_0.5 -fc 0.5 -f -api S1kZZuDHc0d5M7J5vLnUNQ
python3 /Users/songweizhi/PycharmProjects/TreeSAK/TreeSAK/ALE3.py -i1 ALE1_op_dir -i2 ALE2_op_dir -c genome_taxon.txt -color ar_phylum_color_code.txt -o ALE3_op_dir_0.8 -fc 0.8 -f -api S1kZZuDHc0d5M7J5vLnUNQ

'''