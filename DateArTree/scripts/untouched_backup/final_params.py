import os
base_dir = '/mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/add_Chromatiaces'
phy_dir = f'{base_dir}/phy_files'
idir = f'{base_dir}/tree'

# 1. genome list (with/without euk) 
# 2. tree with calibrations B15 with/without euk
# 3. gene folder
# 4. tree without calibrations
# 5. tree with calibrations B15E6 with/without euk (different E6)
strategies_params = {'G18':(f'{base_dir}/phy_files/GnomezR18_nucl/genome.list',
                            f'{base_dir}/system_dating/cal/G18_B5.newick',
                            f'{base_dir}/phy_files/GnomezR19_prot/withEuk/',
                            f'{base_dir}/tree/bac120_AnimalMito.newick',
                            f'{base_dir}/system_dating/cal/G18_B5A2.newick'),
                     'M24':(f'{base_dir}/phy_files/mitoR22_nucl/genome.list',
                            f'{base_dir}/system_dating/cal/M24_B5.newick',
                            f'{base_dir}/phy_files/mitoR24_prot/withEuk/',
                            f'{base_dir}/tree/bac120_plantMito.newick',
                            f'{base_dir}/system_dating/cal/M24_B5A2.newick'),
                     'P39':(f'{base_dir}/phy_files/P39r3_nucl/genome.list',
                            f'{base_dir}/system_dating/cal/P39_B5.newick',
                            f'{base_dir}/phy_files/P39_prot/withEuk/',
                            f'{base_dir}/tree/bac120_Plastid.newick',
                            f'{base_dir}/system_dating/cal/P39_B5P2.newick'),
                     "COG25":(f'{base_dir}/bacOnly.list',
                              f'{base_dir}/system_dating/cal/COG25_B5.newick',
                              f'{base_dir}/phy_files/cog25_prot',
                              f"{base_dir}/tree/bac120_bacOnly.newick",
                              '')}



from glob import glob
all_genes = glob('/mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/add_Chromatiaces/phy_files/*prot/withEuk/*.aln') + glob('/mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/add_Chromatiaces/phy_files/cog25_prot/*.aln')
all_genes = [_.split('/')[-1].replace('.aln','') for _ in all_genes]
genes = [[_ for _ in all_genes if _.startswith('NP')],
        [_ for _ in all_genes if _.startswith('Mito')],
         [_ for _ in all_genes if _.startswith('2')],]
genes.append(list(set(all_genes).difference(set([_ for v in genes for _ in v]))))  
gene_names = ['P39','M24','COG25','G18']
def get_gene_name(gene):
    for gn,g in zip(gene_names,genes):
        if gene in g:
            return gn
setname2genes = dict(zip(gene_names,genes))
gene2setname = {g:setname for setname,genes in setname2genes.items() for g in genes}
        
sn = 'M24'
intree = strategies_params[sn][3]
o_cal_tree = f'{base_dir}/cal/plantMito_B15E6.newick'
old_cal_f = '/mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/dating/cal_latest/B15E6.txt'
a = open(old_cal_f).read()
a.replace('GCA_014859375.1|GCA_002356115.1','GCA_001007875.1|GCA_001469165.1')
with open(f'{base_dir}/cal/plantMito_B15E6.txt','w') as f1:
    f1.write(a)
cmd = f"format_newick.py add-cal -i {intree} -o {o_cal_tree} -c {base_dir}/cal/plantMito_B15E6.txt -f 1"
os.system(cmd)

sn = 'P39'
intree = strategies_params[sn][3]
o_cal_tree = f'{base_dir}/cal/Plastid_B15E6.newick'
old_cal_f = '/mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/dating/plastids_test/calibrations/B15E6.txt'
a = open(old_cal_f).read()
a.replace('GCA_014859375.1|GCA_002356115.1','GCA_001007875.1|GCA_001469165.1')
with open(f'{base_dir}/cal/Plastid_B15E6.txt','w') as f1:
    f1.write(a)
cmd = f"format_newick.py add-cal -i {intree} -o {o_cal_tree} -c {base_dir}/cal/Plastid_B15E6.txt -f 1"
os.system(cmd)

sn = 'G18'
intree = strategies_params[sn][3]
o_cal_tree = f'{base_dir}/cal/AnimalMito_B15NUC2.newick'
old_cal_f = '/mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/dating/nucleus_set/cal/B15NUCE2.txt'
a = open(old_cal_f).read()
a.replace('GCA_014859375.1|GCA_002356115.1','GCA_001007875.1|GCA_001469165.1')
with open(f'{base_dir}/cal/AnimalMito_B15NUCE2.txt','w') as f1:
    f1.write(a)
cmd = f"format_newick.py add-cal -i {intree} -o {o_cal_tree} -c {base_dir}/cal/AnimalMito_B15NUCE2.txt -f 1"
os.system(cmd)        