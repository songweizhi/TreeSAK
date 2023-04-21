import os
import time
from os.path import *
from glob import glob
from ete3 import Tree
from Bio import SeqIO
from tqdm import tqdm
from os.path import exists
from collections import defaultdict
#from pack_up_pipelines.dating.formal.final_params import get_gene_name,strategies_params,genes,gene_names


def mkdir_p(dir):
    """make a directory (dir) if it doesn't exist"""
    if not os.path.exists(dir):
        print(dir)
        os.mkdir(dir)


def batch_iter(iter, batch_size):
    # generating batch according batch_size
    iter = list(iter)
    n_iter = []
    batch_d = 0
    for batch_u in range(0, len(iter), batch_size):
        if batch_u != 0:
            n_iter.append(iter[batch_d:batch_u])
        batch_d = batch_u
    n_iter.append(iter[batch_d : len(iter) + 1])
    return n_iter


def sbatch_all(cmds, reset_workdir=False, thread_per_tasks=1, fixed_cluster="", prefix_name="job", batch_size=1):
    zsh_path = os.popen("which zsh").read().strip("\n ")
    job_directory = f"{os.getcwd()}/.job"

    count_ = 0
    for batch_cmds in tqdm(batch_iter(cmds, batch_size)):
        workdir = realpath(batch_cmds[0].split(";")[0].strip().split(" ")[-1])
        job_file = os.path.join(job_directory, f"{prefix_name}{count_}.job")
        with open(job_file, "w") as fh:
            fh.writelines(f"#!{zsh_path}\n")
            fh.writelines(f"#SBATCH --job-name={prefix_name}{count_}\n")
            fh.writelines(f"#SBATCH --cpus-per-task={thread_per_tasks}\n")
            fh.writelines(
                f"#SBATCH --output={job_directory}/{prefix_name}{count_}.out\n"
            )
            if fixed_cluster.startswith('cl00'):
                fh.writelines(f"#SBATCH -w {fixed_cluster} \n")
            elif fixed_cluster.lower() == 'others':
                fh.writelines(f"#SBATCH --exclude=cl007 \n")
                fh.writelines(f"#SBATCH -N 1 \n")

            if reset_workdir:
                fh.writelines(f"#SBATCH --workdir={workdir}\n")
            for cmd in batch_cmds:
                fh.writelines(cmd + "\n")
        os.system("sbatch %s >/dev/null" % job_file)
        time.sleep(0.5)
        count_ += 1

base_dir = '/mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/add_Chromatiaces'
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


# calculating delta-LL
cmds = []
for contree in glob('/mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/add_Chromatiaces/phy_files/*prot/withEuk/*.treefile') + glob('/mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/add_Chromatiaces/phy_files/cog25_prot/*.treefile'):
    gene = contree.split('/')[-1].replace('.treefile','')
    gn = get_gene_name(gene)
    aln = contree.replace('.treefile','.aln')
    genomes = list([_.id for _ in SeqIO.parse(aln,'fasta')])
    t = ''
    st = strategies_params[gn][3]
    st = Tree(st,1)
    st.prune(genomes)
    t+=st.write()+'\n'
    t+=open(contree).read()+'\n'
    tlist = contree.replace('.treefile','.list')
    with open(tlist,'w') as f1:
        f1.write(t)
    odir = dirname(contree)+'/deltaLL'
    if not exists(odir):
        os.makedirs(odir)
    cmd = f"iqtree -nt 10 -m TESTMERGE -redo -mset WAG,LG,JTT,Dayhoff -mrate E,I,G,I+G -mfreq FU -pre {odir}/{gene} -s {aln} -n 0 -zb 10000 -au -z {tlist}"
    if not exists(f"{odir}/{gene}.iqtree"):
        cmds.append(cmd)
sbatch_all(cmds, thread_per_tasks=4, prefix_name='gt')




# gene2num = {}
# gene2dl = {}
# for iqtree in glob('/mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/add_Chromatiaces/phy_files/*prot/withEuk/deltaLL/*.iqtree') + glob('/mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/add_Chromatiaces/phy_files/cog25_prot/deltaLL/*.iqtree'):
#     gene = iqtree.split('/')[-1].replace('.iqtree','')
#     rows = open(iqtree).read().strip().split('\n')
#     idx = [idx for idx,v in enumerate(rows) if 'deltaL  bp-RELL' in v][0]
#     r1, r2 = rows[idx+2],rows[idx+3]
#     r1 = [_ for _ in r1.strip().split(' ') if _]
#     r2 = [_ for _ in r2.strip().split(' ') if _]
#     if r2[2]=='0':
#         gene2dl[gene] = float(r1[2])
#     else:
#         gene2dl[gene] = float(r2[2])
#     gene2num[gene] = len(Tree(iqtree.replace('.iqtree','.treefile')).get_leaf_names())
# # smallest deltaL indicate fewer topological inconsistence
#
# ###! stepwisely removing genes that with lower likelihood
#
# # be careful which part of genes are removed. those with higher/lower values.
#
# from pack_up_pipelines.dating.smaller_genomes.quick_run import gen_pf_list
# cmds = []
# bdir = '/mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/add_Chromatiaces/refined_genes_DeltaL/stepwise'
# for setname,g in zip(gene_names,
#                      genes):
#     if setname != 'G18':continue
#     sortedg = sorted(g,key=lambda x:gene2dl[x])[::-1]  # ascending
#     #print(setname,gene2dl[sortedg[0]],gene2dl[sortedg[-1]])
#     for i in range(0,30):
#         num_removed_g = i*2
#         remaining_genes = sortedg[num_removed_g:]
#         if len(remaining_genes) <= len(sortedg)//3:break
#         odir = f"{bdir}/{setname}/r{num_removed_g}"
#         if not exists(odir):
#             os.makedirs(odir)
#         with open(f'{odir}/R{setname}.list','w') as f1:
#             f1.write('\n'.join(remaining_genes))
#         phy_file = f'{odir}/R{setname}_pf/1pf.phy'
#         if not exists(phy_file):
#             gen_pf_list(f'{odir}/R{setname}.list',
#                     indir=strategies_params[setname][2] ,
#                     pf_range=[1])
#         caltree = strategies_params[setname][1]
#         cmd = f"python3 ~/script/evol_tk/dating_workflow/bin/dating_pro.py -i {phy_file} -it {caltree} -o {bdir}/mcmc_o/B15_R{setname}.r{num_removed_g} -rg '1 30 1' -sf 30 -c 2 -p 1 "
#         if not exists(f"{bdir}/mcmc_o/B15_R{setname}.r{num_removed_g}"):
#             cmds.append(cmd)
#
#         caltree = strategies_params[setname][4]
#         if not caltree:
#             continue
#         calname = caltree.split('/')[-1].replace('.newick','')
#         cmd = f"python3 ~/script/evol_tk/dating_workflow/bin/dating_pro.py -i {phy_file} -it {caltree} -o {bdir}/mcmc_o/{calname}_R{setname}.r{num_removed_g} -rg '1 30 1' -sf 30 -c 2 -p 1 "
#         if not exists(f"{bdir}/mcmc_o/{calname}_R{setname}.r{num_removed_g}"):
#             cmds.append(cmd)
# sbatch_all(cmds,thread_per_tasks=1,prefix_name='eukDL')
#
#
# cmds = []
# bdir = '/mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/add_Chromatiaces/refined_genes_DeltaL/stepwise'
# for setname,g in zip(gene_names,genes):
#     sortedg = sorted(g,key=lambda x:gene2dl[x])[::-1]  # ascending
#     # for i in [0,1,2,3,4]:
#     i = 2
#     num_removed_g = i*2
#     remaining_genes = sortedg[num_removed_g:]
#     odir = f"{bdir}/{setname}/r{num_removed_g}"
#     for p in [1,2,3,4,5]:
#         phy_file = f'{odir}/R{setname}_pf/{p}pf.phy'
#         if not exists(phy_file):
#             gen_pf_list(f'{odir}/R{setname}.list',
#                 indir=strategies_params[setname][2] ,
#                 pf_range=[p])
#     caltree = strategies_params[setname][1]
#     f_odir = f"{bdir}/mcmc_o/B15_R{setname}.1pf.r{num_removed_g}"
#     cmd = f"python3 ~/script/evol_tk/dating_workflow/bin/dating_pro.py -i {phy_file} -it {caltree} -o {f_odir} -rg '1 30 1' -sf 30 -c 2 -p 1 "
#     if not exists(f_odir):
#         cmds.append(cmd)
#     caltree = strategies_params[setname][4]
#     if not caltree:
#         continue
#     calname = caltree.split('/')[-1].replace('.newick','')
#     f_odir = f"{bdir}/mcmc_o/{calname}_R{setname}.1pf.r{num_removed_g}"
#     cmd = f"python3 ~/script/evol_tk/dating_workflow/bin/dating_pro.py -i {phy_file} -it {caltree} -o {f_odir} -rg '1 30 1' -sf 30 -c 2 -p 1 "
#     if not exists(f_odir):
#         cmds.append(cmd)
# sbatch_all(cmds,thread_per_tasks=1,prefix_name='1pf')

