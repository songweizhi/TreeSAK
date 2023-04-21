# manually test some topology
from ete3 import Tree
import os
from os.path import *
from subprocess import check_call
from bin.multiple_sbatch import sbatch_all

tree1 = Tree("/mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/dating/topology/mixture_models/deno100/phy/deno100.final_TP1.newick",3)
tmp_newick = '/mnt/ivy/thliao/project/AOB/analysis/update_calibrations/dating/96g/dated_tree/96g_C21/96g_C21.newick'
tree2 = Tree(tmp_newick,3)
n = [_ for _ in tree2.traverse() if _.name == 't_n105'][0]
n_p = n.up
n_p.remove_child(n)
n = [_ for _ in tree1.traverse() if _.name == 'I16_S100'][0]
n_p2=n.up
n_p2.remove_child(n)
n_p.add_child(n)
n = [_ for _ in tree2.traverse() if _.name == 't_n103'][0]
n_p = n.up
n_p.remove_child(n)
n_p2.add_child(n)
n_p_1 = n_p2.up
n_p_1.remove_child(n_p2)
n = [_ for _ in tree2.traverse() if _.name == 't_n101'][0]
n_p = n.up
n_p.remove_child(n)
n_p.add_child(n_p2)
n_p_1.add_child(n_p)

n = [_ for _ in tree1.traverse() if _.name == 'I7_S100'][0]
n_up = n.up
n_up.remove_child(n)
n2 = [_ for _ in tree2.traverse() if _.name == 't_n146'][0]
n2_up = n2.up
n2_up.remove_child(n2)
n2 = n2.copy()
n2.prune('GCA_003696005.1|GCA_001803795.1|GCA_005877815.1'.split('|'))

n3 = [_ for _ in n2.traverse() if _.name == 'GCA_005877815.1'][0]
n3_up = n3.up
n3_up.remove_child(n3)
n3_up.add_child(n)
n_up.add_child(n2)

os.system(f"mkdir -p ./dating/topology/mixture_models/deno100/phy")
intree = './dating/topology/mixture_models/deno100/phy/deno100.final_TP1.newick'
with open(intree,'w') as f1:
    f1.write(tree1.write(format=8))
    
et = tree1.copy()
name = 'deno100'
gl = f"./dating/topology/mixture_models/{name}/phy/{name}.list"
all_ids = list(et.get_leaf_names())
with open(gl,'w') as f1:
    f1.write('\n'.join(set(et.get_leaf_names())))
odir = join("dating/cyano_set/cog25_single", name)
if exists(odir):
    os.system(f"rm -rf {odir}/*")
    
cmd = f"""
python3 /home-user/thliao/script/evol_tk/dating_workflow/step_script/extract_cog25.py -in_p /mnt/home-backup/thliao/NCBI/modified_data/direct_protein_files -in_a /mnt/home-backup/thliao/NCBI/modified_data/cog25_annotate -pd /mnt/home-backup/thliao/NCBI/modified_data/prokka_o -o {odir} -evalue 1e-20 -gl  {gl} -ot nucl 
python3 ~/bin/batch_run/batch_mafft.py -i {odir} -s ffn -o {odir} -f -m mafft -gl {gl} ;
"""
check_call(cmd, shell=1)
phy_file = realpath(f"{odir}/../../trees/phy_files/{name}_cog25_nucl.phy")



t,f = ori_tree_dict['euk_cyano']
all_ids = Tree(t,format=f).get_leaf_names()
t,f = ori_tree_dict['no_euk']
cyano_no_euk = Tree(t,format=f).get_leaf_names()
t,f = ori_tree_dict['no_cyano']
no_cyano_euk = Tree(t,format=f).get_leaf_names()
print(len(all_ids),len(cyano_no_euk),len(no_cyano_euk))

from glob import glob
from Bio import SeqIO
from mito_tk import remove_3rd

# cog25 nucl
bac_nc_odir = odir
euk_nc_odir = f'./dating/merged/euk_cog25/nucl'
gene2records = {}
for ofaa in glob(f'{euk_nc_odir}/*.ffn'):
    gene = ofaa.split('/')[-1].replace('.ffn','')
    gene2records[gene] = list(SeqIO.parse(ofaa,'fasta'))
for ofaa in glob(f'{bac_nc_odir}/*.ffn'):
    gene = ofaa.split('/')[-1].replace('.ffn','')
    gene2records[gene] += list([_ for _ in SeqIO.parse(ofaa,'fasta') if _.id in all_ids])

def output_seqs(genomes,p_odir,gene2records):
    if not exists(p_odir):
        os.makedirs(p_odir)
    for gene,records in gene2records.items():
        records = [_ for _ in records if _.id in genomes]
        # with open(f'{p_odir}/{gene}.ffn','w') as f1:
        #     SeqIO.write(records,f1,'fasta-2line')
        with open(f'{p_odir}/{gene}.ffn','w') as f1:
            SeqIO.write([remove_3rd(_) for _ in records],f1,'fasta-2line')
    with open(p_odir+'/used_genomes.list','w') as f1:
        f1.write('\n'.join(list(genomes)))     
def run_concat(odir,ofile,suffix='.faa'):
    cmds = []
    for faa in glob(join(odir,'*'+suffix)):
        if not exists(faa.replace(suffix,'.aln')):
            cmd = f"ginsi {faa} > {faa.replace(suffix,'.aln')}; trimal -in {faa.replace(suffix,'.aln')} -out {faa.replace(suffix,'.trimal')} -automated1 -keepheader"
            cmds.append(cmd)
    cmds.append(f"""python3 /home-user/thliao/script/evol_tk/dating_workflow/bin/concat_aln.py -i {odir} -o {ofile} -s trimal -simple -gl {odir+'/used_genomes.list'} -ct both -no_graph 
    """)
    return cmds   
output_seqs(all_ids,"./dating/cyano_set/nucl/deno100_full",gene2records)
output_seqs(cyano_no_euk,"./dating/cyano_set/nucl/deno100_no_euk",gene2records)
output_seqs(no_cyano_euk,"./dating/cyano_set/nucl/deno100_no_cyano",gene2records)
c1 = run_concat("./dating/cyano_set/nucl/deno100_full","./dating/cyano_set/nucl/phy/deno100_full.trimal",suffix='.ffn')
c2 = run_concat("./dating/cyano_set/nucl/deno100_no_euk","./dating/cyano_set/nucl/phy/deno100_no_euk.trimal",suffix='.ffn')
c3 = run_concat("./dating/cyano_set/nucl/deno100_no_cyano","./dating/cyano_set/nucl/phy/deno100_no_cyano.trimal",suffix='.ffn')
for _ in [c1,c2,c3]:
    sbatch_all(_[:-1])
for _ in [c1[-1],c2[-1],c3[-1]]:
    check_call(_,shell=1)

# cog25 prot
odir = "dating/cyano_set/cog25_single/deno100_prot"
cmd = f"""
python3 /home-user/thliao/script/evol_tk/dating_workflow/step_script/extract_cog25.py -in_p /mnt/home-backup/thliao/NCBI/modified_data/direct_protein_files -in_a /mnt/home-backup/thliao/NCBI/modified_data/cog25_annotate -pd /mnt/home-backup/thliao/NCBI/modified_data/prokka_o -o {odir} -evalue 1e-20 -gl  {gl} """
check_call(cmd,shell=1)

euk_p_odir = f'./dating/merged/euk_cog25/'
bac_p_odir = odir
gene2records = {}
for ofaa in glob(f'{euk_p_odir}/*.faa'):
    gene = ofaa.split('/')[-1].replace('.faa','')
    gene2records[gene] = list(SeqIO.parse(ofaa,'fasta'))
for ofaa in glob(f'{bac_p_odir}/*.faa'):
    gene = ofaa.split('/')[-1].replace('.faa','')
    gene2records[gene] += list([_ for _ in SeqIO.parse(ofaa,'fasta') if _.id in all_ids])
def output_seqs(genomes,p_odir,gene2records):
    if not exists(p_odir):
        os.makedirs(p_odir)
    for gene,records in gene2records.items():
        records = [_ for _ in records if _.id in genomes]
        with open(f'{p_odir}/{gene}.faa','w') as f1:
            SeqIO.write(records,f1,'fasta-2line')
    with open(p_odir+'/used_genomes.list','w') as f1:
        f1.write('\n'.join(list(genomes)))       
output_seqs(all_ids,"./dating/cyano_set/prot/deno100_full",gene2records)
output_seqs(cyano_no_euk,"./dating/cyano_set/prot/deno100_no_euk",gene2records)
output_seqs(no_cyano_euk,"./dating/cyano_set/prot/deno100_no_cyano",gene2records)
c1 = run_concat("./dating/cyano_set/prot/deno100_full","./dating/cyano_set/prot/phy/deno100_full.trimal",suffix='.faa')
c2 = run_concat("./dating/cyano_set/prot/deno100_no_euk","./dating/cyano_set/prot/phy/deno100_no_euk.trimal",suffix='.faa')
c3 = run_concat("./dating/cyano_set/prot/deno100_no_cyano","./dating/cyano_set/prot/phy/deno100_no_cyano.trimal",suffix='.faa')
for _ in [c1,c2,c3]:
    sbatch_all(_[:-1])
for _ in [c1[-1],c2[-1],c3[-1]]:
    check_call(_,shell=1)

## ! convert into diff topolog part (check carefully before using it)

# fix topology within alpha
intree = f'./dating/topology/mixture_models/deno100/phy/deno100.final_TP1.newick'

def get_topology_without_alpha(st):
    _st = st.copy()
    n_replaced = _st.get_common_ancestor(
        ['GCA_000014865.1', 'GCA_000264455.2'])

    if len(n_replaced.get_leaf_names()) == 10:
        n_p = n_replaced.up   
        n_p_p = n_p.up   # All proteobacteria
        n_p.remove_child(n_replaced)
        return _st, n_p_p, n_p
    else:
        beta_gamma = "GCA_018655245.1,GCA_002356115.1".split(',')
        beta_gamma = _st.get_common_ancestor(beta_gamma)
        for _ in n_replaced.children[::]:
            n_replaced.remove_child(_)
        n_replaced.add_child(beta_gamma)
        return _st, all_pro, n_replaced

euk_reference_tree = '/mnt/home-backup/thliao/AOB/analysis/update_calibrations/mito_dating/phylo/manual_topology/euk.tre'
# copy from '/home-user/sswang/project/Mito/results/euk_tree/euk.tre'. The Fig S2A in wang 2021 NC
# rephrase Porphyra purpurea into Porphyra umbilicalis
# add Polysphondylium pallidum manually
mito_usage = [_ for _ in open('/mnt/home-backup/thliao/AOB/analysis/update_calibrations/mito_dating/euk.list').read().split('\n')]
def get_mito(dataset):
    euk_tree = Tree(euk_reference_tree,8)
    euk_tree.prune(dataset)
    et = earse_name(euk_tree)
    return {"Mito": et}

TP_dict = {"TP1": "((((Holo,other_alpha),Rick),Mito),Magneto);",
           "TP2": "(((Holo,other_alpha),(Rick,Mito)),Magneto);",
           "TP3": "((((Holo,Rick),other_alpha),Mito),Magneto);",
           "TP4": "((((Mito,Rick),Holo),other_alpha),Magneto);"}

base_odir = "./dating/topology"
for tp_name, TP in TP_dict.items():
    st = Tree(intree, 8)
    _st, all_pro, n_replaced = get_topology_without_alpha(st.copy())
    Rickettsiales_lineage = "GCA_008189685.1,GCA_003015145.1".split(',')
    Magneto_lineage = ["GCA_000014865.1", "GCA_002109495.1"]
    Holo = "GCA_000469665.2"
    remaining_alpha = 'GCA_000264455.2,GCA_002924445.1'.split(',')

    Rick_n = st.get_common_ancestor(Rickettsiales_lineage)
    Holo_n = [_ for _ in st.traverse() if _.name == Holo][0]
    Mag_n = st.get_common_ancestor(Magneto_lineage)
    other_alpha_n = st.get_common_ancestor(remaining_alpha)
    m_dict = get_mito(mito_usage)
    nodes_dict = {"Rick": Rick_n,
                  "Holo": Holo_n,
                  "Mito":m_dict['Mito'],
                  "other_alpha": other_alpha_n,
                  "Magneto": Mag_n}
    # TP = "(((Holo,other_alpha),Rick),Magneto);"
    for k, n in nodes_dict.items():
        TP = TP.replace(k, n.write(format=3).strip(';'))
    n_replaced.add_child(Tree(TP, format=3))
    

    with open(f"{base_odir}/deno100_{tp_name}_mito_cyano.newick", 'w') as f1:
        f1.write(_st.write(format=8))
        
    _st.remove_child([_ for _ in _st.children if len(_.get_leaf_names())==31][0])
    with open(f"{base_odir}/deno100_{tp_name}_mito.newick", 'w') as f1:
        f1.write(_st.write(format=8))
        
## !
ori_tree_dict = dict(no_cyano=("./dating/topology/deno100_TP1_mito.newick", 8),
                     euk_cyano=(
                         "./dating/topology/deno100_TP1_mito_cyano.newick", 8),
                     no_euk=("./dating/topology/mixture_models/deno100/phy/deno100.final_TP1.newick", 8))
# mito
phy_files_dict1 = dict(euk_cyano="./dating/mito/prot/phy/deno100_full.phy",
                      no_euk="./dating/mito/prot/phy/deno100_no_euk.phy",
                      no_cyano="./dating/mito/prot/phy/deno100_no_cyano.phy")
## for COG25
phy_files_dict2 = dict(euk_cyano="./dating/cyano_set/prot/phy/deno100_full.phy",
                      no_euk="./dating/cyano_set/prot/phy/deno100_no_euk.phy",
                      no_cyano="./dating/cyano_set/prot/phy/deno100_no_cyano.phy")
## for cog25 + mito24
phy_files_dict3 = dict(euk_cyano="./dating/merged/phy/prot/deno100_full.phy",
                      no_euk="./dating/merged/phy/prot/deno100_no_euk.phy",
                      no_cyano="./dating/merged/phy/prot/deno100_no_cyano.phy")

from collections import defaultdict
def get_cmds(cal_,phy_dict,sub_dir_name):
    # cal_ = 'C8'
    # phy_dict = phy_files_dict3
    # sub_dir_name = 'MERGE'
    cmds = defaultdict(list)
    cal_trees_dir = "./dating/calibration_sets/cal_trees"    
    for cal_f in glob(f'./dating/calibration_sets/{cal_}*.txt'):
        scheme = cal_f.split('/')[-1].replace('.txt','').partition('_')[-1]
        if 'no_' in scheme:
            schemes = [scheme,'euk_cyano']
        else:
            schemes = [scheme]
        for s in schemes:
            phy_file = phy_dict[s]
            in_topo,_format = ori_tree_dict[s]
            o_cal_tree = join(cal_trees_dir,f"deno100_{cal_}{scheme}_{s}.newick")
            cmd = f"format_newick.py add-cal -i {in_topo} -o {o_cal_tree} -c {cal_f} -f {_format}"
            #if not exists(o_cal_tree):
            check_call(cmd,shell=1)
            _odir = f"./dating/sys_testing/{sub_dir_name}/deno100_{cal_}{scheme}_Topo{s}_IR_prot"
            if not exists(join(_odir,'mcmc_for','FigTree.tre')):
                cmd = f"""python3 ~/script/evol_tk/dating_workflow/bin/dating_pro.py -i {phy_file} -it {o_cal_tree} -o {_odir} -rg '1 30 1' -sf 20 -c 2 -p 1 """
                cmds['prot IR'].append(cmd)  
            _odir = f"./dating/sys_testing/{sub_dir_name}/deno100_{cal_}{scheme}_Topo{s}_AR_prot"
            if not exists(join(_odir,'mcmc_for','FigTree.tre')):
                cmd = f"""python3 ~/script/evol_tk/dating_workflow/bin/dating_pro.py -i {phy_file} -it {o_cal_tree} -o {_odir} -rg '1 30 1' -sf 20 -c 3 -p 1 """
                cmds['prot AR'].append(cmd)  
            
            phy_file = phy_file.replace('/prot/','/nucl/')
            if not exists(phy_file):
                print(phy_file)
            _odir = f"./dating/sys_testing/{sub_dir_name}/deno100_{cal_}{scheme}_Topo{s}_IR_nucl"    
            if not exists(join(_odir,'mcmc_for','FigTree.tre')):
                cmd = f"""python3 ~/script/evol_tk/dating_workflow/bin/dating_pro.py -i {phy_file} -it {o_cal_tree} -o {_odir} -rg '1 100 1' -sf 20 -c 2 -p 1 -nucl """
                cmds['nucl IR'].append(cmd)
            _odir = f"./dating/sys_testing/{sub_dir_name}/deno100_{cal_}{scheme}_Topo{s}_AR_nucl" 
            if not exists(join(_odir,'mcmc_for','FigTree.tre')):  
                cmd = f"""python3 ~/script/evol_tk/dating_workflow/bin/dating_pro.py -i {phy_file} -it {o_cal_tree} -o {_odir} -rg '1 100 1' -sf 20 -c 3 -p 1 -nucl """
                cmds['nucl AR'].append(cmd)
    return cmds
sub_cmds = []
cmds = get_cmds('C6',phy_files_dict2,'cog25')
sub_cmds.extend([_ for _ in cmds['nucl IR'] if 'no_euk' in _ and 'euk_cyano' not in _])
cmds = get_cmds('C6',phy_files_dict1,'mito')
sub_cmds.extend([_ for _ in cmds['nucl IR'] if 'no_cyano' in _ and 'euk_cyano' not in _])
cmds = get_cmds('C8',phy_files_dict2,'cog25')
sub_cmds.extend([_ for _ in cmds['nucl IR'] if 'no_euk' in _ and 'euk_cyano' not in _])
cmds = get_cmds('C8',phy_files_dict1,'mito')
sub_cmds.extend([_ for _ in cmds['nucl IR'] if 'no_cyano' in _ and 'euk_cyano' not in _])

sbatch_all(sub_cmds,thread_per_tasks=1,prefix_name=f'd100_N_IR')  
    




from tqdm import tqdm
import pandas as pd
from dating_workflow.bin.parse_mcmc import get_posterior_df, get_node_name_from_log
from dating_workflow.debug_for.draw_infinite_plot import fit_line
import re

def get_summary(mcmc, setting, scheme):
    mcmc_df = get_posterior_df(mcmc, burn_in=0, all_col=False)
    coef, r2 = fit_line(
        x=mcmc_df["Posterior mean time (100 Ma)"].values[:-1],   # remove lnL
        y=mcmc_df["CI_width"].values[:-1]  # remove lnL
    )
    t = get_node_name_from_log(
        f'{dirname(dirname(mcmc))}/prior/nodata_mcmctree.log')
    if t is None:
        return
    node2name = {}
    for gname, group in name2group.items():
        raw_name = "t_n%s" % t.get_common_ancestor(group).name
        node2name[raw_name] = gname
    name2t = {}
    for n, name in node2name.items():
        #cis = mcmc_df.loc[n,'CIs']
        mean_time = mcmc_df.loc[n, 'Posterior mean time (100 Ma)']
        name2t[name] = mean_time
    data = [setting, scheme, r2, coef]+[name2t[_]
                                        for _ in name2group]
    return data


def get_final_df():
    sum_df = pd.DataFrame(columns=['setting', 'scheme', 'r-square', 'coef'])
    name2group = {'gamma-AOB': ['GCA_902810445.1', 'GCA_013697045.1'],
                "beta-AOB": ['GCA_001772005.1', 'GCA_900110495.1'],
                "comammox oldest": ['GCA_013140535.1', 'GCA_011090395.1'],
                "gamma-XOB": ['GCA_902810445.1', 'GCA_002072955.1'],
                }
    for name in name2group:
        sum_df[f'{name}'] = ''
    t = glob('/mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/dating/sys_testing/*/deno100*/mcmc_for/mcmc.txt')
    lines=[]
    for mcmc in tqdm(t):
        name = mcmc.split('/')[-3]
        set_name = mcmc.split('/')[-4]
        topo_scheme = '_'.join(name.split('Topo')[-1].split('_')[:2])
        # format_figtree(mcmc,name,topo_scheme)
        data = get_summary(mcmc,name,set_name)
        if not data:
            continue
        sum_df = sum_df.append(pd.DataFrame([data],
                                            columns=sum_df.columns))
        lines.append(int(os.popen(f'wc -l {mcmc}').read().split(' ')[0].strip()))
    sum_df.loc[:,'nSample'] = lines
    sum_df.insert(1, 'clock model', [_.split('_')[-2] for _ in sum_df['setting']])
    sum_df.insert(1, 'seqtype', [_.split('_')[-1] for _ in sum_df['setting']])
    sum_df.insert(1, 'cal', [_.split('_')[1][:2] for _ in sum_df['setting']])
    sum_df.insert(1, 'cal set', [parse_setting(_)[0] for _ in sum_df['setting']])
    sum_df.insert(1, 'topology', [parse_setting(_)[1] for _ in sum_df['setting']])
    s_df = sum_df.drop('setting',axis=1)
    return s_df

s_df = get_final_df().sort_values(['cal','scheme','seqtype', 'clock model','cal set','topology'])


s_df.loc[ (s_df['clock model']=='IR') & (s_df['seqtype']=='nucl') ,:].sort_values(['cal','scheme','seqtype', 'clock model','cal set','topology'])
