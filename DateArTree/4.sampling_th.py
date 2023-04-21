from sklearn.cluster import DBSCAN
from sklearn.neighbors import NearestNeighbors
import re
from bin.multiple_sbatch import sbatch_all
from api_tools import *
from sklearn import cluster
import pandas as pd
from ete3 import Tree
from subprocess import check_call
from os.path import exists
from tqdm import tqdm
from collections import defaultdict
from dating_workflow.step_script.quick_sampling import *
# from pipelines_repeat20210713.genome_sampling import build_phy


os.chdir('/mnt/ivy/thliao/project/AOB/analysis/20210713_repeat')
genome2completeness = pd.read_csv(
    '/mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/xmoA_genome_tree/quality/micomplete_o/output.completeness', comment='#', sep='\t', index_col=0)

def select_repr_from_clusters(g):
    repr_g = sorted(g, key=lambda x: genome2completeness[x])[-1]
    return repr_g

def clustering(dm_df, gids):
    clus = DBSCAN(min_samples=2, metric='precomputed', eps=0.05)
    clus.fit(dm_df.loc[gids, gids])
    t = list(clus.labels_).count(-1)
    t += len(set([_ for _ in list(clus.labels_) if _ != -1]))
    #print('remaining groups', t)
    return clus.labels_

def get_repr(b):
    labels = clustering(dm_df, list(b))
    # clustering into smaller grounps
    genomes = []
    group2genomes = defaultdict(list)
    for genome, group in zip(b, labels):
        if group == -1:
            genomes.append(genome)
        elif group not in group2genomes:
            group2genomes[group].append(genome)
    for g, _genomes in group2genomes.items():
        genomes.append(select_repr_from_clusters(_genomes))
    return genomes


tax_tab = "/home-user/thliao/.cache/ncbi-genome-download/taxonomy.tab"
tax_df = pd.read_csv(tax_tab, sep="\t", index_col=0)
g2name = tax_df['species'].to_dict()
dm_df = pd.read_csv(
    f"/mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/ref_phy/mash_dist/pairwise_merged_db.tab", sep='\t', index_col=0)

sub_g2name = {_:g2name.get(_.split('.')[0]) for _ in dm_df.index}

nbrs = NearestNeighbors(n_neighbors=20, metric='precomputed').fit(dm_df)
result = nbrs.kneighbors(dm_df.loc[c,:],n_neighbors=10)
genomes_d = {}
for idx,r in enumerate(result[1]):
    dist = result[0][idx]
    genomes_d.update(dict(zip(dm_df.columns[r],dist)))
a = list(set(genomes_d))
for n in set(a).difference(set(c)):
    print(sub_g2name[n],genomes_d[n])



# get the defined AOB/MOB in the concatenated gene tree
gene_tre = Tree(
    "/mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/gene_annotations/xmoCBA_remove_archaea/xmoCBA_renamed.newick", 3)
def get_g(nodes):
    _tmp = []
    for n in [n for n in gene_tre.traverse() if n.name in nodes]:
        _tmp.extend([re.findall("GC[FA]_\d*\.\d", l)[0]
                     for l in n.get_leaf_names()
                     if l.startswith('GC')])
    return _tmp

beta_amo = get_g(['I141_S100'])   # beta
nitro_amo = get_g(['I142_S100'])  # nitrospira
gamma_amo = get_g(['I53_S100'])  # gamma
pmo_lineage = get_g(['I52_S100'])  # gamma
# pxm_lineage = get_g(["I514_S100"]) # gamma
info_df = pd.read_excel('xmoA_containing_genomes.xlsx', index_col=0)
genome2tax = info_df[['class', 'phylum']].to_dict(orient='index')

all_amo = [_ for _ in gamma_amo+beta_amo+nitro_amo if _]
nan_amo = [_ for _ in all_amo if 'nan' == str(genome2tax[_]['class'])] # without enough taxonomic information
gamma_amo = [_ for _ in all_amo if 'Gamma' in str(genome2tax[_]['class'])]
beta_amo = [_ for _ in all_amo if 'Beta' in str(genome2tax[_]['class'])]
nitro_amo = [_ for _ in all_amo if 'Nitrospirae' ==
             str(genome2tax[_]['phylum'])]
k = {'gamma': gamma_amo, 
     'beta': beta_amo, 
     'nitro': nitro_amo}

# resolve the aob without taxonomic information with k-mer based distance matrix 
for nan_aob in nan_amo:
    if nan_aob in dm_df.index:
        min_s = sorted(k, 
                       key=lambda x: dm_df.loc[nan_aob, k[x]].mean())[0]
        x = k[min_s]
        k[min_s].append(nan_aob)

# remove weird genome (otherwise it will affect the sample below)
weird_phylogenetic_positional_geomes = ['GCA_013816325.1',  # in beta but distant
                                        'GCA_015709635.1',  # in gamma but distant
                                        'GCA_015709615.1', 
                                        'GCA_018665985.1',  # not AOB since it lack cycAB, but it embedded in beta
                                        'GCA_007280175.1'  # nitrospira
                                        ]

def remove_weird(list_g):
    return [_ for _ in list_g if _ not in weird_phylogenetic_positional_geomes]

gamma_amo = remove_weird(gamma_amo)
beta_amo += ['GCA_015709635.1', 'GCA_015709615.1']
beta_amo = remove_weird(beta_amo)
nitro_amo = remove_weird(nitro_amo)

final_genomes = []
cyano_gids = open("/mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/cyano_gids.list").read().split('\n')

# get the nearest and most distant functionals genomes 
for a in [beta_amo, nitro_amo, gamma_amo, pmo_lineage]:
    gids = [k for k in set(a) if k in dm_df.index]
    new_dis = {_a: dm_df.loc[_a, cyano_gids].mean() for _a in gids}
    genomes = get_repr(gids)
    final_genomes += list(sorted(genomes, key=lambda x: (new_dis[x],x)))[:3]
    final_genomes += list(sorted(genomes, key=lambda x: (new_dis[x],x)))[-2:]

def get_nodes(gids):
    # get the shallowest node containing all given gids
    _gids = [_ for _ in gids if _ in st.get_leaf_names()]
    last_n = ''
    for n in st.traverse(strategy='postorder'):
        if set(n.get_leaf_names()).intersection(set(_gids)) == set(_gids):
            last_n = n
            # break
            return last_n.get_leaf_names()
        
def get_by_quantile(genomes,genomes_dis):
    got_genomes = []
    for q in [0,0.01,0.03,0.05,0.1,0.2,0.3]:
        # use quantile 
        d = np.quantile(genomes_dis,q)
        g = sorted(genomes,
                   key=lambda x: (abs(dict(zip(genomes,genomes_dis))[x]-d),x) )[0]
        # get the quantile distance and the nearest genome
        got_genomes.append(g)
    return got_genomes

num_neighbours = 5
for a in [beta_amo, nitro_amo, gamma_amo, pmo_lineage]:
    gids = [k for k in set(a) if k in dm_df.index]

    dis, idx = neigh.kneighbors(dm_df.loc[gids, :].mean().values.reshape(1, -1),
                                n_neighbors=len(gids)+500)
    b = {_: _dis
         for _, _dis in zip(dm_df.columns[idx][0], dis[0])}
    genomes = get_repr(list(b))
    genomes = [_ for _ in genomes if  (_ not in gids)]
    same_lineage_genomes = get_nodes(gids)
    genomes = [_ for _ in genomes if _ not in same_lineage_genomes]
    # get genome outside the same lineage
    genomes = list(sorted(genomes, key=lambda x: (b[x],x) ))
    genomes_dis = [b[_] for _ in genomes]
    
    if genome2tax[gids[0]]['phylum']=='Nitrospirae':
        genomes = get_by_quantile(genomes,genomes_dis)
        final_genomes += genomes[:num_neighbours]
    else:
        final_genomes += genomes[:num_neighbours]

# manual curate the genomes
too_shallow_genomes = ['GCF_009828925.1',
                       'GCF_003932755.1',
                       'GCF_000341735.1',
                       'GCA_003994235.2',
                       'GCF_000384075.1',
                       'GCA_000963695.1',
                       'GCF_000372865.1',
                       'GCF_000527095.1',
                       'GCF_000733855.1']
fg = [_ 
      for _ in final_genomes 
      if _ not in too_shallow_genomes] + \
    weird_phylogenetic_positional_geomes

# alpha genomes for mito-based dating analysis
alpha_gids = ['GCA_000014865.1',
 'GCA_002109495.1',
 'GCA_003015145.1',
 'GCA_008189685.1',
 'GCA_000469665.2',
 'GCA_000264455.2',
 'GCA_900107585.1',
 'GCA_006385135.1',
 'GCA_002924445.1',
 'GCA_007197755.1']
    
outgroups_manual = "GCA_002007425.1 GCA_001566965.1 GCA_001187785.1 GCA_001697225.1 GCA_001772005.1 GCA_900119825.1 GCA_003576955.1 GCA_005877815.1".split(' ')

beta_outgroups = ['GCA_006496635.1','GCA_010027455.1']
fg += outgroups_manual + beta_outgroups

fg = remove_weird(fg)
fg += cyano_gids + ['GCA_015709635.1', 'GCA_015709615.1'] + alpha_gids

fg += ['GCA_001803795.1','GCA_001303205.1']
# outside the phylum nitrospiral
extra_gids = """GCA_004421255.1
GCA_000019665.1
GCA_000019965.1
GCA_000020225.1
GCA_003551785.1
GCA_003696005.1
GCA_002781445.1""".split('\n')
fg+=extra_gids

with open('./dating/topology/dating_genomes.list', 'w') as f1:
    f1.write('\n'.join(set([_ for _ in set(fg) if _])))
