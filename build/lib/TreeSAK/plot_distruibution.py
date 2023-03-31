import io
import pandas as pd
from tqdm import tqdm
from ete3 import Tree
from glob import glob
from os.path import *
import plotly.express as px
import plotly.graph_objects as go
import plotly.figure_factory as ff


def read_mcmc(mcmc, all_col=False):
    if type(mcmc) != str:
        return mcmc
    if all_col:
        mcmc_df = pd.read_csv(mcmc, sep='\t', index_col=0)
    else:
        f1 = open(mcmc)
        header = [_ for _ in next(f1).strip().split('\t')]
        r_header = [_ for _ in header if not _.startswith('r_g')]
        # normally it need to iterate rows and ignore the columns representing rates
        text = '\t'.join(r_header) + '\n'
        r_header = set(r_header)
        for row in f1:
            text += '\t'.join([r for r, h in zip(row.strip().split('\t'), header) if h in r_header]) + '\n'
        mcmc_df = pd.read_csv(io.StringIO(text), sep='\t', index_col=0)
    return mcmc_df


def get_node_name_from_log(f):
    # f should be the *.log file
    rows = open(f).read().split('\n')
    idx = [_ for _, r in enumerate(rows) if r == 'Species tree']
    if not idx:
        print("prior not complete")
        return
    idx = idx[0]
    start_idx = idx + 3
    end_idx = 0
    for _ in range(idx, 100000):
        if rows[_] == '':
            end_idx = _
            break
    tree_idx1 = end_idx + 1
    tree_idx2 = end_idx + 2
    # find the index
    n2father = {}
    for i in range(start_idx, end_idx):
        row = [_ for _ in rows[i].split(' ') if _]
        father, n, name = row[0], row[1], row[2]
        n2father[name if len(row) == 4 else n] = father

    t = Tree(rows[tree_idx2], format=8)
    for l in t.traverse('postorder'):
        if l.up is None:
            break
        if not l.up.name:
            l.up.name = n2father[l.name]
    return t


indir                   = '/Users/songweizhi/Desktop/DateArTree/plot_distruibution/stepwise'
tree_dir                = '/Users/songweizhi/Desktop/DateArTree/plot_distruibution/treefile_dir'
plot_dir                = '/Users/songweizhi/Desktop/DateArTree/plot_distruibution'

gene_names              = ['M24', 'COG25']
M24_gene_list           = ['MitoCOG0043', 'MitoCOG0040', 'MitoCOG0055', 'MitoCOG0052', 'MitoCOG0053', 'MitoCOG0133', 'MitoCOG0008', 'MitoCOG0009', 'MitoCOG0027', 'MitoCOG0031', 'MitoCOG0030', 'MitoCOG0001', 'MitoCOG0003', 'MitoCOG0012', 'MitoCOG0010', 'MitoCOG0004', 'MitoCOG0005', 'MitoCOG0011', 'MitoCOG0039', 'MitoCOG0060', 'MitoCOG0071', 'MitoCOG0059', 'MitoCOG0067', 'MitoCOG0066']
COG25_gene_list         = ['223163', '223176', '223175', '223607', '223159', '223165', '223170', '223164', '223158', '223172', '223128', '223665', '223275', '223328', '223280', '223127', '223279', '273102', '223130', '223181', '223180', '223168', '223178', '223596', '223556']


setname2genes = dict()
setname2genes['M24']   = M24_gene_list
setname2genes['COG25'] = COG25_gene_list


gene2num = {}
gene2dl = {}
for gene_id in (M24_gene_list + COG25_gene_list):
    pwd_tree_file  = '%s/%s.treefile' % (tree_dir, gene_id)
    pwd_iqtree_log = '%s/%s.iqtree'   % (tree_dir, gene_id)
    rows = open(pwd_iqtree_log).read().strip().split("\n")
    idx = [idx for idx, v in enumerate(rows) if "deltaL  bp-RELL" in v][0]
    r1, r2 = rows[idx + 2], rows[idx + 3]
    r1 = [_ for _ in r1.strip().split(" ") if _]
    r2 = [_ for _ in r2.strip().split(" ") if _]
    if r2[2] == "0":
        gene2dl[gene_id] = float(r1[2])
    else:
        gene2dl[gene_id] = float(r2[2])
    gene2num[gene_id] = len(Tree(pwd_tree_file).get_leaf_names())

# plot 1
for setname, genes in setname2genes.items():
    dl_list = [gene2dl[_] for _ in genes]
    dl_list = sorted(dl_list, reverse=True)
    fig = go.Figure()
    fig.add_bar(y=dl_list)
    fig.update_layout(title_text=setname,title_x=0.5,title_y=1,width=700,height=100,template='simple_white',
                     margin_b=10,margin_l=10,margin_r=10,margin_t=10)
    fig.write_image('%s/Plot_1_%s.pdf' % (plot_dir, setname))

for gene_set in gene_names:
    for _model in ['LG']:  # C60
        t = []
        for f in glob(f'{indir}/{gene_set}/r*/1pf_{_model}/mcmctree/mcmc.txt'):
            if exists(f.replace('mcmc.txt', 'FigTree.tre')):
                t.append((f.split('/')[-4] + ' MCMC', f))
        t = sorted(t, key=lambda x: int(x[0].split(' ')[0][1:]))

        dfs = []
        targets = []
        for cal, mcmc in tqdm(t):
            tre = get_node_name_from_log(mcmc.replace('mcmc.txt','03_mcmctree.log'))
            df  = read_mcmc(mcmc)
            try:
                df = df.sample(5000)
            except:
                print(mcmc)
            for lca, name in [('GCA_001828545.1,GCA_005524015.1', 'Anammox'), ('GCA_013697045.1,GCA_002356115.1', 'Gamma-AOB'),
                              ('GCA_001772005.1,GCA_013521015.1', 'Beta-AOB'), ('GCA_017879665.1,GCA_013140535.1', 'Comammox'),
                              ('Acanthamoeba_castellanii,Andalucia_godoyi', 'Euk'), ('Andalucia_godoyi,Ostreococcus_tauri', 'Euk'),
                              ('Cyanophora_paradoxa,NC_002186.1', 'Euk')]:
                try:
                    n = tre.get_common_ancestor(lca.split(',')).name
                    targets.append(str(n))
                    n = 't_n' + str(n)
                    times = df[[n]]
                except:
                    continue

                times.columns = ['time']
                times.loc[:, 'group name'] = name
                times.loc[:, 'cal'] = cal
                dfs.append(times)

        # plot 2
        _df = pd.concat(dfs, axis=0)
        g2color = {"Gamma-AOB": "#78fce0", "Beta-AOB": "#956bb4", "Comammox": "#edc21a", "Anammox": "#ff8000"}
        _df = _df.loc[_df["group name"].isin(list(g2color)), :]
        _fig = px.violin( _df, y="cal", x="time", color="group name", color_discrete_map=g2color, points=False, orientation="h")
        _fig.update_traces(side="positive", fillcolor='rgba(0,0,0,0)', width=1.8)
        _fig.update_traces(showlegend=False)
        num_y = len(_df["cal"].unique())
        _fig.layout.template = "simple_white"
        _fig.layout.width = 700
        _fig.layout.height = 750
        _fig.update_xaxes(range=[40, 0])
        _fig.update_layout(margin_t=10, title_text=f'{gene_set} {_model}', title_x=0.5)
        _fig.write_image(f'{plot_dir}/Plot_2_{gene_set}_gradient_{_model}.pdf')

        # plot 3
        xs = []
        ys = []
        for ng, subdf in sorted(_df.groupby('cal'),key=lambda x: int(x[0].split(' ')[0].replace('r', ''))):
            t1 = subdf.loc[subdf['group name'] == 'Gamma-AOB', 'time'].median()
            t2 = subdf.loc[subdf['group name'] == 'Anammox', 'time'].median()
            deltaT = t2-t1
            ys.append(deltaT)
            xs.append(int(ng.split(' ')[0].replace('r', '')))
        fig = go.Figure()
        fig.add_scatter(x=xs, y=ys, mode='markers+lines', showlegend=False)
        fig.update_layout(width=300, height=300, margin_t=30, margin_l=10, margin_b=10, margin_r=10,
                          template='simple_white', title_text=f'{gene_set} {_model}', title_x=0.5)
        fig.write_image('%s/Plot_3_%s_%s.pdf' % (plot_dir, gene_set, _model))

