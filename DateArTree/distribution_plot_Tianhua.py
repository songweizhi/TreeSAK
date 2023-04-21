import io
import glob
import tqdm
import argparse
import arviz as az
from PIL import Image
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from os.path import exists


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
        text = '\t'.join(r_header)+'\n'
        r_header=set(r_header)
        for row in f1:
            text += '\t'.join([r for r,h in zip(row.strip().split('\t'),header) if h in r_header])+'\n'
        mcmc_df = pd.read_csv(io.StringIO(text), sep='\t', index_col=0)
    return mcmc_df


t = []
for f in glob(f'/mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/add_Chromatiaces/refined_genes_DeltaL/stepwise/G18/r*/1pf_C60/mcmctree/mcmc.txt'):
#             if not exists(dirname(f)+'/M24_B5A2.newick'):
#                 print(f)
#                 continue
    if exists(f.replace('mcmc.txt','FigTree.tre')):
        t.append((f.split('/')[-4]+ ' MCMC',f))
t = sorted(t,key=lambda x:int(x[0].split(' ')[0][1:]))

dfs = []
targets = []
for cal,mcmc in tqdm(t):
#     if re.findall(f"_(B[1-9]+)",cal)[0] not in ['B1','B2','B3','B5']:
#         continue
    tre= get_node_name_from_log(mcmc.replace('mcmc.txt','03_mcmctree.log'))
    df = read_mcmc(mcmc)
    try:
        df = df.sample(5000)
    except:
        print(mcmc)
    for lca, name in [('GCA_001828545.1,GCA_005524015.1','Anammox'),
                      ('GCA_013697045.1,GCA_002356115.1','Gamma-AOB'),
                      ('GCA_001772005.1,GCA_013521015.1','Beta-AOB'),
                      ('GCA_017879665.1,GCA_013140535.1','Comammox'),
                      ('Acanthamoeba_castellanii,Andalucia_godoyi','Euk'),   # Gomez19
                      ('Andalucia_godoyi,Ostreococcus_tauri','Euk'),   # mito24
                      ('Cyanophora_paradoxa,NC_002186.1','Euk'),    # plastid
                ]:
        try:
            n = tre.get_common_ancestor(lca.split(',')).name
            targets.append(str(n))
            n = 't_n' + str(n)
            times = df[[n]]
            #print(cal,n,name)
        except:
            #print(n,name)
            continue

        times.columns = ['time']
        times.loc[:,'group name'] = name
        times.loc[:,'cal'] = cal
        dfs.append(times)

_df = pd.concat(dfs,axis=0)
g2color = {
    "Gamma-AOB": "#78fce0",
    "Beta-AOB": "#956bb4",
    "Comammox": "#edc21a",
    "Anammox": "#ff8000",
   # 'Euk':"#66bb6a"
}

_df = _df.loc[_df["group name"].isin(list(g2color)), :]
_fig = px.violin(
    _df,
    y="cal",
    x="time",
    color="group name",
    color_discrete_map=g2color,
    points=False,
    orientation="h",
    # color='calibration sets',
    #showlegend=False
)
_fig.update_traces(
    side="positive",
    fillcolor='rgba(0,0,0,0)',
    width=1.8,
)

_fig.update_traces(showlegend=False)
num_y = len(_df["cal"].unique())
_fig.layout.template = "simple_white"

_fig.layout.width = 700
_fig.layout.height = 750
_fig.update_xaxes(range=[40, 0])
_fig.update_layout(margin_t=10)
_fig.show()

_fig.write_image(f'/mnt/ivy/thliao/project/AOB/analysis/20210713_repeat/add_Chromatiaces/refined_genes_DeltaL/stepwise/G18/G18_gradient_C60.pdf')

xs = []
ys = []
for ng,subdf in sorted(_df.groupby('cal'),key=lambda x:int(x[0].split(' ')[0].replace('r','')) ):
    #subdf = _df.loc[idx,:]
    t1 = subdf.loc[subdf['group name'] == 'Gamma-AOB','time'].median()
    t2 = subdf.loc[subdf['group name'] == 'Anammox','time'].median()
    deltaT = t2-t1
    ys.append(deltaT)
    xs.append(int(ng.split(' ')[0].replace('r','')))
fig = go.Figure()
fig.add_scatter(x=xs,y=ys,mode='markers+lines',showlegend=False)
fig.update_layout(width=300,height=300,margin_t=30,margin_l=10,margin_b=10,margin_r=10,template='simple_white',
                 title_text = f'G18 C60',title_x=0.5)
#display(Image(fig.to_image()))
