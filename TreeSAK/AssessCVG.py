import io
import argparse
import arviz as az
import pandas as pd
import plotly.graph_objects as go


AssessCVG_usage = '''
================================= AssessCVG example commands =================================

TreeSAK AssessCVG -m1 r1_mcmc.txt -m2 r1_mcmc.txt -o convergence_plot.png

# This script was modified based on the script from Tianhua Liao: 
https://github.com/444thLiao/evol_tk/blob/master/dating_workflow/vis/assess_convergence.py

==============================================================================================
'''


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


def cal_HPD_CI(df,burn_in=2000):
    """
    get HPD CI through mcmc.txt directly instead of reading the log/out file.
    Only calculate high density probility 95%.
    Args:
        df (pd.DataFrame): [description]
        burn_in (int, optional): [description]. Defaults to 2000.
    """
    col2CI = {}
    for colname,col in df.iteritems():
        vals = col.values[burn_in:]
        col2CI[colname] = az.hdi(vals, hdi_prob=.95)
    return col2CI


def get_posterior_df(mcmc, burn_in=2000, scale=1, all_col=True):
    mcmc_df = read_mcmc(mcmc, all_col=all_col)
    if pd.isna(mcmc_df.iloc[-1, -1]):
        # if not completed
        mcmc_df = mcmc_df.drop(mcmc_df.index[-1])
    mcmc_df = mcmc_df.loc[~mcmc_df.isna().any(1), :]
    node_names = [_ for _ in mcmc_df.columns if _.startswith('t_n')]
    rates = [_ for _ in mcmc_df.columns if _.startswith('r_g')]
    paras = [_ for _ in mcmc_df.columns if _.startswith('mu') or _.startswith('sigma2')]

    post_df = pd.DataFrame(columns=['Posterior mean time (100 Ma)',
                                    'CI_width', 'CIs'],
                           index=node_names)
    raw_n2CI = cal_HPD_CI(mcmc_df, burn_in=burn_in)
    if 'lnL' in mcmc_df.columns:
        post_df.loc['lnL', :] = 'NA'
        post_df.loc['lnL', :] = [round(mcmc_df.loc[:, 'lnL'].mean(), 2),
                                 round(raw_n2CI['lnL'][1] - raw_n2CI['lnL'][0], 2),
                                 f"{round(raw_n2CI['lnL'][0], 2)} - {round(raw_n2CI['lnL'][1], 2)}",
                                 ]

    n2CI = {k: f"{round(v[0] * scale, 2)} - {round(v[1] * scale, 2)}"
            for k, v in raw_n2CI.items()}
    n2mean_time = {k: round(v * scale, 2)
                   for k, v in mcmc_df.mean().to_dict().items()}

    post_df.loc[node_names, 'Posterior mean time (100 Ma)'] = [n2mean_time[_]
                                                               for _ in post_df.index
                                                               if _ != 'lnL']
    post_df.loc[node_names, 'CIs'] = [n2CI[_]
                                      for _ in post_df.index
                                      if _ != 'lnL']
    post_df.loc[node_names, 'CI_width'] = [raw_n2CI[_][1] * scale - raw_n2CI[_][0] * scale
                                           for _ in post_df.index
                                           if _ != 'lnL']
    return post_df


def AssessCVG(args):

    mcmc_txt_1  = args['m1']
    mcmc_txt_2  = args['m2']
    output_plot = args['o']

    CI_1 = get_posterior_df(mcmc_txt_1)
    CI_2 = get_posterior_df(mcmc_txt_2)

    # remove lnL row
    CI_1 = CI_1.iloc[:-1, :]
    CI_2 = CI_2.iloc[:-1, :]

    dis1 = list(CI_1['Posterior mean time (100 Ma)'])
    dis2 = list(CI_2['Posterior mean time (100 Ma)'])

    fig = go.Figure()
    fig.add_scatter(x=dis1, y=dis2, name='compared', mode='markers')
    fig.add_scatter(x=[min(dis1 + dis2), max(dis1 + dis2)],
                    y=[min(dis1 + dis2), max(dis1 + dis2)],
                    mode='lines', name='y=x')

    fig.layout.width = 750
    fig.layout.height = 750
    fig.layout.xaxis.title = "run2 posterior mean time"
    fig.layout.yaxis.title = "run1 posterior mean time"
    fig.write_image(output_plot)


if __name__ == '__main__':

    # initialize the options parser
    parser = argparse.ArgumentParser()
    parser.add_argument('-m1',      required=True,   help='mcmc.txt from run 1')
    parser.add_argument('-m2',      required=True,   help='mcmc.txt from run 2')
    parser.add_argument('-o',       required=True,   help='output convergence plot')
    args = vars(parser.parse_args())
    AssessCVG(args)
