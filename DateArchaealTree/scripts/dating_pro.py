"""
20191202 Liao Tianhua
This script is mainly wrriten for debugging the mcmctree problems.
Of course it could run mainly mcmctree from input tree, alignment, aaRatemodel.

For now, it is designed for both protein and nucleotide.

For better understanding the usage of this script, see usage_dating.md
"""

import io
import multiprocessing as mp
import os
import re
import subprocess
from glob import glob
from os.path import *
from subprocess import check_call
import click
from Bio import SeqIO
from tqdm import tqdm

# __file__ = '/home-user/thliao/script/evolution_relative/dating_workflow/step_script/dating_pro.py'
# template file dir
template_dir = join(dirname(dirname(__file__)), 'ctl_template')
mcmc_ctl = join(template_dir, 'mcmctree.ctl')
codeml_ctl = join(template_dir, 'codeml.ctl')
baseml_ctl = join(template_dir, 'baseml.ctl')
aaRatefile = join(template_dir, 'lg.dat')


def env_exe(name):
    p = os.environ.get('PATH')
    p = p.split(':')
    for _p in p:
        if exists(join(_p,name)):
            return join(_p,name)
    return 'None'
    #
    # bin_dir = dirname(sys.executable)
    # f = join(bin_dir, name)
    # if exists(f):
    #     return f
    # f = popen(f'which {name} 2>&1 ').read().strip('\n')
    # return f


paml_bin = dirname(env_exe("mcmctree"))


def modify(file, **kwargs):
    text = open(file).read()
    text = text.split('\n')
    new_text = []
    for row in text:
        key = row.split('=')[0].strip()
        if key in kwargs:
            new_text.append(f"{key} = {kwargs[key]}")
        else:
            new_text.append(row)
    return '\n'.join(new_text)


def run(args):
    if isinstance(args, str):
        cmd = args
        log = '/dev/null'
    else:
        cmd, log = args
    try:
        # subprocess.check_output(cmd, shell=1)
        check_call(cmd,
                   shell=True,
                   stdout=open(log, 'w'))

    except subprocess.CalledProcessError as e:
        pass
        # print('error', e.output)
    if log != '/dev/null':
        t = open(log, 'r', newline='\n').read().replace('\r', '\n')
        with open(log, 'w') as f1:
            f1.write(t)


def get_num_phy_file(in_phyfile):
    ndata = 0
    for row in open(in_phyfile):
        row = row.strip('\n').strip().split(' ')
        row = [_ for _ in row if _]
        if all([_.isnumeric() for _ in row]):
            ndata += 1
    return ndata


def generate_tmp(in_phyfile, in_treefile, odir, ndata, template_ctl=mcmc_ctl, use_nucl=False):
    # template_ctl_01 = './01_mcmctree.ctl'
    if not exists(odir):
        os.makedirs(odir)
    new_01_ctl = join(odir, '01_mcmctree_modify.ctl')
    params = {'ndata': ndata,
              'seqfile': in_phyfile,
              'treefile': in_treefile,
              'outfile': './01_out'}
    if use_nucl:
        params['seqtype'] = '0'
    else:
        params['seqtype'] = '2'
    text = modify(template_ctl, **params)
    with open(new_01_ctl, 'w') as f1:
        f1.write(text)
    run(f"export PATH=''; cd {odir}; {paml_bin}/mcmctree 01_mcmctree_modify.ctl 2>&1")

def rename_tmp(ali_dir, ctl_file):
    phy_f = ctl_file.replace('.ctl', '.txt')
    rows = open(phy_f).read().split('\n')
    rows = [_ for _ in rows if _.startswith('GC')]
    rows = ['>' + re.sub('\ {2,}', '\n', row).replace(' ', '') for row in rows]
    now_f = list(SeqIO.parse(io.StringIO('\n'.join(rows)), format='fasta'))
    non_empty_seqs = len([_ for _ in now_f if set(_.seq) != {'-'}])
    match_n = []
    for f in glob(join(ali_dir, '*.trimal')):
        name = basename(f).replace('.trimal', '')
        tmp_f = list(SeqIO.parse(f, format='fasta'))
        _non_empty_seqs = len([_ for _ in tmp_f if set(_.seq) != {'-'}])
        if len(tmp_f[0].seq) == len(now_f[0].seq) and non_empty_seqs == _non_empty_seqs:
            match_n.append(name)
    if len(match_n) == 1:
        return match_n[0]
    else:
        raise Exception(f'multiple...or non,for {ctl_file},have {len(match_n)} file')


def collecting_tmp(tmp_indir, odir, ali_dir=None):
    if not exists(odir):
        os.makedirs(odir)
    ctls = glob(join(tmp_indir, 'tmp0*.ctl'))
    for ctl in ctls:
        name = basename(ctl).replace('.ctl', '')
        if ali_dir is not None:
            # TODO: modify it into it original name/ gene name
            # hard to find a one....
            new_name = name
        else:
            new_name = name

        os.makedirs(join(odir, f'{new_name}'), exist_ok=1)
        os.system(f'mv {tmp_indir}/{name}.* {odir}/{new_name}/')


def run_each_tmp(tmp_indir, odir, aaRatefile=aaRatefile, extra_cmd=None, use_nucl=False):
    if not exists(odir):
        os.makedirs(odir)
    params = []
    ctls = glob(join(tmp_indir, '*', 'tmp0*.ctl'))
    ctls = [_ for _ in ctls if 'modify' not in _]  # remove modify
    for ctl in ctls:
        _params = {'model': 2,
                   'aaRatefile': abspath(aaRatefile),
                   'fix_alpha': 0,
                   'alpha': 0.5,
                   'ncatG': 4}
        if use_nucl:
            _params.pop('aaRatefile')
        new_text = modify(ctl,
                          **_params)
        new_file = ctl.replace('.ctl', '.modify.ctl')
        with open(new_file, 'w') as f1:
            f1.write(new_text)
        if use_nucl:
            bin_ml = f"{paml_bin}/baseml"
        else:
            bin_ml = f"{paml_bin}/codeml"
        params.append((
            f"cd {dirname(new_file)}; {bin_ml} {basename(new_file)}",
            new_file.replace('.modify.ctl', '.log')))

    if extra_cmd is not None:
        params.extend(extra_cmd)
    with mp.Pool(processes=30) as tp:
        _ = list(tqdm((tp.imap(run, params)), total=len(params)))

    # cat all rst2 into in.BV
    rst2_list = glob(join(tmp_indir, '*', 'rst2'))
    text = ''
    for rst2 in sorted(rst2_list,
                       key=lambda x: int(basename(dirname(x)).replace('tmp', ''))):
        text += open(rst2).read()
    with open(join(odir, 'in.BV'), 'w') as f1:
        f1.write(text)
    os.system(f"cp {join(odir, 'in.BV')} {tmp_indir}/out.BV")

default_params = {'seqfile': '',
             'treefile': '',
             'ndata': '',
             'model':0,
             'seqtype': 0,
             'usedata': "2 in.BV 1",
             'outfile': './03_mcmctree.out',
             'clock': 2,
             'BDparas': "1 1 0.1",
             'rgene_gamma': '1 35 1',
             'sigma2_gamma': '1 10 1',
             'burnin': '2000',
             'sampfreq': '2',
             'nsample': '20000',
             'alpha': 0.5}

def final_mcmctree(inBV, in_phyfile, in_treefile, odir, ndata, template_ctl=mcmc_ctl, params_dict={},
                   use_nucl=False,extra_cmd=None,
                   model='0'):
    # for final mcmctree
    if not exists(odir):
        os.makedirs(odir)
    ## default parameters, it will be overwritten by the `main` function
    bd_paras = '1 1 0.1'
    rgene_gamma = '1 35 1'
    sigma2_gamma = '1 10 1'
    burnin = '2000'
    sampfreq = '2'
    nsample = '20000'
    seqfile_b = in_phyfile
    treefile_b = in_treefile
    outfile = './03_mcmctree.out'
    seqtype = 0 if use_nucl else 2
    clock = 2
    param = {'seqfile': seqfile_b,
             'treefile': treefile_b,
             'ndata': ndata,
             'model':str(model),
             'seqtype': seqtype,
             'usedata': "2 in.BV 1",
             'outfile': outfile,
             'clock': clock,
             'BDparas': bd_paras,
             'rgene_gamma': rgene_gamma,
             'sigma2_gamma': sigma2_gamma,
             'burnin': burnin,
             'sampfreq': sampfreq,
             'nsample': nsample,
             'print':'1',
             'alpha': 0.5}
    if params_dict:
        param.update(params_dict)
    text = modify(template_ctl,
                  **param)
    if not exists(f'{odir}/in.BV'):
        os.system(f'cp {inBV} {odir}')
    ofile = join(odir, '03_mcmctree.ctl')
    with open(ofile, 'w') as f1:
        f1.write(text)
    tqdm.write("start running the final mcmctree. ")

    params = [(f"cd {dirname(ofile)}; {paml_bin}/mcmctree 03_mcmctree.ctl 2>&1",
              ofile.replace('.ctl', '.log'))]

    if extra_cmd is not None:
        params.extend(extra_cmd)

    while 1:
        os.system(params[0][0]+' > ' +ofile.replace('.ctl', '.log'))
        if exists(join(dirname(ofile),'FigTree.tre')):
            break


def run_nodata_prior(in_phyfile, in_treefile, odir, ndata, template_ctl=mcmc_ctl, params_dict={}, use_nucl=False):
    if not exists(odir):
        os.makedirs(odir)
    bd_paras = '1 1 0.1'
    rgene_gamma = '1 35 1'
    sigma2_gamma = '1 10 1'
    burnin = '2000'
    sampfreq = '2'
    nsample = '20000'
    seqfile_b = in_phyfile
    treefile_b = in_treefile
    outfile = './nodata_mcmctree.out'
    seqtype = 0 if use_nucl else 2
    clock = 2
    param = {'seqfile': seqfile_b,
             'treefile': treefile_b,
             'ndata': ndata,
             'seqtype': seqtype,
             'usedata': "0",
             'outfile': outfile,
             'clock': clock,
             'BDparas': bd_paras,
             'rgene_gamma': rgene_gamma,
             'sigma2_gamma': sigma2_gamma,
             'burnin': burnin,
             'sampfreq': sampfreq,
             'nsample': nsample,
             'alpha': 0.5}
    if params_dict:
        param.update(params_dict)
    text = modify(template_ctl,
                  **param)
    ofile = join(odir, 'nodata_mcmctree.ctl')
    with open(ofile, 'w') as f1:
        f1.write(text)

    return (f"cd {dirname(ofile)}; {paml_bin}/mcmctree nodata_mcmctree.ctl ",
            ofile.replace('.ctl', '.log'))


def main(in_phyfile, in_treefile, total_odir,
         use_nucl=False, ali_dir=None,
         run_tmp=True, run_prior_only=True, params_dict={}):
    if not exists(total_odir):
        os.makedirs(total_odir)
    ndata = get_num_phy_file(in_phyfile)
    nodata_dir = join(total_odir, 'prior')
    mcmc_for_dir = join(total_odir, 'mcmc_for')
    tmp_odir = join(total_odir, 'tmp_files')
    for f in [nodata_dir,mcmc_for_dir,tmp_odir]:
        if not exists(f): os.makedirs(f)
    prior_cmd = run_nodata_prior(in_phyfile=in_phyfile,
                                 in_treefile=in_treefile,
                                 odir=nodata_dir,
                                 ndata=ndata,
                                 params_dict=params_dict,
                                 use_nucl=use_nucl)
    if run_prior_only:
        run(prior_cmd)
        return
    if run_tmp and not isinstance(run_tmp, str):
        generate_tmp(in_phyfile,
                     in_treefile,
                     tmp_odir,
                     ndata, use_nucl=use_nucl)
        collecting_tmp(tmp_odir,
                       tmp_odir,
                       ali_dir=ali_dir)

        run_each_tmp(tmp_odir,
                     mcmc_for_dir,
                    #  extra_cmd=[prior_cmd],
                     use_nucl=use_nucl)
    elif isinstance(run_tmp, str):
        os.system(f"cp -rf {run_tmp}/* {tmp_odir}")
        os.system(f"cp -rf {run_tmp}/out.BV {mcmc_for_dir}/in.BV")
        # collecting_tmp(run_tmp,
        #                tmp_odir,
        #                ali_dir=ali_dir)
    final_mcmctree(inBV=join(mcmc_for_dir, 'in.BV'),
                   in_phyfile=in_phyfile,
                   in_treefile=in_treefile,
                   odir=mcmc_for_dir,
                   ndata=ndata,
                   params_dict=params_dict,
                   use_nucl=use_nucl,
                   extra_cmd=[prior_cmd]
                   )


def process_path(path):
    if not '/' in path:
        path = './' + path
    if path.startswith('~'):
        path = expanduser(path)
    if path.startswith('.'):
        path = abspath(path)
    return path


def change_parameters(mcmc_for_dir,odir=None,**kwargs):
    if odir is None:
        parent_dir = dirname(realpath(mcmc_for_dir))
        odir = join(parent_dir,'new_odir')
    if not exists(odir):
        os.makedirs(odir)
    cmd = f'ln -sf `realpath {mcmc_for_dir}/in.BV` {odir}/'
    check_call(cmd,shell=1)
    existed_mcmc = glob(join(mcmc_for_dir,'*.ctl'))[0]
    text = modify(existed_mcmc,
                  **kwargs)
    ofile = join(odir, 'mcmctree.ctl')
    with open(ofile, 'w') as f1:
        f1.write(text)

    p = realpath(dirname(ofile))
    return (f"cd {p}; {paml_bin}/mcmctree mcmctree.ctl ",
            ofile.replace('.ctl', '.log'))

@click.command()
@click.option('-i', '--in_phy', 'in_phyfile',help='phylip format alignment file')
@click.option('-it', '--in_tree', 'in_treefile',help='newick format tree with calibration information')
@click.option('-id', '--in_ali_dir', 'in_ali_dir')
@click.option('-o', 'odir',help='')
@click.option('-nucl', 'use_nucl', is_flag=True, default=False)
@click.option('-no_tmp', 'run_tmp', default=True,help='Provided a pre-existed `tmp_files` directory which can avoid generating `in.BV` again. ')
@click.option('-only_prior', 'only_prior', is_flag=True, default=False,help="Not run posterior estimating step")
@click.option('-sf', 'sampfreq', default='20',help="sample frequency  [20]")
@click.option('-p', 'print_f', default='2',help="verbose of print  [2]")
@click.option('-rg', 'rgene_gamma', default='1 35 1',help="rgene_gamma: prior on mutation rate   [1 35 1]")
@click.option('-sg', 'sigma2_gamma', default='1 10 1',help="sigma2_gamma: shape and scale parameters  [1 10 1]")
@click.option('-bd', 'bdparse', default='1 1 0.1',help="verbose of print  [2]")
@click.option('-c', 'clock', default='2',help="2 indicate using IR clock model, while 3 denote AR clock model")
def cli(in_phyfile, in_treefile, in_ali_dir,
        odir, use_nucl,
        run_tmp, only_prior, sampfreq,
        print_f, rgene_gamma, sigma2_gamma,
        bdparse,
        clock):
    in_phyfile = process_path(in_phyfile)
    in_treefile = process_path(in_treefile)
    params_dict = {'sampfreq': str(sampfreq),
                   'print': str(print_f),
                   'rgene_gamma': rgene_gamma,
                   'sigma2_gamma': sigma2_gamma,
                   "BDparas":bdparse,
                   'clock': clock}
    main(in_phyfile, in_treefile,
         use_nucl=use_nucl,
         ali_dir=in_ali_dir,
         total_odir=odir, run_tmp=run_tmp,
         run_prior_only=only_prior, params_dict=params_dict)


if __name__ == "__main__":
    cli()
