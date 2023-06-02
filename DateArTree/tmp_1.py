import os
import sys
import click
import time
from tqdm import tqdm
from os.path import realpath


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

    count_ = 0
    zsh_path = os.popen("which zsh").read().strip("\n ")

    for batch_cmds in tqdm(batch_iter(cmds, batch_size)):
        workdir = realpath(batch_cmds[0].split(";")[0].strip().split(" ")[-1])
        job_file = f"{prefix_name}{count_}.job"
        with open(job_file, "w") as fh:
            fh.writelines(f"#!{zsh_path}\n")
            fh.writelines(f"#SBATCH --job-name={prefix_name}{count_}\n")
            fh.writelines(f"#SBATCH --cpus-per-task={thread_per_tasks}\n")
            #fh.writelines(f"#SBATCH --output={job_directory}/{prefix_name}{count_}.out\n")

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

