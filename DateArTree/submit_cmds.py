import os
import argparse


def submit_cmds(args):

    js_prefix = args['p']
    cmd_txt   = args['c']
    js_cpu    = args['t']

    js_index = 1
    for each_cmd in open(cmd_txt):
        cmd_to_submit = each_cmd.strip()
        submit_cmd = 'submitHPC.sh --cmd "%s" -n %s -c %s_%s' % (cmd_to_submit, js_cpu, js_prefix, js_index)
        js_index += 1
        os.system(submit_cmd)


if __name__ == '__main__':

    # initialize the options parser
    parser = argparse.ArgumentParser()
    parser.add_argument('-p',   required=False, default='js',        help='job script prefix')
    parser.add_argument('-c',   required=True,                       help='cmds file')
    parser.add_argument('-t',   required=False, default=1, type=int, help='number of threads')
    args = vars(parser.parse_args())
    submit_cmds(args)
