"""
Batch upload the tree to itol and output the pdf

details are referred to https://github.com/albertyw/itolapi

and
api of itol http://itol.embl.de/help.cgi#batch


"""

from glob import glob
from itolapi import Itol
import pandas as pd
from tqdm import tqdm
import os

ko2treeID = {}
df = pd.read_csv('/mnt/home-backup/thliao/plancto/gene_GandL/genetrees/target_kos.tsv',sep='\t')
ko2name = dict(zip(df.iloc[:,0],df.iloc[:,2]))
for tpath in tqdm(glob(f"gene_GandL/genetrees/iqtrees/*_MV_rooted.newick")):
    ko = tpath.split('/')[-1].split('_')[0]
    name = ko2name[ko]
    os.system(f"cp -r {tpath} ./tmp.tree ")
    itol_uploader = Itol()
    itol_uploader.params['projectName'] = 'batch_upload_anammox'
    # better to create a project with a unique name.
    itol_uploader.params['APIkey'] = 'S1kZZuDHc0d5M7J5vLnUNQ'  # sine we are the same account, we can use the same APIkey
    itol_uploader.params['treeName'] = f"{name}_MV_rooted"
    itol_uploader.add_file('./tmp.tree')
    # add tree file first
    itol_uploader.add_file('./gene_GandL/genetrees/all_names.txt')
    itol_uploader.add_file('./gene_GandL/genetrees/db_phy_color2.txt')
    itol_uploader.add_file('./gene_GandL/genetrees/Lineage_colorrange.txt')
    # add annotation files
    status = itol_uploader.upload()
    assert status != False
    ko2treeID[ko] = itol_uploader.comm.tree_id
    itol_exporter = itol_uploader.get_itol_export()
    # following parameters are optional, refer to https://itol.embl.de/help.cgi#batchExp
    itol_exporter.set_export_param_value('datasets_visible','0')
    itol_exporter.set_export_param_value('display_mode','2')
    itol_exporter.set_export_param_value('range_mode','2')
    itol_exporter.set_export_param_value('dashed_lines','0')
    itol_exporter.set_export_param_value('format', 'pdf')
    itol_exporter.export(f"/mnt/home-backup/thliao/plancto/gene_GandL/genetrees/iqtree_o_pdf/{name}_MV_rooted.pdf")

# following are some useful methods.
# itol_uploader.comm.upload_output
# # SUCCESS: 1234567890
# itol_uploader.comm.tree_id
# # 1234567890
# itol_uploader.get_webpage()
# # http://itol.embl.de/external.cgi?tree=1234567890&restore_saved=1
