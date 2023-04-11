import os
import argparse
import itertools
from ete3 import Tree


tree_file = '/Users/songweizhi/Desktop/test.tree'
tre_sub = Tree(tree_file, format=1, quoted_node_names=True)
print(tre_sub)

