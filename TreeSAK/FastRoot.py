import os
import re
import sys
import numpy
import cvxopt
import logging
import argparse
from numpy import *
from treeswift import *
from sys import stdin, stdout, argv, exit, stderr


########################################################################################################################

PROGRAM_NAME = "FastRoot"
PROGRAM_AUTHOR = ["Uyen Mai", "Merve Kilic", "Erfan Sayyari", "Siavash Mirarab"]
PROGRAM_LICENSE = "MIT License"
PROGRAM_VERSION = "1.5"
PROGRAM_YEAR = "2017"
PROGRAM_INSTITUTE = "University of California at San Diego"


def new_logger(myName, myLevel=logging.INFO, myStream=stdout):
    logger = logging.getLogger(myName)
    logger.setLevel(myLevel)
    handler = logging.StreamHandler(myStream)
    formatter = logging.Formatter('%(levelname)s:%(name)s:%(message)s')
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    logger.propagate = False

    return logger


########################################################################################################################

'''
logger = logging.getLogger("Tree_extend.py")
logger.setLevel(logging.INFO)
handler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(levelname)s:%(name)s:%(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.propagate = False
'''


class Tree_extend(object):
    def __init__(self, ddpTree=None, tree_file=None, schema="newick"):  # ,logger_id=1,logger_stream=sys.stderr):
        # self.logger = new_logger(__name__+ "_" + str(logger_id),myStream=logger_stream)
        if tree_file:
            self.ddpTree = read_tree(tree_file, schema)
        else:
            self.ddpTree = ddpTree

    def Bottomup_label(self):
        # assign each node a label so that we can later relate to it
        i = 0
        for node in self.ddpTree.traverse_postorder():
            if node.is_leaf():
                node.name = 'L' + str(i)
            else:
                node.name = 'I' + str(i)
            i += 1

    def Topdown_label(self, label_type="all"):
        # assign each node a label so that we can later relate to it
        i = 0

        for node in self.ddpTree.traverse_preorder():
            if node.is_leaf():
                if label_type == "all" or label_type == "leaves":
                    node.name = 'L' + str(i)
                else:
                    node.name = node.label
            else:
                if label_type == "all" or label_type == "internal":
                    node.name = 'I' + str(i)
                else:
                    node.name = node.label
            i += 1

    def Bottomup_update(self):
        for node in self.ddpTree.traverse_postorder():
            self.Node_init(node)
            self.bUp_update(node)

    def Topdown_update(self):
        for node in self.ddpTree.traverse_preorder():
            self.tDown_update(node, self.Opt_function)

    def compute_distances(self):
        D = {}

        def __compute_dRoot__(node, cumm_l):
            if node.is_leaf():
                D[node.name] = cumm_l
            else:
                for child in node.child_nodes():
                    __compute_dRoot__(child, cumm_l + child.edge_length)

        __compute_dRoot__(self.ddpTree.root, 0)
        return D

    def compute_ingroup_distances(self):
        D = []

        def __compute_dLeaf__(node, cumm_l):
            if node.is_leaf():
                D.append(cumm_l)
            else:
                for child in node.child_nodes():
                    __compute_dLeaf__(child, cumm_l + child.edge_length)

        children = self.ddpTree.root.child_nodes()
        crowded_child = None
        maxleaf = -1

        for node in children:
            if node.nleaf > maxleaf:
                maxleaf = node.nleaf
                crowded_child = node

        __compute_dLeaf__(children[1], 0)

        return D

    def filter_branch(self, threshold=None):
        # filter out abnormally long branches
        i = 1
        self.logger.info("Iteration: " + str(i))
        self.Reroot()
        while 1:
            check = self.filter_by_threshold(threshold=threshold)
            if (not check):
                self.logger.info("I could not remove anything more! I stop here!")
                break
            i += 1
            self.logger.info("Iteration: " + str(i))
            self.reset()
            self.Reroot()

    def filter_by_threshold(self, threshold=None, k=3.5):
        if threshold is None:
            threshold = self.compute_threshold(k=k)

        def __filter__(node, cumm_l):
            removed = False
            node.child_removed = False
            for child in node.child_nodes():
                check = __filter__(child, cumm_l + child.edge_length)
                removed = removed or check

            p = node.parent_node
            # if ( cumm_l > threshold ) or ( node.child_removed and len(node.child_nodes()) == 0 ):
            if (cumm_l > threshold) or (node.child_removed and node.num_children() == 0):
                # remove node
                p.remove_child(node)
                # update parent node
                p.child_removed = True
                removed = True
                try:
                    self.logger.info(node.label + " removed")
                except:
                    self.logger.info(node.name + " removed")
            # elif len(node.child_nodes()) == 1:
            elif node.num_child_nodes() == 1:
                # remove node and attach its only child to its parent
                e1 = node.edge_length
                child = node.child_nodes()[0]
                e2 = child.edge_length
                p.remove_child(node)
                node.remove_child(child)
                p.add_child(child)
                child.edge_length = e1 + e2
            return removed

        return __filter__(self.get_root(), 0)

    def compute_threhold(self, k=3.5):
        self.logger.warning("Abstract class! Should never be called")
        return 0

    def reset(self):
        self.logger.warning("Abstract class! Should never be called")

    def find_root(self):
        self.Topdown_label()  # temporarily included for debugging
        self.Bottomup_update()
        self.prepare_root()
        self.Topdown_update()

    def opt_score(self):
        self.logger.warning("Abstract class! Should never be called")

    def report_score(self):
        self.logger.warning("Abstract class! Should never be called")

    def Reroot(self):
        self.find_root()
        # self.report_score()
        # d2currRoot = 0
        # br2currRoot = 0
        if self.opt_root != self.ddpTree.root:
            # d2currRoot,br2currRoot = self.reroot_at_edge(self.opt_root.edge, self.opt_root.edge_length-self.opt_x, self.opt_x)
            self.reroot_at_edge(self.opt_root, self.opt_x)
            # self.ddpTree.reroot(self.opt_root,self.opt_x)

        # return head_id, tail_id, edge_length, self.opt_x
        # return d2currRoot,br2currRoot

    def Opt_function(self, node):
        self.logger.warning("Abstract method! Should never be called")

    def tree_as_newick(self, outstream=sys.stdout, label_by_name=False):
        # dendropy's method to write newick seems to have problem ...
        self.__write_newick(self.ddpTree.root, outstream, label_by_name=label_by_name)
        outstream.write(";\n")

    #            outstream.write(bytes(";\n", "ascii"))

    def __write_newick(self, node, outstream, label_by_name=False):
        if node.is_leaf():
            if label_by_name:
                outstream.write(str(node.name))
            #                    outstream.write(bytes(str(node.name), "ascii"))
            else:
                try:
                    outstream.write(node.label)
                #                        outstream.write(bytes(node.label, "ascii"))
                except:
                    outstream.write(node.label)
        #                        outstream.write(bytes(str(node.label), "ascii"))
        else:
            outstream.write('(')
            # outstream.write(bytes('(', "ascii"))
            is_first_child = True
            for child in node.child_nodes():
                if is_first_child:
                    is_first_child = False
                else:
                    outstream.write(',')
                #                        outstream.write(bytes(',', "ascii"))
                self.__write_newick(child, outstream, label_by_name=label_by_name)
            outstream.write(')')
        #                outstream.write(bytes(')', "ascii"))
        if not node.is_leaf():
            if label_by_name:
                outstream.write(str(node.name))
            #                    outstream.write(bytes(str(node.name), "ascii"))
            elif node.label is not None:
                outstream.write(str(node.label))
        #                    outstream.write(bytes(str(node.label), "ascii"))

        if not node.edge_length is None:
            outstream.write(":" + str(node.edge_length))

    #                outstream.write(bytes(":" + str(node.edge_length), "ascii"))

    def reroot_at_edge(self, node, length):
        # the method provided by dendropy DOESN'T seem to work ...
        # change edge to opt_root
        length1 = node.edge_length - length
        length2 = length
        if not node:
            return
        head = node  # opt_root = v = node
        tail = node.parent  # u parent of opt_root
        if not tail:
            return

        if (length2 == 0) and head.is_leaf():
            return 0, 0

        # new_root = self.ddpTree.node_factory()
        new_root = Node()

        tail.remove_child(head)

        new_root.add_child(head)
        head.edge_length = length2

        p = tail.parent
        l = tail.edge_length

        new_root.add_child(tail)
        tail.edge_length = length1

        br2currRoot = 0
        d2currRoot = length1

        #            if tail.label == self.ddpTree.root.label:
        if (tail is self.ddpTree.root):
            head = new_root

        while tail is not self.ddpTree.root:
            # MAD@ add
            # q = tail.parent #tail should have 2 parents right now: new_root and its old parent
            q = head.parent
            # End MAD@ add
            head = tail
            tail = p
            p = tail.parent

            br2currRoot += 1
            d2currRoot += l

            l1 = tail.edge_length
            tail.remove_child(head)
            # MAD@ add
            head.parent = q
            # End MAD@ add

            head.add_child(tail)
            tail.edge_length = l
            l = l1

        # out of while loop: tail IS now tree.root
        if tail.num_children() == 1:
            # merge the 2 branches of the old root and adjust the branch length
            # sis = [child for child in tail.child_nodes()][0]
            sis = tail.child_nodes()[0]
            l = sis.edge_length
            tail.remove_child(sis)
            head.add_child(sis)
            sis.edge_length = l + tail.edge_length
            head.remove_child(tail)
            # tail.remove_child(head)

        new_root.name = self.ddpTree.root.name
        self.ddpTree.root.name = "OLD"
        self.ddpTree.root = new_root

        ### MAD@ add
        #            for node in self.ddpTree.traverse_postorder():
        #                for child in node.child_nodes():
        #                    if child.parent_node is not node:
        #                        logger.info("Error found!")
        #                        child.parent_node = node
        ### MAD@ add

        return d2currRoot, br2currRoot

    def get_root(self):
        return self.ddpTree.root


class OGR_Tree(Tree_extend):
    # supportive class to implement outgroup-reroot (OGR = outgroup reroot, hence the name)
    # this rooting method solve the difficulty in finding the root when there are mulitple outgroups
    # and they are not monophyletic. It seeks for the rooting place that maximizes the triplet score
    # of the specified outgroups.
    def __init__(self, outgroups, ddpTree=None, tree_file=None, schema="newick", logger_id=1, logger_stream=sys.stderr):
        super(OGR_Tree, self).__init__(ddpTree, tree_file, schema)
        self.logger = new_logger("OGR_Tree_" + str(logger_id), myStream=logger_stream)
        # L = self.ddpTree.leaf_nodes()
        L = []
        for leaf in self.ddpTree.traverse_leaves():
            L.append(leaf)
        self.OGs = set([x.label for x in L if x.label in set(outgroups)])
        self.nOGs = len(self.OGs)
        self.nIGs = len(L) - self.nOGs
        self.max_nTrpls = self.nIGs * self.nOGs * (self.nOGs - 1) / 2 + self.nOGs * self.nIGs * (self.nIGs - 1) / 2
        self.reset()

    def reset(self):
        self.opt_root = self.ddpTree.root
        self.opt_nTrpls = 0

    def Node_init(self, node, nTrpl_in=0, nTrpl_out=0, nOGs=0, nIGs=0):
        node.nTrpl_in = nTrpl_in
        node.nTrpl_out = nTrpl_out
        node.nOGs = nOGs
        node.nIGs = nIGs

    def Opt_function(self, node):
        curr_nTrpls = node.nTrpl_in + node.nTrpl_out
        if curr_nTrpls > self.opt_nTrpls:
            self.opt_nTrpls = curr_nTrpls
            self.opt_root = node
            self.opt_x = node.edge_length / 2  # NOTE: this method does not consider branch length, the *middle point* of the edge is just arbitrarily chosen

    def bUp_update(self, node):
        if node.is_leaf():
            node.nOGs = 1 if node.label in self.OGs else 0
            node.nIGs = 1 if node.nOGs == 0 else 0
        else:
            C = node.child_nodes()

            node.nOGs = sum([c.nOGs for c in C])
            node.nIGs = sum([c.nIGs for c in C])

            node.nTrpl_in = sum([c.nTrpl_in for c in C])

            for i, c1 in enumerate(C):
                for c2 in C[i + 1:]:
                    IG_trpls = c1.nIGs * c2.nIGs * (self.nOGs - node.nOGs)
                    OG_trpls = c1.nOGs * c2.nOGs * (self.nIGs - node.nIGs)
                    node.nTrpl_in += IG_trpls + OG_trpls

    def tDown_update(self, node, opt_function):
        C = node.child_nodes()

        for child in C:
            C1 = [c for c in C if c is not child]
            child.nTrpl_out = node.nTrpl_out

            for i, c1 in enumerate(C1):
                child.nTrpl_out += c1.nTrpl_in
                child.nTrpl_out += (self.nIGs - node.nIGs) * c1.nIGs * child.nOGs
                child.nTrpl_out += (self.nOGs - node.nOGs) * c1.nOGs * child.nIGs

                for c2 in C1[i + 1:]:
                    IG_trpls = c1.nIGs * c2.nIGs * child.nOGs
                    OG_trpls = c1.nOGs * c2.nOGs * child.nIGs

                    child.nTrpl_out += IG_trpls + OG_trpls

            opt_function(child)

    def prepare_root(self):
        pass

    def opt_score(self):
        return self.opt_nTrpls / float(self.max_nTrpls) if self.max_nTrpls != 0 else None

    def report_score(self):
        myScore = self.opt_score()
        if myScore is None:
            self.logger.warning("OG rooting failed because the tree has no outgroup")
        return "Triplet score: " + str(self.opt_score())


class MPR_Tree(Tree_extend):
    # supportive class to implement midpoint-reroot (mpr = mid point reroot, hence the name)G
    def __init__(self, ddpTree=None, tree_file=None, schema="newick", logger_id=1, logger_stream=sys.stderr):
        super(MPR_Tree, self).__init__(ddpTree, tree_file, schema)
        self.logger = new_logger("MPR_Tree_" + str(logger_id), myStream=logger_stream)
        self.reset()

    def reset(self):
        self.max_distance = -1
        self.opt_root = self.ddpTree.root
        self.opt_x = 0

    def Node_init(self, node, max_in=None, max_out=-1):
        node.max_in = max_in if max_in else [0, 0]
        node.max_out = max_out

    def Opt_function(self, node):
        m = max(node.max_in)
        curr_max_distance = m + node.max_out
        x = (node.max_out - m) / 2
        if curr_max_distance > self.max_distance and x >= 0 and x <= node.edge_length:
            self.max_distance = curr_max_distance
            self.opt_x = x
            self.opt_root = node

    def bUp_update(self, node):
        if not node.is_leaf():
            node.max_in = []
            for child in node.child_nodes():
                node.max_in.append(max(child.max_in) + child.edge_length)

    def tDown_update(self, node, opt_function):
        child_idx = 0
        for child in node.child_nodes():
            child.max_out = max([node.max_out] + [node.max_in[k] for k in range(len(node.max_in))
                                                  if k != child_idx]) + child.edge_length
            opt_function(child)
            child_idx += 1

    def prepare_root(self):
        pass

    def compute_threhold(self, k=3.5):
        self.logger.warning("Trying to compute threshold for MPR_Tree, which is not supported.")
        return 0

    def opt_score(self):
        return self.max_distance / 2

    def report_score(self):
        return "Tree height: " + str(self.opt_score())


########################################################################################################################

'''
logger = logging.getLogger("MinVar")
logger.setLevel(logging.INFO)
handler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter('%(levelname)s:%(name)s:%(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.propagate = False
'''

class minVAR_Base_Tree(Tree_extend):
    # supportive base class to implement VAR-reroot, hence the name
    def __init__(self, ddpTree=None, tree_file=None, schema="newick",logger_id=1,logger_stream=sys.stderr):
        super(minVAR_Base_Tree, self).__init__(ddpTree, tree_file, schema)
        self.logger = new_logger("MinVar_Tree_" + str(logger_id),myStream=logger_stream)
        self.reset()

    def reset(self):
        self.minVAR = None
        self.opt_root = self.ddpTree.root
        self.opt_x = 0

    def Node_init(self, node, nleaf=1, sum_in=0, sum_total=0, var=-1):
        node.sum_in = sum_in
        node.sum_total = sum_total
        node.nleaf = nleaf
        node.var = var

    def Opt_function(self, node, a, b, c):
        self.logger.info("Abstract method! Should never be called")

    def compute_dRoot_VAR(self):
        cumm = {'ssq': 0, 'sum': 0}

        def compute_dRoot(node, cumm_l):
            if node.is_leaf():
                cumm['ssq'] += cumm_l ** 2
                cumm['sum'] += cumm_l
            else:
                for child in node.child_nodes():
                    compute_dRoot(child, cumm_l + child.edge_length)

        compute_dRoot(self.get_root(), 0)
        N = self.get_root().nleaf
        root_var = cumm['ssq'] / N - (cumm['sum'] / N) ** 2
        self.get_root().var = root_var

    def bUp_update(self, node):
        if node.is_leaf():
            node.nleaf = 1
            node.sum_in = 0
        else:
            node.nleaf = 0
            node.sum_in = 0
            for child in node.child_nodes():
                node.nleaf += child.nleaf
                node.sum_in += child.sum_in + child.nleaf * child.edge_length

    def Update_var(self, child, node, edge_length):
        alpha = 2 * (node.sum_total - 2 * (child.sum_in + child.nleaf * edge_length)) / self.total_leaves
        beta = 1 - 2 * float(child.nleaf) / self.total_leaves
        a = 1 - beta * beta
        b = alpha - 2 * node.sum_total * beta / self.total_leaves
        c = node.var
        child.var = a * edge_length * edge_length + b * edge_length + c
        return a, b, c

    def tDown_update(self, node, opt_function):
        for child in node.child_nodes():
            child.sum_total = node.sum_total + (self.total_leaves - 2 * child.nleaf) * child.edge_length
            a, b, c = self.Update_var(child, node, child.edge_length)
            opt_function(child, a, b, c)

    def prepare_root(self):
        root = self.get_root()
        root.sum_total = root.sum_in
        self.compute_dRoot_VAR()
        self.total_leaves = root.nleaf

    def opt_score(self):
        return self.minVAR

    def report_score(self):
        return "MinVar score: " + str(self.opt_score())


class MVDF_Tree(minVAR_Base_Tree):
    # supportive class to implement VAR-reroot + deepest node + factorization
    def __init__(self, ddpTree=None, tree_file=None, schema="newick"):
        super(MVDF_Tree, self).__init__(ddpTree, tree_file, schema)
        self.deep_node = None

    def reset(self):
        super(MVDF_Tree, self).reset()
        self.deep_node = None

    def Opt_function(self, node, a, b, c):
        x = -b / (2 * a)
        if x >= 0 and x <= node.edge_length:
            #                curr_minVAR = a*x*x + b*x + c
            factor = float(node.nleaf) / self.total_leaves
            factor = factor * (1 - factor)
            curr_minVAR = (a * x * x + b * x + c) / factor

            if node.var < node.parent_node.var:
                deep_node = node
            else:
                deep_node = node.parent_node

            updateNeed = False
            if (self.deep_node is None) or (deep_node.var < self.deep_node.var):
                self.deep_node = deep_node
                self.minVAR = curr_minVAR
                updateNeed = True
            elif (self.deep_node is deep_node) and (curr_minVAR < self.minVAR):
                self.minVAR = curr_minVAR
                updateNeed = True

            if updateNeed:
                self.opt_root = node
                self.opt_x = node.edge_length - x

    #                self.logger.info(str(curr_minVAR) + "\t" + node.label
    #                      + "\t" + str(node.edge_length-x) + "\t" + str(self.Tree_records[node.idx].var)
    #                      + "\t" + (str(node.parent_node.label) if node.parent_node else "None")
    #                      + "\t" + str(self.Tree_records[node.parent_node.idx].var))

    def compute_threshold(self, k=3.5):
        # should be called only AFTER the MV root was found
        mean = (self.opt_root.sum_total - self.opt_x *
                (self.total_leaves - 2 * self.opt_root.nleaf)) / self.total_leaves
        factor = float(self.opt_root.nleaf) / self.total_leaves
        factor = factor * (1 - factor)
        rootVar = self.minVAR * factor
        self.logger.info(mean)
        self.logger.info(rootVar)
        std = math.sqrt(rootVar)
        return mean + k * std


class MVD0_Tree(minVAR_Base_Tree):
    # supportive class to implement VAR-reroot + deepest node + no factorization
    def __init__(self, ddpTree=None, tree_file=None, schema="newick"):
        super(MVD0_Tree, self).__init__(ddpTree, tree_file, schema)
        self.deep_node = None

    def reset(self):
        super(MVD0_Tree, self).reset()
        self.deep_node = None

    def Opt_function(self, node, a, b, c):
        x = -b / (2 * a)
        if x >= 0 and x <= node.edge_length:
            curr_minVAR = a * x * x + b * x + c

            if node.var < node.parent_node.var:
                deep_node = node
            else:
                deep_node = node.parent_node

            updateNeed = False
            if (self.deep_node is None) or (deep_node.var < self.deep_node.var):
                self.deep_node = deep_node
                self.minVAR = curr_minVAR
                updateNeed = True
            elif (self.deep_node is deep_node) and (curr_minVAR < self.minVAR):
                self.minVAR = curr_minVAR
                updateNeed = True

            if updateNeed:
                self.opt_root = node
                self.opt_x = node.edge_length - x

    #                self.logger.info(str(curr_minVAR) + "\t" + node.label
    #                      + "\t" + str(node.edge_length-x) + "\t" + str(self.Tree_records[node.idx].var)
    #                      + "\t" + (str(node.parent_node.label) if node.parent_node else "None")
    #                      + "\t" + str(self.Tree_records[node.parent_node.idx].var))

    def compute_threshold(self, k=3.5):
        # should be called only AFTER the MV root was found
        mean = (self.opt_root.sum_total - self.opt_x *
                (self.total_leaves - 2 * self.opt_root.nleaf)) / self.total_leaves
        self.logger.info(mean)
        self.logger.info(self.minVAR)
        std = math.sqrt(self.minVAR)
        return mean + k * std


class MV0F_Tree(minVAR_Base_Tree):
    # supportive class to implement VAR-reroot + no deepest node + factorization
    #        def __init__(self, ddpTree = None, tree_file = None, schema = "newick"):
    #            super().__init__(ddpTree, tree_file, schema)

    def Opt_function(self, node, a, b, c):
        x = -b / (2 * a)
        if x >= 0 and x <= node.edge_length:
            #                curr_minVAR = a*x*x + b*x + c
            factor = float(node.nleaf) / self.total_leaves
            factor = factor * (1 - factor)
            curr_minVAR = (a * x * x + b * x + c) / factor
            if self.minVAR is None or curr_minVAR < self.minVAR:
                self.minVAR = curr_minVAR
                self.opt_root = node
                self.opt_x = node.edge_length - x

    #                self.logger.info(str(curr_minVAR) + "\t" + node.label
    #                      + "\t" + str(node.edge_length-x) + "\t" + str(self.Tree_records[node.idx].var)
    #                      + "\t" + (str(node.parent_node.label) if node.parent_node else "None")
    #                      + "\t" + str(self.Tree_records[node.parent_node.idx].var))

    def compute_threshold(self, k=3.5):
        # should be called only AFTER the MV root was found
        mean = (self.opt_root.sum_total - self.opt_x *
                (self.total_leaves - 2 * self.opt_root.nleaf)) / self.total_leaves
        factor = float(self.opt_root.nleaf) / self.total_leaves
        factor = factor * (1 - factor)
        rootVar = self.minVAR * factor
        self.logger.info(mean)
        self.logger.info(rootVar)
        std = math.sqrt(rootVar)
        return mean + k * std


class MV00_Tree(minVAR_Base_Tree):
    # supportive class to implement VAR-reroot + no deepest node + no factorization
    #        def __init__(self, ddpTree = None, tree_file = None, schema = "newick"):
    #            super().__init__(ddpTree, tree_file, schema)

    def Opt_function(self, node, a, b, c):
        x = -b / (2 * a)
        if x >= 0 and x <= node.edge_length:
            curr_minVAR = a * x * x + b * x + c
            if self.minVAR is None or curr_minVAR < self.minVAR:
                self.minVAR = curr_minVAR
                self.opt_root = node
                self.opt_x = node.edge_length - x

    def compute_threshold(self, k=3.5):
        # should be called only AFTER the MV root was found
        mean = (self.opt_root.sum_total - self.opt_x *
                (self.total_leaves - 2 * self.opt_root.nleaf)) / self.total_leaves
        self.logger.info(mean)
        self.logger.info(self.minVAR)
        std = math.sqrt(self.minVAR)
        return mean + k * std


class MBR_Tree(Tree_extend):
    # supportive class to implement midpoint balance root
    def __init__(self, ddpTree=None, tree_file=None, schema="newick"):
        super(MBR_Tree, self).__init__(ddpTree, tree_file, schema)

        self.BPs = []  # BPs : balance points
        self.opt_root = self.ddpTree.root
        self.opt_x = 0

    def Node_init(self, node, nleaf=1, sum_in=0, sum_out=-1):
        self.nleaf = nleaf
        self.sum_in = sum_in
        self.sum_out = sum_out

    def Opt_function(self, node):
        nleaf = node.nleaf
        mean_in = node.sum_in / nleaf
        mean_out = node.sum_out / (self.total_leaves - nleaf)
        x = (mean_out - mean_in) / 2
        if x >= 0 and x <= node.edge_length:
            self.BPs.append((node, x, mean_in + x))
            node.x = x
            node.mean = mean_in + x
        else:
            node.x = None
            node.mean = None

    def bUp_update(self, node):
        node.sum_in = 0
        if node.is_leaf():
            node.nleaf = 1
        else:
            node.nleaf = 0
            for child in node.child_nodes():
                node.nleaf += child.nleaf
                node.sum_in += child.sum_in + child.nleaf * child.edge_length

    def tDown_update(self, node, opt_function):
        child_idx = 0
        for child in node.child_nodes():
            child.sum_out = (node.sum_out + node.sum_in + child.edge_length *
                             (self.total_leaves - 2 * child.nleaf) - child.sum_in)
            opt_function(child)
            child_idx += 1

    def prepare_root(self):
        root = self.get_root()
        root.sum_out = 0
        self.total_leaves = root.nleaf
        root.x = None
        root.mean = None

    def list_balance_points(self):
        self.Topdown_label()
        self.Bottomup_update()
        self.prepare_root()
        self.Topdown_update()

        for (node, x, mean) in self.BPs:
            if node.is_leaf():
                #                  self.logger.info(node.label + "\t" + str(x) + "\t" + str(mean))
                self.logger.info(node.label + "\t" + str(x) + "\t" + str(mean))
            else:
                self.logger.info(node.label + "\t" + str(x) + "\t" + str(mean))

    def build_balance_tree(self):
        self.Topdown_label()  # keep this step for now for debugging purpose
        self.Bottomup_update()
        self.prepare_root()
        self.Topdown_update()

        # self.list_balance_points()

        self.balance_tree = self.ddpTree.extract_tree()

        # bottom up pruning
        for node in self.balance_tree.traverse_postorder():
            node.type = "real"
            node.BPbelow = False

            '''if node.is_leaf():
                self.logger.info("parent: " + node.label)# + "\t" + str(node.extraction_source.x))
            else:
                self.logger.info("parent: " + node.label)#+ "\t" + str(node.extraction_source.x))'''

            for ch in node.child_nodes():
                '''try:
                    self.logger.info("child: " + ch.label)# + "\t" + str(ch.extraction_source.x))
                except:
                    self.logger.info("child: " + ch.label) #+ "\t" + str(ch.extraction_source.x))'''

                if ch.BPbelow or (ch.extraction_source.x is not None):
                    node.BPbelow = True
                # node.BPbelow = node.BPbelow or ch.BPbelow or (ch.extraction_source.x is not None)

                if not ch.BPbelow:
                    # remove the whole clade under ch
                    # for ch1 in ch.child_nodes():
                    #    ch.remove_child(ch1)
                    edgelen = ch.edge_length
                    node.remove_child(ch)

                    if ch.extraction_source.x is not None:
                        # add a new node p at the balance point
                        # set p to be a child of node (edge length ch.edge_length - x)
                        # add a new node ch1 to be another child of p (edge length ch.mean)
                        edgelen = ch.edge_length

                        # p = self.ddpTree.node_factory()
                        # ch1 = self.ddpTree.node_factory()
                        p = Node()
                        ch1 = Node()

                        p.type = "bp"  # bp: balance-point
                        p.ref_child = ch.extraction_source  # link p to the original tree (for later use after finding midpoint)
                        ch1.type = "dm"  # dm: dummy

                        # node.remove_child(ch)
                        node.add_child(p)
                        p.add_child(ch1)

                        p.edge_length = edgelen - ch.extraction_source.x
                        ch1.edge_length = ch.extraction_source.mean

                elif ch.extraction_source.x is not None:
                    # add a new node p at the balance point
                    # set p to be a child of node (edge length ch.edge_length - x)
                    # set ch to be a child of p (edge length x)
                    # add a new node ch1 to be another child of p (edge length ch.mean)

                    edgelen = ch.edge_length

                    # p = self.ddpTree.node_factory()
                    p = Node()
                    # ch1 = self.ddpTree.node_factory()
                    ch1 = Node()

                    p.type = "bp"
                    p.ref_child = ch.extraction_source  # link p to the original tree (for later use after finding midpoint)
                    ch1.type = "dm"

                    node.remove_child(ch)
                    node.add_child(p)
                    p.add_child(ch)
                    p.add_child(ch1)

                    ch.edge_length = ch.extraction_source.x
                    p.edge_length = edgelen - ch.extraction_source.x
                    ch1.edge_length = ch.extraction_source.mean

                    # topdown pruning
        node = self.balance_tree.root
        nchild = len(node.child_nodes())
        while nchild > 0 and nchild < 2:
            # node has less than 2 children
            temp = node
            node = node.child_nodes()[0]
            temp.remove_child(node)
            if node.type == "dm":
                node = temp
                break
            nchild = len(node.child_nodes())

        self.balance_tree.root = node
        self.balance_tree.root.edge_length = None
        # balance_tree.root = None

        # mptre = MPR_Tree(ddpTree=balance_tree)
        # mptre.tree_as_newick()

        # return balance_tree

    def find_root(self):
        self.build_balance_tree()
        mptre = MPR_Tree(ddpTree=self.balance_tree)
        mptre.tree_as_newick()
        mptre.find_root()

        self.logger.info(mptre.opt_root.type)

        if mptre.opt_root.type == "bp":
            self.opt_root = mptre.opt_root.ref_child
            self.opt_x = mptre.opt_root.ref_child.x + mptre.opt_x
        elif mptre.opt_root.type == "dm":
            self.logger.info("Hmm... Is it possible that a dummy was found as the opt_root?")
        else:
            self.opt_root = mptre.opt_root.extraction_source
            self.opt_x = mptre.opt_x

        self.logger.info(self.opt_root.label)
        self.logger.info(self.opt_x)


########################################################################################################################

logger = logging.getLogger("quadprog_solvers")
logger.setLevel(logging.INFO)
handler = logging.StreamHandler(stdout)
formatter = logging.Formatter('%(levelname)s:%(name)s:%(message)s')
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.propagate = False

def cvxopt_solve_qp(P, q, G=None, h=None, A=None, b=None, maxIter=1000):
    P = .5 * (P + P.T)  # make sure P is symmetric
    args = [cvxopt.matrix(P), cvxopt.matrix(q)]
    if G is not None:
        args.extend([cvxopt.matrix(G), cvxopt.matrix(h)])
        if A is not None:
            args.extend([cvxopt.matrix(A), cvxopt.matrix(b)])
    sol = cvxopt.solvers.qp(*args,options={'show_progress':False,'maxiters':maxIter})
    if 'optimal' not in sol['status']:
        if "unknown" in sol['status']:
            logger.warning("Couldn't find optimal solution on one branch. Perhaps due to maximum iterations exceeded. Consider increasing the maximum iterations via -x.")
        else:
            logger.warning("Couldn't find optimal solution on one branch. Solution status: " + sol['status'])
    #return None
    return numpy.array(sol['x']).reshape((P.shape[1],))


########################################################################################################################

EPSILON = 1e-5

class RTT_Tree(Tree_extend):
    # supportive base class to implement RTT-reroot, hence the name
    def __init__(self, smplTimes, ddpTree=None, tree_file=None, schema="newick", logger_id=1, logger_stream=stderr,
                 maxIter=1000):
        super(RTT_Tree, self).__init__(ddpTree, tree_file, schema)
        self.logger = new_logger("RTT_Tree_" + str(logger_id), myStream=logger_stream)
        self.smplTimes = smplTimes
        self.reset()
        self.maxIter = maxIter

    def reset(self):
        self.RTT = None
        self.opt_root = self.ddpTree.root
        self.opt_y = 0
        self.opt_x = 0
        self.opt_mu = 0
        self.tmin = min(self.smplTimes.values())

    def Node_init(self, node, nleaf=1, SDI=0, SD=0, ST=0, SDT=0, SSD=0):
        node.SDI = SDI
        node.SD = SD
        node.nleaf = nleaf
        node.ST = ST
        node.SDT = SDT
        node.SSD = SSD

    def Opt_function(self, node, SST, deltaT, deltaD, SDT, SSD, ST, SD):
        n = self.total_leaves
        a, b, c, d, e, f = n, SST, (-2 * deltaT), (2 * deltaD), (-2 * SDT), SSD
        k, m, r = 2 * (n - 2 * node.nleaf), -2 * ST, 2 * SD

        tmin = self.tmin

        # use quadprog to compute mu_star, y_star, and x_star
        P = array([[a, k / 2, c / 2.], [k / 2, n, m / 2], [c / 2, m / 2, b]])
        q = array([d / 2., r / 2, e / 2])
        G = array([[-1., 0., 0.], [0., 0., -1.], [1., 0., 0.], [0., 1., -tmin]])
        h = array([0., EPSILON, node.edge_length, 0]).reshape((4,))
        solution = cvxopt_solve_qp(P, q, G, h, maxIter=self.maxIter)
        x_star = solution[0]
        y_star = solution[1]
        mu_star = solution[2]
        curr_RTT = a * x_star * x_star + b * mu_star * mu_star + c * x_star * mu_star + d * x_star + e * mu_star + f + n * y_star * y_star + k * x_star * y_star + m * mu_star * y_star + r * y_star

        if self.RTT is None or (curr_RTT - self.RTT < -EPSILON):
            self.RTT = curr_RTT
            self.opt_root = node
            self.opt_x = node.edge_length - x_star
            self.opt_y = y_star
            self.opt_mu = mu_star

    def bUp_update(self, node):
        if node.is_leaf():
            node.nleaf = 1
            node.SDI = 0
            node.ST = self.smplTimes[node.label]
        else:
            node.nleaf = 0
            node.SDI = 0
            node.ST = 0
            for child in node.child_nodes():
                node.nleaf += child.nleaf
                node.SDI += child.SDI + child.nleaf * child.edge_length
                node.ST += child.ST

    def Update_var(self, child, node, edge_length):
        SST = self.SST
        deltaT = self.ddpTree.root.ST - 2 * child.ST
        deltaD = -2 * child.nleaf * edge_length - 2 * child.SDI + node.SD
        SDT = node.SDT
        SSD = node.SSD
        ST = self.ST
        SD = node.SD
        return SST, deltaT, deltaD, SDT, SSD, ST, SD

    def tDown_update(self, node, opt_function):
        for child in node.child_nodes():
            child.SD = node.SD + (self.total_leaves - 2 * child.nleaf) * child.edge_length
            child.SDT = node.SDT + child.edge_length * (self.ddpTree.root.ST - 2 * child.ST)
            child.SSD = node.SSD + (self.total_leaves - 4 * child.nleaf) * (child.edge_length ** 2) + 2 * (
                        node.SD - 2 * child.SDI) * child.edge_length
            SST, deltaT, deltaD, SDT, SSD, ST, SD = self.Update_var(child, node, child.edge_length)
            opt_function(child, SST, deltaT, deltaD, SDT, SSD, ST, SD)

    def prepare_root(self):
        root = self.get_root()
        root.SD = root.SDI
        self.total_leaves = root.nleaf
        self.ST = root.ST
        self.ddpTree.root.droot = 0
        self.ddpTree.root.troot = 0
        root.SD, root.SSD, root.SDT, self.SST = 0, 0, 0, 0
        for v in self.ddpTree.traverse_preorder():
            if not v.is_root():
                # must have defined edge lengths
                v.droot = v.parent.droot + v.edge_length
                if v.is_leaf():
                    root.SSD += (v.droot ** 2)
                    self.SST += (self.smplTimes[v.label] ** 2)
                    root.SD += v.droot
                    root.SDT += (v.droot * self.smplTimes[v.label])

    def opt_score(self):
        return self.RTT

    def return_values(self):
        # returns RTT score, mu, t0
        return (self.opt_score() / self.total_leaves), (self.opt_mu), (self.opt_y / self.opt_mu)

    def report_score(self):
        return "RTT=" + str(self.opt_score() / self.total_leaves) + "\tmu=" + str(self.opt_mu) + "\tt0=" + str(
            self.opt_y / self.opt_mu)


########################################################################################################################


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input',        required=False, type=argparse.FileType('r'), default=stdin,     help="Input File (default is STDIN)")
    parser.add_argument('-m', '--method',       required=False, type=str, default="MV",                         help="Method (MP for midpoint, MV for minVAR, OG for outgroup, RTT for root-to-tip) (default is MV)")
    parser.add_argument('-g', '--outgroups',    required=False, type=str,                                       help="Specify the outgroups. If specifying a list of outgroups, put them between quotes (i.e. \"). Otherwise, specifying a file which containts all the outgroups. Can only be used with -m OG")
    parser.add_argument('-t', '--smplTimes',    required=False, type=argparse.FileType('r'),                    help="The file containing the sampling times at leaves; to be used with -m RTT")
    parser.add_argument('-o', '--outfile',      required=False, type=argparse.FileType('w'), default=stdout,    help="Output File (default is STDOUT)")
    parser.add_argument('-s', '--schema',       required=False, type=str, default="newick",                     help="Schema of your input treefile (default is newick)")
    parser.add_argument('-f', '--infofile',     required=False, type=argparse.FileType('w'), default=None,      help="Save all the logging to this file. Default: print to stderr")
    parser.add_argument("-v", "--version",      action='version', version=PROGRAM_NAME + " " + PROGRAM_VERSION, help="Show FastRoot version and exit")
    parser.add_argument("-x", "--maxIter",      required=False, type=int, default=None,                         help="Maximum number of iterations to run cvxopt")

    # print help message if no argument is given
    if len(argv) == 1:
        logger = new_logger(__name__)
        logger.info("Running " + PROGRAM_NAME + " version " + PROGRAM_VERSION)
        parser.print_help()
        exit(0)

    args = parser.parse_args()
    stream = args.infofile if args.infofile else stderr
    logger = new_logger(__name__, myStream=stream)
    logger.info("Running " + PROGRAM_NAME + " version " + PROGRAM_VERSION)
    new_argv = re.sub(' +', ' ', " ".join(argv).replace("\n", " "))
    logger.info(PROGRAM_NAME + " was called as follows: " + new_argv)

    METHOD2FUNC = {'MP': MPR_Tree, 'MV': MV00_Tree, 'OG': OGR_Tree, 'RTT': RTT_Tree}
    METHOD2DESC = {'MP': "Midpoint", 'MV': "MinVar", 'OG': "Outgroup", 'RTT': "Root-to-Tip"}
    method = args.method.upper()

    # reading outgroups
    if args.outgroups:
        if method != 'OG':
            method = 'OG'
            logger.warning("The rooting method is set to outgroup rooting (OG) due to the presence of outgroups")
        if os.path.exists(args.outgroups):
            OGs = []
            for line in open(args.outgroups, 'r'):
                OGs.append(line.strip())
        else:
            OGs = args.outgroups.split()
    else:
        OGs = None

    # reading sampling times
    if args.smplTimes:
        smplTimes = {}
        for line in args.smplTimes:
            sp, t = line.strip().split()
            smplTimes[sp] = float(t)
        if method != 'RTT':
            method = 'RTT'
            logger.warning(
                "The rooting method is set to root-to-tip rooting (RTT) due to the presence of sampling times")

    if method == 'RTT' and args.smplTimes is None:
        logger.error("Need sampling times for root-to-tip rooting")
        exit()
    elif method == 'OG' and args.outgroups is None:
        logger.error("Need outgroups for outgroup rooting")
        exit()

    assert method in METHOD2FUNC, "Invalid method! Valid options: MP for midpoint, MV for minVAR, OG for outgroups, RTT for root-to-tip"
    logger.info("Rooting Method: " + METHOD2DESC[method] + " Rooting")

    maxIter = 0
    if method == 'RTT':
        if args.maxIter and args.maxIter < 1000:
            logger.warning(
                "Invalid number of maximum iterations (-x). Must be at least 1000. Set back to 1000 by default.")
        maxIter = max(1000, args.maxIter) if args.maxIter else 1000
        logger.info("Maximum iterations: " + str(maxIter))
    elif args.maxIter is not None:
        logger.warning("The maximum number of iterations (-x) is only used with root-to-tip rooting (RTT)")

    # read and root each tree
    for i, line in enumerate(args.input):
        tree = read_tree(line, schema=args.schema.lower())
        if method == 'OG':
            a_tree = OGR_Tree(OGs, ddpTree=tree, logger_id=i + 1, logger_stream=stream)
        elif method == 'RTT':
            a_tree = RTT_Tree(smplTimes, ddpTree=tree, logger_id=i + 1, logger_stream=stream, maxIter=maxIter)
        else:
            a_tree = METHOD2FUNC[method](ddpTree=tree, logger_id=i + 1, logger_stream=stream)

        a_tree.Reroot()
        logger.info("Tree " + str(i + 1) + " " + a_tree.report_score())
        a_tree.tree_as_newick(outstream=args.outfile)


if __name__ == "__main__":
    main()


'''


'''