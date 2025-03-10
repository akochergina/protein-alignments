import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from modules.tree import Tree
from modules.UPGMA import UPGMA

def test_create_and_print_tree():
    tree_example_l = Tree(val=[1])
    tree_example_r = Tree(val=[2])
    tree_example= Tree([0,1], tree_example_l, tree_example_r)
    tree_example.print_tree()


def test_ugpma():
    list_of_seq=["CHAT","CAT","HER"]
    identity_score=1
    substitution_score=-1
    gap_opening_score=-10
    gap_extension_score = -2
    test_tree=UPGMA(list_of_seq, False, gap_opening_score, gap_extension_score, identity_score, substitution_score)
    #test_tree.print_tree()
