""" Here, we will use UPGMA and needleman_wunsh_multiple to align multiple sequences"""

import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from modules.nw_multiple import needleman_wunsch_multiple as nwm
from modules.UPGMA import UPGMA

print_results=False


def rec_align_tree(tree_of_indexes,list_of_seq, blosum_m, gap_opening_score, gap_extension_score,print_result, identity_score, substitution_score):
    """Recursively aligns a set of sequences based on a hierarchical tree structure.

    Parameters:
    - tree_of_indexes: A tree structure representing the hierarchical clustering of sequences.
    - list_of_seq: A list of sequences to be aligned.
    - blosum_m: A substitution matrix (e.g., BLOSUM) used for scoring alignments.
    - gap_opening_score: The penalty for opening a gap in the alignment.
    - gap_extension_score: The penalty for extending an existing gap.
    - print_result: Boolean indicating whether to print alignment results.
    - identity_score: Score assigned for identical matches.
    - substitution_score: Score assigned for substitutions.

    Returns:
    - A single aligned sequence or alignment block.
    """
    if tree_of_indexes.val is None:
        return []
    if len(tree_of_indexes.val)==1:
        return [list_of_seq[tree_of_indexes.val[0]]]
    else:
        block1=rec_align_tree(tree_of_indexes.left, list_of_seq, blosum_m, gap_opening_score, gap_extension_score,print_result, identity_score, substitution_score)
        block2=rec_align_tree(tree_of_indexes.right, list_of_seq, blosum_m, gap_opening_score, gap_extension_score,print_result, identity_score, substitution_score)
        print( "aligning : ", block1, block2)
        return nwm(block1, block2, blosum_m, gap_opening_score, gap_extension_score,print_result, identity_score, substitution_score)[1]

def align_multiple_sequences(list_of_seq, blosum_m, gap_opening_score=-10, gap_extension_score=-2, identity_score=1, substitution_score=-1):
    """
    Aligns multiple sequences using progressive alignment guided by a UPGMA tree.

    Parameters:
    - list_of_seq: A list of sequences to be aligned.
    - blosum_m: A substitution matrix (e.g., BLOSUM) used for scoring alignments.
    - gap_opening_score: The penalty for opening a gap in the alignment.
    - gap_extension_score: The penalty for extending an existing gap.
    - identity_score: Score assigned for identical matches (default: 1).
    - substitution_score: Score assigned for substitutions (default: -1).

    Returns:
    - The final multiple sequence alignment.
    """
    upgma_tree= UPGMA(list_of_seq, blosum_m, gap_opening_score, gap_extension_score, identity_score, substitution_score)
    alignement = rec_align_tree(upgma_tree, list_of_seq,blosum_m, gap_opening_score, gap_extension_score, print_results, identity_score, substitution_score)

    print("Multiple Alignement :", alignement)
    return alignement

    

