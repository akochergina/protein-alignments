'''
UGPMA Module

This module aims at constructing a hierarchical tree of clusters of sequences.
Its main function takes as a argument a list of sequences, and returns a Tree of the cluster.
The tree structure is defined in tree.py
There is an intermediate step of computing a distance matrix between each sequence

'''

import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from modules.needleman_wunsch import needleman_wunsch as nw
from modules.tree import Tree

def mat_distances(liste_de_sequences, blosum_m, gap_opening_score, gap_extension_score, identity_score=1, substitution_score=-1):
  """
    Computes a pairwise distance matrix for a list of sequences using the Needleman-Wunsch algorithm.

    Parameters:
    ----------
    sequence_list : list of str
        A list of sequences to compare.
    blosum_m : bool
        If True, uses the BLOSUM62 substitution matrix.
    gap_opening_score : int
        The penalty for opening a gap in the alignment.
    gap_extension_score : int
        The penalty for extending an existing gap.
    identity_score : int, optional
        The score for matching identical characters (default is 1).
    substitution_score : int, optional
        The score for aligning non-identical characters (default is -1).

    Returns:
    -------
    distance_matrix : list of list of int
        A symmetric matrix containing pairwise distances between sequences.

    Notes:
    ------
    - The distance is computed using the Needleman-Wunsch algorithm.
    - The matrix is symmetric, meaning `distance_matrix[i][j] == distance_matrix[j][i]`.
    - The function uses the `nw` function to compute alignment scores.
"""
  matD=[[0 for j in range(len(liste_de_sequences))] for i in range(len(liste_de_sequences))]
  for i in range(len(liste_de_sequences)):
    for j in range(len(liste_de_sequences)):
      if j>=i:
        score=nw([liste_de_sequences[i], liste_de_sequences[j]], blosum_m, gap_opening_score, gap_extension_score, identity_score=identity_score, substitution_score=substitution_score)[0]
        matD[i][j]=score
        matD[j][i]=matD[i][j]
  return matD


def inter_cluster_dist(list_indices1,list_indices2,distanceMatrix):
    """
    Computes the distance between two clusters.

    Parameters:
    ----------
    list_indices1 : list of int
        A list of indices representing the first cluster.
    list_indices2 : list of int
        A list of indices representing the second cluster.
    distance_matrix : list of list of float
        A distance matrix where `distance_matrix[i][j]` represents 
        the distance between the corresponding sequences.

    Returns:
    -------
    float
        The computed distance between the two clusters.

    Notes:
    ------
    - The distance is calculated as the average pairwise distance 
      between all elements in `list_indices1` and `list_indices2`.
    - The function assumes that `distance_matrix` is symmetric.
    - The result is computed as the sum of all pairwise distances 
      divided by the total number of pairs.
    """
    dist=0
    for i in list_indices1:
        for j in list_indices2:
            dist+=distanceMatrix[i][j]
    return dist/(len(list_indices1)*len(list_indices2))


def UPGMA(liste_de_sequences, blosum_m=False, gap_opening_score=-10, gap_extension_score=-2, identity_score=1, substitution_score=-1):
    """
    Constructs a hierarchical clustering tree using the UPGMA (Unweighted Pair Group Method with Arithmetic Mean) algorithm.

    Parameters:
    ----------
    sequence_list : list of str
        A list of sequences to cluster.
    blosum_m : bool, optional
        If True, uses the BLOSUM62 substitution matrix (default is False).
    gap_opening_score : int, optional
        The penalty for opening a gap in the alignment (default is -10).
    gap_extension_score : int, optional
        The penalty for extending an existing gap (default is -2).
    identity_score : int, optional
        The score for matching identical characters (default is 1).
    substitution_score : int, optional
        The score for aligning non-identical characters (default is -1).

    Returns:
    -------
    Tree
        A hierarchical clustering tree representing the relationships between sequences.

    Notes:
    ------
    - The algorithm starts with each sequence as an individual cluster.
    - The pairwise distance matrix is computed using the Needleman-Wunsch alignment.
    - At each step, the two closest clusters are merged into a new cluster.
    - The process continues until a single tree structure is formed.
    - The final tree can be visualized using `print_tree_sequences()`.
    """

    matD=mat_distances(liste_de_sequences, blosum_m, gap_opening_score, gap_extension_score, identity_score, substitution_score) #matrice de distance entre chaque paire de séquence
    clusters=[Tree(val=[i]) for i in range(len(liste_de_sequences))] #liste de cluster, chaque cluster étant un arbre
    maxd = 0
    while len(clusters)>1:
        max_i=0
        max_j=1
        for i in range(len(clusters)):
            for j in range(i+1,len(clusters)):
                dij=inter_cluster_dist(clusters[i].val,clusters[j].val,matD)
                if dij>maxd :
                    maxd=dij
                    max_i=i
                    max_j=j
        indices_new_cluster=clusters[max_i].val+clusters[max_j].val
        new_cluster=Tree(indices_new_cluster, left=clusters[max_i], right=clusters[max_j])
        clusters.pop(max_j)
        clusters.pop(max_i)
        clusters.append(new_cluster)
    tree_of_seq=clusters[0]
    # tree_of_seq.print_tree_sequences(liste_de_sequences)
    return tree_of_seq

#UPGMA(["CAT", "CHAT", "HER"], True, -8, -8)

