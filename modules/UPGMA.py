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
  """renvoie une matrice de distances deux à deux des séquences de la liste
  Entrée : une liste de string qui sont les séquences à comparer
  Sortie : une matrice de distances deux à deux des séquences de la liste
            la distance est calculée en utilisant needelman wunsh"""
  matD=[[0 for j in range(len(liste_de_sequences))] for i in range(len(liste_de_sequences))]
  for i in range(len(liste_de_sequences)):
    for j in range(len(liste_de_sequences)):
      if j>=i:
        score=nw([liste_de_sequences[i], liste_de_sequences[j]], blosum_m, gap_opening_score, gap_extension_score, identity_score=identity_score, substitution_score=substitution_score)[0]
        matD[i][j]=score
        matD[j][i]=matD[i][j]
  return matD


def inter_cluster_dist(list_indices1,list_indices2,distanceMatrix):
    """renvoie la distance entre le cluster 1 et 2
    Entrée : deux clusters, étant des listes d'indices
            une matrice de distance qui associe aux indices concernées la distance entre les séquences correspondantes
    Sortie : float, la distance entre les deux clusters"""
    dist=0
    for i in list_indices1:
        for j in list_indices2:
            dist+=distanceMatrix[i][j]
    return dist/(len(list_indices1)*len(list_indices2))


def UPGMA(liste_de_sequences, blosum_m, gap_opening_score, gap_extension_score, identity_score=1, substitution_score=-1):
    """renvoie l'arbre de clusters associé à la liste de séquence
    entrée : une liste de séquence et la fonction de cout de l'alignement
    sortie : un arbre de clusters"""

    matD=mat_distances(liste_de_sequences, blosum_m, gap_opening_score, gap_extension_score, identity_score, substitution_score) #matrice de distance entre chaque paire de séquence
    clusters=[Tree(val=[i]) for i in range(len(liste_de_sequences))] #liste de cluster, chaque cluster étant un arbre
    mind = float('inf')
    while len(clusters)>1:
        min_i=0
        min_j=1
        for i in range(len(clusters)):
            for j in range(i+1,len(clusters)):
                dij=inter_cluster_dist(clusters[i].val,clusters[j].val,matD)
                if dij<mind :
                    mind=dij
                    min_i=i
                    min_j=j
        indices_new_cluster=clusters[min_i].val+clusters[min_j].val
        new_cluster=Tree(indices_new_cluster, left=clusters[min_i], right=clusters[min_j])
        clusters.pop(min_j)
        clusters.pop(min_i)
        clusters.append(new_cluster)
    tree_of_seq=clusters[0]
    tree_of_seq.print_tree()
    return tree_of_seq

