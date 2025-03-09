import pandas
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from Bio.Align import substitution_matrices


""" 
Here the goal is to align two blocs of strings already aligned. Its the needleman wunsh algorithm adapted.
"""


def cost_function_block(list1, list2, blosum_m=False, identity_score=1, substitution_score=-1):
    """ 
    returns the cost of aligning two character lists
    Character lists represent the already aligned characters of the blocks
    Parameters:
    ------------
    list1 : list of already aligned caracters 1
    list2 : list of already aligned caracters 2
    blosum_m : bool
        If True, we use BLOSUM62 matrix.
    identity_score : int
        Score for aligning identical characters.
    substitution_score : int
        Score for aligning non-identical characters.

    Returns:
    Score : int, associated score
    """
    score=0
    n=len(list1)
    m=len(list2)

    if blosum_m:
        matrix = substitution_matrices.load("BLOSUM62")

    for i in range(n):
        for j in range(m):
            if blosum_m:
                if list1[i]!='-' and list2[j]!='-':
                    score+= int(matrix[(list1[i],list2[j])])
                if list1[i]=='-' and list2[j]=='-':
                    score+=identity_score
                else :
                    score+=0
            else :
                if list1[i] == list2[j]:
                    score+= identity_score
                else:
                    score+= substitution_score
    return score/(n*m)

#print(cost_function_block(['H', '-'], ['E'], True))


def plot_nw_matrix_multiple(matrix, arrow_matrix, block1, block2):
    """
    Visualize the Needleman-Wunsch matrix with arrows.

    Parameters:
    ----------
    matrix : pandas.DataFrame
        The filled matrix with scores.
    arrow_matrix : pandas.DataFrame
        The matrix of arrows indicating traceback paths.
    block1 : list of str
        List containing sequences already aligned
    block 2 : list of str
        List containing sequences already aligned. Objective is now to align block1 and block2
    """
    fig, ax = plt.subplots(figsize=(6, 4))
    
    # Convert matrix to numpy array, replacing None with 0
    matrix_np = matrix.fillna(0).to_numpy()
    
    # Define labels for the axes (insert '-' at the beginning to account for initial gap)
    row_labels = ['-'] + [ [block1[k][i] for k in range (len(block1))] for i in range (len(block1[0])) ]
    col_labels = ['-'] + [ [block2[k][i] for k in range (len(block2))] for i in range (len(block2[0])) ]
    
    # Create a heatmap with custom labels
    
    sns.heatmap(matrix_np.astype(float), annot=True, fmt=".0f", cmap="Blues", linewidths=0.5, 
                ax=ax, cbar=False, xticklabels=col_labels, yticklabels=row_labels)

    ax.xaxis.set_ticks_position("top")
    ax.xaxis.set_label_position("top")

    traceback_indexes = [(matrix.shape[0]-1, matrix.shape[1]-1)]
    while traceback_indexes[-1] != (0, 0):
        for prev_i, prev_j in arrow_matrix.at[traceback_indexes[-1][0], traceback_indexes[-1][1]]:
            traceback_indexes.append((prev_i, prev_j))

    # Draw arrows based on arrow_matrix
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            if arrow_matrix.at[i, j] is not None:
                for prev_i, prev_j in arrow_matrix.at[i, j]:
                    dx = prev_j - j  # X direction
                    dy = prev_i - i  # Y direction

                    # color is red if arrow is in traceback_indexes
                    if ((prev_i, prev_j) in traceback_indexes and (i, j) in traceback_indexes):
                        color = 'red'
                    else:
                        color = 'black'
                    if dx == 0:
                        ax.arrow(j + 0.5, i + 0.35, dx * 0.5, dy * 0.5, head_width=0.1, 
                                 head_length=0.1, fc=color, ec=color)
                    elif dy == 0:
                        ax.arrow(j + 0.25, i + 0.5, dx * 0.4, dy * 0.5, head_width=0.1, 
                                 head_length=0.1, fc=color, ec=color)
                    else:
                        ax.arrow(j + 0.35, i + 0.35, dx * 0.6, dy * 0.6, head_width=0.1, 
                             head_length=0.1, fc=color, ec=color)
    
    ax.set_title("Needleman-Wunsch Alignment Matrix with Traceback Arrows")
    plt.show()

def print_alignments(alignments):
    """
    Print the aligned sequences.

    Parameters:
    ----------
    alignments : list of str
        List containing the aligned sequences.
    """
    print("Alignments:")
    for alignment in alignments:
        print(alignment)
    return

def print_nw_result_multiple(matrix, arrow_matrix, score, alignments, block1, block2):
    '''
    Print the Needleman-Wunsch result.

    Parameters:
    ----------
    matrix : pandas.DataFrame
        The filled matrix.
    arrow_matrix : pandas.DataFrame
        The matrix of arrows. In matrix_arrows[i, j] we store the indexes of the cells where we can go from cell i, j.
    score : int
        The score of the alignment.
    alignments : list of str
        The aligned sequences.
    block1 : list of str
        List containing sequences already aligned
    block 2 : list of str
        List containing sequences already aligned. Objective is now to align block1 and block2
    '''
    print(f"Alignment was made with Needleman-Wunsch algorithm. Score is {score}.")
    print_alignments(alignments)
    plot_nw_matrix_multiple(matrix, arrow_matrix, block1, block2)
    return

def fill_needleman_wunsch_matrix_multiple(block1, block2, blosum_m, gap_opening_score, gap_extension_score, identity_score=1, substitution_score=-1):
    '''
    Fill the Needleman-Wunsch matrix for multiple sequences, regrouped in two already aligned blocks
    and store the indexes of the arrows with possibility to have up to 3 arrows.

    Parameters:
    ----------
    block1 : list of str
        List containing sequences already aligned
    block 2 : list of str
        List containing sequences already aligned. Objective is now to align block1 and block2
    blosum_m : bool
        If True, we use BLOSUM62 matrix.
    gap_opening_score : int
        Score for opening a gap.
    gap_extension_score : int
        Score for extending a gap.
    identity_score : int
        Score for aligning identical characters.
    substitution_score : int
        Score for aligning non-identical characters.

    Returns:
    -------
    matrix : pandas.DataFrame
        The filled matrix.
    '''
    # Initialize the matrix
    rows = len(block1[0]) + 1
    cols = len(block2[0]) + 1
    matrix = pandas.DataFrame(index=range(rows), columns=range(cols))
    arrow_matrix = pandas.DataFrame(index=range(rows), columns=range(cols))
    gap_matrix = pandas.DataFrame(index=range(rows), columns=range(cols))

    arrow_matrix.at[0, 0] = None
    matrix.at[0, 0] = 0
    gap_matrix.at[0, 0] = 0
 

    # Initialize the first row and column
    for i in range(1, rows):
        matrix.at[i, 0] = gap_opening_score + (i-1) * gap_extension_score
        arrow_matrix.at[i, 0] = [(i-1, 0)]
        gap_matrix.at[i, 0] = 1
    for j in range(1, cols):
        matrix.at[0, j] = gap_opening_score + (j-1) * gap_extension_score
        arrow_matrix.at[0, j] = [(0, j-1)]
        gap_matrix.at[0, j] = 1

    # Fill the matrix 
    for i in range(1, rows):
        for j in range(1, cols):
            # Calculate the scores
            match = matrix.at[i-1, j-1] + cost_function_block([(block1[k][i-1]) for k in range(len(block1))],[(block2[k][j-1]) for k in range (len(block2))], blosum_m, identity_score, substitution_score)
            if gap_matrix.at[i-1, j] == 1:
                delete = matrix.at[i-1, j] + gap_extension_score
            else:
                delete = matrix.at[i-1, j] + gap_opening_score
            if gap_matrix.at[i, j-1] == 1:
                insert = matrix.at[i, j-1] + gap_extension_score
            else:
                insert = matrix.at[i, j-1] + gap_opening_score

            # Update the matrix
            max_score = max(match, delete, insert)
            arrow_matrix.at[i, j] = []
            matrix.at[i, j] = max_score
            if max_score == match:
                arrow_matrix.at[i, j].append((i-1, j-1))
                gap_matrix.at[i, j] = 0
            if max_score == delete:
                arrow_matrix.at[i, j].append((i-1, j))
                gap_matrix.at[i, j] = 1
            if max_score == insert:
                arrow_matrix.at[i, j].append((i, j-1))
                gap_matrix.at[i, j] = 1

    return matrix, arrow_matrix

def needleman_wunsch_multiple(block1, block2, blosum_m, gap_opening_score=-10, gap_extension_score=-2, print_result=False, identity_score=1, substitution_score=-1):
    """
    Perform Needleman-Wunsch alignment, to align two blocks of already aligned sequences

    Parameters:
    ----------
    block1 : list of str
        List containing sequences already aligned
    block 2 : list of str
        List containing sequences already aligned. Objective is now to align block1 and block2
    blosum_m : bool
        If True, we use BLOSUM62 matrix.
    gap_opening_score : int
        Score for opening a gap.
    gap_extension_score : int
        Score for extending a gap.
    print_result : bool
        If True, print the matrix.
    identity_score : int
        Score for aligning identical characters.
    substitution_score : int
        Score for aligning non-identical characters.

    Returns:
    -------
    alignment : tuple
        A tuple containing the aligned sequences and a score.
    """
    i = len(block1[0])
    j = len(block2[0])
    if i < j :
        block1, block2 = block2, block1
        i, j = j, i

    if blosum_m:
        matrix, arrow_matrix = fill_needleman_wunsch_matrix_multiple(block1, block2, blosum_m, gap_opening_score, gap_extension_score)
    else:
        matrix, arrow_matrix = fill_needleman_wunsch_matrix_multiple(block1, block2, blosum_m, gap_opening_score, gap_extension_score, identity_score, substitution_score)
    
    score = matrix.at[len(block1[0]), len(block2[0])]


    aa = [("") for i in range (len(block1))]
    bb = [("") for i in range (len(block2))]

    while i > 0 or j > 0:
        prev_i, prev_j = arrow_matrix.at[i, j][0]
        if i - prev_i == 1 and j - prev_j == 1:
            for k in range(len(block1)):
                if i > 0:  
                    aa[k] += block1[k][i-1]  
                else:
                    aa[k] += '-'
            for k in range(len(block2)):
                if j > 0:
                    bb[k] += block2[k][j-1]
                else:
                    bb[k] += '-'
        elif i - prev_i == 1:
            for k in range (len(block1)):
                aa[k]+=(block1[k][i-1])
            for k in range(len(block2)):
                bb[k]+='-'
        else:
            for k in range (len(block1)):
                aa[k]+='-'
            for k in range(len(block2)):
                bb[k]+=(block2[k][j-1])
        i = prev_i
        j = prev_j
    
    
    for k in range (len(block1)):
      aa[k]=aa[k][::-1]
    for k in range(len(block2)):
      bb[k]=bb[k][::-1]
    alignments = aa+bb

    if print_result:
        print_nw_result_multiple(matrix, arrow_matrix, score, alignments, block1, block2)

    return score, alignments
