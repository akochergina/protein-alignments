import pandas
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from Bio.Align import substitution_matrices


def score_i_j_alignment(i: chr, j: chr, blosum_m: bool, identity_score=1, substitution_score=-1):
    """
    Calculate the score of aligning characters i and j. If blosum_m is True, we use BLOSUM62 matrix.

    Parameters:
    ----------
    i : chr
        First character to align.
    j : chr
        Second character to align.
    blosum_m : bool
        If True, we use BLOSUM62 matrix.
    identity_score : int
        Score for aligning identical characters.
    substitution_score : int
        Score for aligning non-identical characters.

    Returns:
    -------
    score : int
        The score of aligning i and j.
    """
    if blosum_m:
        matrix = substitution_matrices.load("BLOSUM62")
        return int(matrix[(i, j)])

    if i == j:
        return identity_score
    else:
        return substitution_score


def plot_nw_matrix(matrix, arrow_matrix, seq1, seq2):
    """
    Visualize the Needleman-Wunsch matrix with arrows.

    Parameters:
    ----------
    matrix : pandas.DataFrame
        The filled matrix with scores.
    arrow_matrix : pandas.DataFrame
        The matrix of arrows indicating traceback paths.
    seq1 : str
        The first sequence (used for row labels).
    seq2 : str
        The second sequence (used for column labels).
    """
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Convert matrix to numpy array, replacing None with 0
    matrix_np = matrix.fillna(0).to_numpy()
    
    # Define labels for the axes (insert '-' at the beginning to account for initial gap)
    row_labels = ['-'] + list(seq1)
    col_labels = ['-'] + list(seq2)
    
    # Create a heatmap with custom labels
    sns.heatmap(matrix_np.astype(float), annot=True, fmt=".0f", cmap="Blues", linewidths=0.5, 
                ax=ax, cbar=False, xticklabels=col_labels, yticklabels=row_labels)
    
    ax.xaxis.set_ticks_position("top")
    ax.xaxis.set_label_position("top")

    #traceback index[0] = n, m
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

def print_nw_result(matrix, arrow_matrix, score, alignment1, alignment2, seq1, seq2):
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
    alignment1 : str
        The first aligned sequence.
    alignment2 : str
        The second aligned sequence.
    '''
    print(f"Alignment was made with Needleman-Wunsch algorithm. Score is {score}.")
    print(f"Alignment is: \n{alignment1}\n{alignment2}")
    plot_nw_matrix(matrix, arrow_matrix, seq1, seq2)
    return

def fill_needleman_wunsch_matrix(seq1, seq2, blosum_m, gap_opening_score, gap_extension_score, identity_score=1, substitution_score=-1):
    '''
    Fill the Needleman-Wunsch matrix and store the indexes of the arrows with possibility to have up to 3 arrows.

    Parameters:
    ----------
    seq1 : str
        First sequence to align.
    seq2 : str
        Second sequence to align.
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
    rows = len(seq1) + 1
    cols = len(seq2) + 1
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
            if blosum_m:
                match = matrix.at[i-1, j-1] + score_i_j_alignment(seq1[i-1], seq2[j-1], blosum_m)
            else:
                match = matrix.at[i-1, j-1] + score_i_j_alignment(seq1[i-1], seq2[j-1], blosum_m, identity_score, substitution_score)
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

def needleman_wunsch(seq1, seq2, blosum_m, gap_opening_score, gap_extension_score, identity_score=1, substitution_score=-1):
    """
    Perform Needleman-Wunsch alignment.

    Parameters:
    ----------
    seq1 : str
        First sequence to align.
    seq2 : str
        Second sequence to align.
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
    alignment : tuple
        A tuple containing the aligned sequences and a score.
    """
    if blosum_m:
        matrix, arrow_matrix = fill_needleman_wunsch_matrix(seq1, seq2, blosum_m, gap_opening_score, gap_extension_score)
    else:
        matrix, arrow_matrix = fill_needleman_wunsch_matrix(seq1, seq2, blosum_m, gap_opening_score, gap_extension_score, identity_score, substitution_score)
    
    score = matrix.at[len(seq1), len(seq2)]

    alignement1 = ''
    alignement2 = ''
    i = len(seq1)
    j = len(seq2)

    while i > 0 or j > 0:
        for prev_i, prev_j in arrow_matrix.at[i, j]:
            if i - prev_i == 1 and j - prev_j == 1:
                alignement1 = seq1[i-1] + alignement1
                alignement2 = seq2[j-1] + alignement2
            elif i - prev_i == 1:
                alignement1 = seq1[i-1] + alignement1
                alignement2 = '-' + alignement2
            else:
                alignement1 = '-' + alignement1
                alignement2 = seq2[j-1] + alignement2
            i = prev_i
            j = prev_j
    
    print_nw_result(matrix, arrow_matrix, score, alignement1, alignement2, seq1, seq2)

    return score, alignement1, alignement2

