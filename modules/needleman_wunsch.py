import pandas

def score_i_j_alignment(i: chr, j: chr, identity_score, substitution_score, gap_score):
    """
    Calculate the score of aligning characters i and j.

    Parameters:
    ----------
    i : chr
        First character to align.
    j : chr
        Second character to align.
    identity_score : int
        Score for aligning identical characters.
    substitution_score : int
        Score for aligning non-identical characters.
    gap_score : int
        Score for introducing a gap.

    Returns:
    -------
    score : int
        The score of aligning i and j.
    """
    if i == j:
        return identity_score
    elif i == '-' or j == '-':
        return gap_score
    else:
        return substitution_score

def fill_needleman_wunsch_matrix(seq1, seq2, identity_score, substitution_score, gap_score):
    '''
    Fill the Needleman-Wunsch matrix and store the indexes of the arrows with possibility to have up to 3 arrows.

    Parameters:
    ----------
    seq1 : str
        First sequence to align.
    seq2 : str
        Second sequence to align.
    identity_score : int
        Score for aligning identical characters.
    substitution_score : int
        Score for aligning non-identical characters.
    gap_score : int
        Score for introducing a gap.

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

    arrow_matrix.at[0, 0] = None
    matrix.at[0, 0] = 0

    # Initialize the first row and column
    for i in range(1, rows):
        matrix.at[i, 0] = i * gap_score
        arrow_matrix.at[i, 0] = [(i-1, 0)]
    for j in range(1, cols):
        matrix.at[0, j] = j * gap_score
        arrow_matrix.at[0, j] = [(0, j-1)]

    # Fill the matrix 
    for i in range(1, rows):
        for j in range(1, cols):
            # Calculate the scores
            match = matrix.at[i-1, j-1] + score_i_j_alignment(seq1[i-1], seq2[j-1], identity_score, substitution_score, gap_score)
            delete = matrix.at[i-1, j] + gap_score
            insert = matrix.at[i, j-1] + gap_score

            # Update the matrix
            max_score = max(match, delete, insert)
            arrow_matrix.at[i, j] = []
            matrix.at[i, j] = max_score
            if max_score == match:
                arrow_matrix.at[i, j].append((i-1, j-1))
            if max_score == delete:
                arrow_matrix.at[i, j].append((i-1, j))
            if max_score == insert:
                arrow_matrix.at[i, j].append((i, j-1))

    return matrix, arrow_matrix

def needleman_wunsch(seq1, seq2, identity_score, substitution_score, gap_score):
    """
    Perform Needleman-Wunsch alignment.

    Parameters:
    ----------
    seq1 : str
        First sequence to align.
    seq2 : str
        Second sequence to align.
    identity_score : int
        Score for aligning identical characters.
    substitution_score : int
        Score for aligning non-identical characters.
    gap_score : int
        Score for introducing a gap.

    Returns:
    -------
    alignment : tuple
        A tuple containing the aligned sequences and a score.
    """
    matrix, arrow_matrix = fill_needleman_wunsch_matrix(seq1, seq2, identity_score, substitution_score, gap_score)

    print(matrix)
    print(arrow_matrix)
    
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

    return score, alignement1, alignement2

# Test the function
seq1 = 'ACGT'
seq2 = 'AGT'
identity_score = 1
substitution_score = -1
gap_score = -1

needleman_wunsch(seq1, seq2, identity_score, substitution_score, gap_score)