import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from modules.needleman_wunsch import needleman_wunsch, score_i_j_alignment, score_i_j_alignment_multidim, print_alignments
from Bio import pairwise2 
from Bio.Align import substitution_matrices

def test_alignment_score_multdim():
    print("Running tests for score_i_j_alignment_multidim...")
    gap_score = -10
    assert score_i_j_alignment_multidim('A', ['A'], False, gap_score) == 1, f"Expected 1, got {score_i_j_alignment_multidim('A', ['A'], False, gap_score)}"
    assert score_i_j_alignment_multidim('A', ['A', 'B'], False, gap_score, substitution_score=0) == 0.5, f"Expected 0.5, got {score_i_j_alignment_multidim('A', ['A', 'B'], False, gap_score)}"
    assert score_i_j_alignment_multidim('A', ['A', 'B', 'C', 'D'], False, gap_score, substitution_score=0) == 0.25, f"Expected 0.25, got {score_i_j_alignment_multidim('A', ['A', 'B', 'C', 'D'], False, gap_score)}"
    
    assert score_i_j_alignment_multidim('A', ['A'], True, gap_score) == score_i_j_alignment('A', 'A', True), f"Expected {score_i_j_alignment('A', 'A', True)}, got {score_i_j_alignment_multidim('A', ['A'], True, gap_score)}"
    assert score_i_j_alignment_multidim('A', ['A', 'B'], True, gap_score) == (score_i_j_alignment('A', 'A', True) + score_i_j_alignment('A', 'B', True)) / 2, f"Expected {(score_i_j_alignment('A', 'A', True) + score_i_j_alignment('A', 'B', True)) / 2}, got {score_i_j_alignment_multidim('A', ['A', 'B'], True, gap_score)}"
    assert score_i_j_alignment_multidim('A', ['A', 'B', 'C', 'D'], True, gap_score) == (score_i_j_alignment('A', 'A', True) + score_i_j_alignment('A', 'B', True) + score_i_j_alignment('A', 'C', True) + score_i_j_alignment('A', 'D', True)) / 4, f"Expected {(score_i_j_alignment('A', 'A', True) + score_i_j_alignment('A', 'B', True) + score_i_j_alignment('A', 'C', True) + score_i_j_alignment('A', 'D', True)) / 4}, got {score_i_j_alignment_multidim('A', ['A', 'B', 'C', 'D'], True, gap_score)}"

    assert score_i_j_alignment_multidim('A', ['-'], False, gap_score) == -10
    assert score_i_j_alignment_multidim('A', ['-'], True, gap_score) == -10
    print("All tests passed for score_i_j_alignment_multidim \n")


def test_needleman_wunsch_basic():
    """
    Check Needleman-Wunsch alignment on a basic example.
    """
    print("Running tests for needleman_wunsch_basic...")
    sequences = ["ACGT", "AGT"]
    identity_score = 1
    substitution_score = -1
    gap_opening_score = -2
    gap_extension_score = -1

    score, alignments = needleman_wunsch(sequences, False, gap_opening_score, gap_extension_score, identity_score=identity_score, substitution_score=substitution_score)

    assert score == 1, f"Expected 1, got {score}"  
    assert len(alignments[0]) == len(alignments[1])
    print_alignments(alignments)
    print("All tests passed for needleman_wunsch_basic \n")

def test_needleman_wunsch_vs_biopython():
    """
    Compare our implementation of Needleman-Wunsch with Biopython.
    """
    print("Running tests for needleman_wunsch_vs_biopython...")
    sequences = ["ACGT", "AGT"]
    identity_score = 1
    substitution_score = -1
    gap_opening_score = -2
    gap_extension_score = -1

    
    score, alignments_nw = needleman_wunsch(sequences, False, gap_opening_score, gap_extension_score, identity_score=identity_score, substitution_score=substitution_score)

    # Biopython
    alignments = pairwise2.align.globalms(sequences[0], sequences[1], identity_score, substitution_score, gap_opening_score, gap_extension_score)
    biopython_score = alignments[0].score

    assert abs(score - biopython_score) < 1e-6
    assert alignments_nw[0] == alignments[0].seqA
    assert alignments_nw[1] == alignments[0].seqB
    print_alignments(alignments_nw)
    print(f"Biopython alignment is: \n{alignments[0].seqA}\n{alignments[0].seqB}")
    print("All tests passed for needleman_wunsch_vs_biopython \n")

def test_needleman_wunsch_vs_biopython_blossum():
    """
    Compare our implementation of Needleman-Wunsch with Biopython using BLOSUM62.
    """
    print("Running tests for needleman_wunsch_vs_biopython_blossum...")
    sequences = ["ACGT", "AGT"]
    gap_opening_score = -2
    gap_extension_score = -1
    blosum62 = substitution_matrices.load("BLOSUM62")

    score, alignments_nw = needleman_wunsch(sequences, True, gap_opening_score, gap_extension_score)

    # Biopython
    alignments = pairwise2.align.globalds(sequences[0], sequences[1], blosum62, gap_opening_score, gap_extension_score)
    biopython_score = alignments[0].score

    assert abs(score - biopython_score) < 1e-6
    assert alignments_nw[0] == alignments[0].seqA
    assert alignments_nw[1] == alignments[0].seqB
    print_alignments(alignments_nw)
    print(f"Biopython alignment is: \n{alignments[0].seqA}\n{alignments[0].seqB}")
    print("All tests passed for needleman_wunsch_vs_biopython_blossum \n")
    

test_alignment_score_multdim()
test_needleman_wunsch_basic()
test_needleman_wunsch_vs_biopython()
test_needleman_wunsch_vs_biopython_blossum()