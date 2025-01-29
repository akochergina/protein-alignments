import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from modules.needleman_wunsch import needleman_wunsch
from Bio import pairwise2 


def test_needleman_wunsch_basic():
    """
    Check Needleman-Wunsch alignment on a basic example.
    """
    seq1 = "ACGT"
    seq2 = "AGT"
    identity_score = 1
    substitution_score = -1
    gap_score = -2

    score, align1, align2 = needleman_wunsch(seq1, seq2, identity_score, substitution_score, gap_score)

    assert score == 1  
    assert len(align1) == len(align2)
    print(f"Alignment is: \n{align1}\n{align2}")

def test_needleman_wunsch_vs_biopython():
    """
    Compare our implementation of Needleman-Wunsch with Biopython.
    """
    seq1 = "ACGT"
    seq2 = "AGT"
    identity_score = 1
    substitution_score = -1
    gap_score = -2

    
    score, align1, align2 = needleman_wunsch(seq1, seq2, identity_score, substitution_score, gap_score)

    # Biopython
    alignments = pairwise2.align.globalms(seq1, seq2, identity_score, substitution_score, gap_score, gap_score)
    biopython_score = alignments[0].score

    assert abs(score - biopython_score) < 1e-6
    assert align1 == alignments[0].seqA
    assert align2 == alignments[0].seqB
    print(f"Alignment is: \n{align1}\n{align2}")
    print(f"Biopython alignment is: \n{alignments[0].seqA}\n{alignments[0].seqB}")


test_needleman_wunsch_basic()