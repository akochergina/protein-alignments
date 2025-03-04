import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from modules.multiple_align import align_multiple_sequences as ams
from Bio import pairwise2 

def test_multiple_align_basic():
    """
    Check Needleman-Wunsch alignment on a basic example.
    """
    list_of_seq=["CAT","CHAT","HER"]
    identity_score = 1
    substitution_score = -1
    gap_score = -2
    blosum_m=False
    gap_opening_score = gap_score
    gap_extension_score = gap_score

    alignments = ams(list_of_seq, blosum_m, gap_opening_score, gap_extension_score, identity_score, substitution_score)

    print(alignments)
    assert ("CHAT" in alignments and "C-AT" in alignments and "-HER" in alignments)
    print("All tests passed for needleman_wunsch_basic \n")

