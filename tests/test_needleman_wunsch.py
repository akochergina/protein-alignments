import pytest
import sys
import os
sys.path.append(os.path.abspath(os.path.join('..')))
from modules.needleman_wunsch import needleman_wunsch
from Bio import pairwise2  


def test_needleman_wunsch_basic():
    """
    Проверяем, корректно ли работает алгоритм на простом примере.
    """
    seq1 = "ACGT"
    seq2 = "AGT"
    identity_score = 1
    substitution_score = -1
    gap_score = -2

    score, align1, align2 = needleman_wunsch(seq1, seq2, identity_score, substitution_score, gap_score)

    # Ожидаемые результаты
    assert score == 1  # Примерное значение, зависит от логики твоего алгоритма
    assert len(align1) == len(align2)

def test_needleman_wunsch_vs_biopython():
    """
    Сравниваем наш алгоритм с Biopython.
    """
    seq1 = "ACGT"
    seq2 = "AGT"
    identity_score = 1
    substitution_score = -1
    gap_score = -2

    # Наш алгоритм
    score, align1, align2 = needleman_wunsch(seq1, seq2, identity_score, substitution_score, gap_score)

    # Biopython
    alignments = pairwise2.align.globalms(seq1, seq2, identity_score, substitution_score, gap_score, gap_score)
    biopython_score = alignments[0].score

    # Проверяем, что наши результаты совпадают
    assert abs(score - biopython_score) < 1e-6
