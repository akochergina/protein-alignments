import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from modules.needleman_wunsch import nw_multiple
from Bio import pairwise2 

print(nw_multiple.cost_function_block(['A','A','T'],['A']))