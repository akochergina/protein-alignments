# protein-alignments
# Multiple Protein Sequence Alignment Project

## Overview
This project implements several multiple protein sequence alignment algorithms, including:
- **Needleman-Wunsch algorithm** for pairwise alignment.
- **Extensions of Needleman-Wunsch** for N+1 and N+M sequence alignments.
- **UPGMA clustering algorithm** for hierarchical clustering.
- **Evaluation using BALIbase benchmark and Sum-of-Pairs (SP) scoring.**

## Getting Started
### 1. Running the Example Notebook
To quickly explore the main functionalities of the project, run the **Jupyter notebook** `using_example.ipynb`.

### 2. Setting Up the BALIbase Dataset
For final testing and evaluation, place the **BALIbase dataset** in the root directory of the project:
```
/project_root/
    ├── balibase/
    ├── modules/
    ├── notebooks/
    ├── tests/
    ├── README.md
```

## Project Structure
```
/project_root/
    ├── modules/                    # Core Python modules for sequence alignment
    │   ├── needleman_wunsch.py    # Needleman-Wunsch and its extension N+1
    │   ├── nw_multiple.py          # Needleman-Wunsch extension N+M
    │   ├── multiple_align.py      # Core alignment algorithm 
    │   ├── upgma.py                # UPGMA clustering implementation
    │   ├── tree.py                 # tree class
    ├── balibase/          # BALIbase dataset (required for final testing)
    ├── notebooks/ 
    │   ├── using_example.ipynb # Example Jupyter Notebook
    ├── tests/              # Some test functions
    ├── README.md          # Project documentation
```

## Dependencies
Ensure the following Python libraries are installed:
```
Biopython
numpy
pandas
matplotlib
seaborn
```
Install them via pip:
```
pip install biopython numpy pandas matplotlib seaborn
```

## Evaluation
The accuracy of the alignment is evaluated using:
1. **BALIbase dataset** for reference alignments.
2. **Sum-of-Pairs (SP) score** to measure alignment similarity.


