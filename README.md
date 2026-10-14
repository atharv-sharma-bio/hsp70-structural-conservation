# HSP70 Structural Conservation Analysis Pipeline

[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![BioPython](https://img.shields.io/badge/library-BioPython-green.svg)](https://biopython.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview
This repository contains an automated computational biology pipeline designed to evaluate sequence divergence and 3D backbone structural variation (RMSD) across evolutionary homologs of the **HSP70 Nucleotide-Binding Domain (NBD)**.

The research manuscript based on this codebase is currently under preparation for submission to the *Journal of Emerging Investigators* (JEI).

## Key Research Findings
* **Hypothesis (RMSD < 2.5 Å):** **Rejected.** 
* **Biological Cause:** Elevated global RMSD values (9.915 Å - 19.878 Å) stem from intrinsic NBD conformational dynamics—specifically subdomain hinge flexing during ATP hydrolysis and open-state trapping induced by nucleotide exchange factors—rather than permanent evolutionary fold degradation.

## Summary Results

| Organism | Taxa | PDB ID | Sequence Identity (%) | Aligned Cα Pairs | 3D Backbone RMSD (Å) |
| :--- | :--- | :--- | :--- | :--- | :--- |
| **Human** | Eukaryote | `1S3X` | 100.00% (Ref) | 382 | 0.000 |
| **Yeast** | Eukaryote | `3D2F` | 99.74% | 377 | 9.915 |
| **Bacterium** | Prokaryote | `1DKG` | 54.57% | 82 | 19.571 |
| **Bovine** | Eukaryote | `1CPR` | 19.37% | 129 | 19.878 |

## Mathematical & Algorithmic Implementation

1. **Global Pairwise Sequence Alignment:** Implemented via the **Needleman–Wunsch algorithm** using `Bio.Align.PairwiseAligner`.
2. **Sequence-Guided 3D Rigid-Body Superposition:** Executed using the **Kabsch algorithm** (`Bio.PDB.Superimposer`) across matched Cα atom pairs.

## Installation & Usage

1. **Clone the repository:**
   ```bash
   git clone [https://github.com/YOUR_USERNAME/hsp70-structural-conservation.git](https://github.com/YOUR_USERNAME/hsp70-structural-conservation.git)
   cd hsp70-structural-conservation
