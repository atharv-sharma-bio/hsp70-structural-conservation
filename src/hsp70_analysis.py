import os
import urllib.request
import numpy as np
from Bio.PDB import PDBParser, Superimposer
from Bio.Align import PairwiseAligner
from Bio.SeqUtils import seq1

PDB_IDS = {
    "Human": "1S3X",
    "Yeast": "3D2F",
    "E_coli": "1DKG",
    "Bovine": "1CPR"
}

DATA_DIR = "data"

os.makedirs(DATA_DIR, exist_ok=True)

def fetch_pdb(pdb_id):
    """Fetch PDB structure file from RCSB PDB."""
    file_path = os.path.join(DATA_DIR, f"{pdb_id}.pdb")
    if not os.path.exists(file_path):
        url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
        print(f"Downloading {pdb_id} from RCSB...")
        urllib.request.urlretrieve(url, file_path)
    return file_path

def get_chain_ca_and_seq(pdb_path):
    """Parse C-alpha atoms and amino acid sequence from PDB file."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_path)
    model = structure[0]
    chain = list(model.get_chains())[0]

    ca_atoms = []
    seq_chars = []

    for residue in chain:
        if residue.get_id()[0] == ' ' and 'CA' in residue:
            ca_atoms.append(residue['CA'])
            try:
                seq_chars.append(seq1(residue.get_resname()))
            except Exception:
                seq_chars.append('X')

    return ca_atoms, "".join(seq_chars)

def run_analysis():
    print("=== HSP70 Structural Conservation Analysis Pipeline ===")

    ref_path = fetch_pdb(PDB_IDS["Human"])
    ref_ca, ref_seq = get_chain_ca_and_seq(ref_path)

    print(f"\nReference (Human 1S3X): {len(ref_seq)} residues, {len(ref_ca)} Cα atoms.")

    aligner = PairwiseAligner()
    aligner.mode = 'global'
    aligner.match_score = 1.0
    aligner.mismatch_score = 0.0
    aligner.open_gap_score = -1.0

    print("\nTarget | Sequence Identity (%) | Aligned Cα Pairs | RMSD (Å)")
    print("-" * 60)

    for organism, pdb_id in PDB_IDS.items():
        if organism == "Human":
            continue

        target_path = fetch_pdb(pdb_id)
        target_ca, target_seq = get_chain_ca_and_seq(target_path)

        alignments = aligner.align(ref_seq, target_seq)
        best_alignment = alignments[0]

        matches = sum(1 for a, b in zip(best_alignment[0], best_alignment[1]) if a == b and a != '-')
        max_len = max(len(ref_seq), len(target_seq))
        seq_identity = (matches / max_len) * 100

        ref_matched_ca = []
        target_matched_ca = []

        ref_idx = 0
        target_idx = 0

        for char_ref, char_target in zip(best_alignment[0], best_alignment[1]):
            if char_ref != '-' and char_target != '-':
                if ref_idx < len(ref_ca) and target_idx < len(target_ca):
                    ref_matched_ca.append(ref_ca[ref_idx])
                    target_matched_ca.append(target_ca[target_idx])
            if char_ref != '-':
                ref_idx += 1
            if char_target != '-':
                target_idx += 1

        superimposer = Superimposer()
        superimposer.set_atoms(ref_matched_ca, target_matched_ca)
        superimposer.apply(target_matched_ca)
        rmsd = superimposer.rms

        print(f"{organism:<8} ({pdb_id}) | {seq_identity:>20.2f}% | {len(ref_matched_ca):>16} | {rmsd:>8.3f} Å")

if __name__ == "__main__":
    run_analysis()
