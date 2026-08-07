from Bio.Align import PairwiseAligner
from Bio.Align import substitution_matrices
import csv


class SmithWatermanScore:

    def __init__(self):
        self.self_scores = {}
        self.smith_waterman_scores = {}
        self.aligner = PairwiseAligner()

        # Smith–Waterman = allineamento locale
        self.aligner.mode = "local"

        # Matrice di sostituzione
        self.aligner.substitution_matrix = substitution_matrices.load("BLOSUM50")

        # Gap penalties
        self.aligner.open_gap_score = -12
        self.aligner.extend_gap_score = -2

    def sequence_identity(self, dataset_file_path):

        with open(dataset_file_path, mode='r') as file:
            file_rows = list(csv.DictReader(file))
            # Per ogni sequenza del dataset eseguo il confronto con tutte le altre
            for row in file_rows:
                seq1 = row['Sequence']
                pdb_id1 = row['PDB_ID']
                for record in file_rows:
                    seq2 = record['Sequence']
                    pdb_id2 = record['PDB_ID']
                    if pdb_id1 == pdb_id2:
                        # self score ha come chiave solo il pdbid1 perchè poi risulta utile così nel calcolo dello score normalizzato
                        self.self_scores[pdb_id1] = self.aligner.score(seq1, seq2)
                    else:
                        # Utilizzo tuple ordinate per evitare record duplicati che hanno stesso valore ma chiave invertita
                        self.smith_waterman_scores[tuple(sorted([pdb_id1, pdb_id2]))] = self.aligner.score(seq1, seq2)
        return self.self_scores, self.smith_waterman_scores

    @staticmethod
    def normalized_score(score_dictionary, self_scores):
        normalized_scores = {}
        for (pdb1, pdb2), sw_score in score_dictionary.items():
            nsws = (2 * sw_score) / (self_scores[pdb1] + self_scores[pdb2])

            normalized_scores[tuple(sorted([pdb1, pdb2]))] = nsws
        return normalized_scores
