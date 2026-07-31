from Bio.Align import PairwiseAligner
from Bio.Align import substitution_matrices


class SmithWatermanScore:
    @staticmethod
    def sequence_identity(seq1, seq2):
        aligner = PairwiseAligner()

        # Smith–Waterman = allineamento locale
        aligner.mode = "local"

        # Matrice di sostituzione
        aligner.substitution_matrix = substitution_matrices.load("BLOSUM50")

        # Gap penalties
        aligner.open_gap_score = -12
        aligner.extend_gap_score = -2
        score = aligner.score(seq1, seq2)
        return score

    @staticmethod
    def normalized_score(score_dictionary, self_scores):
        normalized_scores = {}
        for (pdb1, pdb2), sw_score in score_dictionary.items():
            nsws = (2 * sw_score) / (self_scores[pdb1] + self_scores[pdb2])

            normalized_scores[tuple(sorted([pdb1, pdb2]))] = nsws
        return normalized_scores
