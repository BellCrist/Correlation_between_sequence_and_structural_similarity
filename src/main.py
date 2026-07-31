import sys
import csv
from FastaDatasetFileConversion import FastaDatasetConversion
from Alignment_based import SmithWatermanScore
from Alignment_free import SequenceSimilarityAlignmentFree
import pandas as pd

if __name__ == '__main__':
    # convertitore = FastaDatasetConversion()
    # convertitore.fasta_to_csv_conversion()

    # Algoritmo Smith Waterman all-against-all
    smithWaterman = SmithWatermanScore()
    self_scores = {}
    scores = {}

    df = pd.read_csv('Storage/Files/wood_pearson_dataset.csv')
    alignment_free = SequenceSimilarityAlignmentFree()
    alignment_free.calcola_distanza_alfpy(df)
    sys.exit()

    with open('Storage/Files/wood_pearson_dataset.csv', mode='r') as file:
        fileRows = list(csv.DictReader(file))
        for row in fileRows:
            seq1 = row['Sequence']
            pdbId1 = row['PDB_ID']
            for record in fileRows:
                seq2 = record['Sequence']
                pdbId2 = record['PDB_ID']
                if pdbId1 == pdbId2:
                    print('Stessa sequenza')
                    self_scores[pdbId1] = smithWaterman.sequence_identity(seq1, seq2)
                else:
                    # Utilizzo tuple ordinate per evitare record duplicati che hanno stesso valore ma chiave invertita
                    scores[tuple(sorted([pdbId1, pdbId2]))] = smithWaterman.sequence_identity(seq1, seq2)

        normalized_scores = smithWaterman.normalized_score(scores, self_scores)

        # Preparazione algoritmo kmer
