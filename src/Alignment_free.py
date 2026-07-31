import pandas as pd
import numpy as np
import sys
from itertools import product
from alfpy.utils.seqrecords import SeqRecords
from alfpy import word_pattern, word_vector
from alfpy import word_distance
from alfpy import word_d2


class SequenceSimilarityAlignmentFree:

    @staticmethod
    def calcola_distanza_alfpy(df_csv, k=3):
        """
        Calcola la matrice delle distanze D2 tra le sequenze contenute nel DataFrame
        utilizzando le classi WordPattern e D2 della libreria ALFPY.
        """

        sequence = ''
        sequence_list = []
        labels = []

        # 1. Parsing delle sequenze e creazione degli oggetti Seq di alfpy
        for idx, row in df_csv.iterrows():
            pdb_id = (row['PDB_ID'])
            sequence = str(row['Sequence']).strip().upper()

            # Istanza della classe Seq di alfpy
            sequence_list.append(sequence)
            labels.append(pdb_id)
        # seq_obj = SeqRecords(labels, sequence_list)

        seq_lengths = [len(s) for s in sequence_list]
        # 2. Utilizzo della classe WordPattern per l'estrazione dei k-mer
        # k: lunghezza della parola amminoacidica (es. 3 per tripeptidi)
        wp = word_pattern.create(sequence_list, k)
        word_vector_object = word_vector.Counts(seq_lengths, wp)

        # 3. Utilizzo della classe D2 per il calcolo della matrice di distanza
        d2_obj = word_d2.Distance(word_vector_object.data)
        #TODO ciclare il dataset e per ogni coppia di sequenze (da passare al metodo pwdist_d2) ci salviamo lo score
        alignment_free_scores = {}
        n = len(labels)
        for i in range(n):
            for j in range(i+1, n):
                score = d2_obj.pwdist_d2(i, j)
                alignment_free_scores[tuple(sorted([labels[i], labels[j]]))] = score
                print(alignment_free_scores)
                sys.exit()

        # # 4. Organizzazione dei risultati in un DataFrame Pandas con etichette
        # df_dist = pd.DataFrame(dist_matrix, index=labels, columns=labels)
        # return df_dist
