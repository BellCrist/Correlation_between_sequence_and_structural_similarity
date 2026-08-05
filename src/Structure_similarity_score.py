import os
import sys

from Structure_similarity_method import StructureSimilarityMethods
from src.utils.Unique_key_generator import unique_key_generator


class StructureSimilarity:

    @staticmethod
    def calculcate_structural_similarity_score(labels):
        # Qui prendiamo tutti i path dei files dentro PDB_filese li utilizziamo per ottenere lo score
        # di similarità tramite TM-align
        structure_similarity_method = StructureSimilarityMethods()
        tm_align_data = {}

        for i in range(len(labels)):
            for j in range(i + 1, len(labels)):
                pdb_1 = labels[i]
                pdb_2 = labels[j]

            # Se la coppia esiste già, anche invertita, non rieseguo il tm-align
                key_pair = unique_key_generator(pdb_1, pdb_2)
                if key_pair in tm_align_data:
                    continue

                pdb_1, pdb_2 = key_pair
                file_path1 = f"Storage/PDB_files/{pdb_1}.pdb"
                file_path2 = f"Storage/PDB_files/{pdb_2}.pdb"

                if os.path.exists(file_path1) and os.path.exists(file_path2):
                    tm_score1, tm_score2 = structure_similarity_method.tm_align(file_path1, file_path2)
                    if tm_score1 is not None and tm_score2 is not None:
                        # Salviamo entrambi i valori nel dizionario nidificato
                        tm_align_data[(pdb_1, pdb_2)] = {
                            "tm_chain1": tm_score1,
                            "tm_chain2": tm_score2
                        }

        return tm_align_data
