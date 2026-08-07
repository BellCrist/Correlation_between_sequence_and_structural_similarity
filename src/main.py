import sys
import csv
from Fasta_Dataset_File_Conversion import FastaDatasetConversion
from Download_protein_structure_files import StructureFileDownloader
from Alignment_based import SmithWatermanScore
from Alignment_free import SequenceSimilarityAlignmentFree
import pandas as pd
from Structure_similarity_score import StructureSimilarity

if __name__ == '__main__':
    # print("Conversione file Fasta al csv...")
    # convertitore = FastaDatasetConversion()
    # if convertitore.fasta_to_csv_conversion():
    #     print("Conversione completata con successo")
    # else:
    #     print("Conversione non riuscita")

    smithWaterman = SmithWatermanScore()
    d2_distance = {}

# Calcolo della distanza d2 tra le sequenze del dataset
    print("Esecuzione algoritmo d2 distance...\n")
    df = pd.read_csv('Storage/Files/wood_pearson_dataset.csv')
    alignment_free = SequenceSimilarityAlignmentFree()
    d2_distance_scores = alignment_free.calcola_distanza_d2(df)

# Calcolo della similarità di sequenza tramite NSWS (All-against-all)
    print("Esecuzione algoritmo NSWS...\n")
    smithWaterman.sequence_identity()
    dataset_file_path = 'Storage/Files/wood_pearson_dataset.csv'
    self_scores, smith_waterman_scores = smithWaterman.sequence_identity(dataset_file_path)
    normalized_scores = smithWaterman.normalized_score(smith_waterman_scores, self_scores)

# Calcolo del TM-score per la similarità strutturale
    print("Calcolo similarità strutturale...")
    structure_similarity = StructureSimilarity()
    tm_scores = structure_similarity.calculcate_structural_similarity_score(alignment_free.get_labels())

    # Unione dei risultati ottenuti dentro un unico file csv
    print("Salvataggio risultati su file csv\n")
    rows = []
    for pair, nsws_score in normalized_scores.items():
        d2_score = d2_distance_scores.get(pair)
        tm_score = tm_scores.get(pair)
        rows.append({
            "PDB_ID_1": pair[0],
            "PDB_ID_2": pair[1],
            "NSWS": nsws_score,
            "D2": d2_score,
            "TM_score1": tm_score['tm_chain1'],
            "TM_score2": tm_score['tm_chain2']
        })

    df_results = pd.DataFrame(rows)
    df_results.to_csv("Storage/Files/sequence_similarity_results.csv", index=False)
    print("File di comparazione completato.")

# Download dei file .pdb dagli endpoint pdb, utili per la similarità strutturale
    # print("Inizio download file strutture...")
    # labels = alignment_free.get_labels()
    # structure_downloader = StructureFileDownloader()
    # structure_downloader.download_structure_file(labels)
    # print("Completato")
