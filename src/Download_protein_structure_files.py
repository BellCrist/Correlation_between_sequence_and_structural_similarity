import os
import urllib.request


class StructureFileDownloader:

    @staticmethod
    def download_structure_file(pdb_id_list):
        cartella_pdb = "Storage/PDB_files"
        os.makedirs(cartella_pdb, exist_ok=True)

        # 2. Cicla sulla tua lista di pdb_id (la variabile 'labels' che hai creato prima)
        for pdb_id in pdb_id_list:
            # Assicurati che l'ID sia pulito (es. da spaziature) e di solito in maiuscolo
            pdb_id = str(pdb_id).strip().upper()

            # Il link ufficiale del database RCSB PDB per scaricare i file .pdb
            url = f"https://files.rcsb.org/download/{pdb_id}.pdb"

            # Il percorso dove verrà salvato il file nel tuo computer
            file_path = os.path.join(cartella_pdb, f"{pdb_id}.pdb")

            # 3. Scarica il file solo se non lo hai già fatto
            if not os.path.exists(file_path):
                try:
                    # Scarica il file dall'URL e salvalo
                    urllib.request.urlretrieve(url, file_path)
                    print(f"Scaricato con successo: {pdb_id}")
                except Exception as e:
                    # Cattura eventuali errori (es. se un PDB_ID non esiste o è obsoleto)
                    print(f"ERRORE: Impossibile scaricare {pdb_id}. Dettaglio: {e}")
            else:
                print(f"Già presente in locale: {pdb_id} (salto il download)")
