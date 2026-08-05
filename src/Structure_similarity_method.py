import subprocess
import re
import sys


class StructureSimilarityMethods:

    @staticmethod
    def tm_align(path_pdb_id_1, path_pdb_id_2):
        comando = ["./TMalign", path_pdb_id_1, path_pdb_id_2]

        try:
            # 2. Eseguiamo il processo catturando l'output testuale
            risultato = subprocess.run(comando, capture_output=True, text=True, check=True)
            text_output = risultato.stdout

            # 3. Cerchiamo la riga del TM-score usando le regex
            # L'output di TM-align solitamente contiene righe come: "TM-score= 0.52345 (if normalized by...)"
            # Questa regex cerca "TM-score=", eventuali spazi, e poi cattura il numero decimale
            match_c1 = re.search(r"TM-score=\s*([0-9.]+)\s*\(if normalized by length of Chain_1", text_output)
            match_c2 = re.search(r"TM-score=\s*([0-9.]+)\s*\(if normalized by length of Chain_2", text_output)

            if match_c1 and match_c2:
                # Estraggo il numero e lo converto in float
                tm1 = float(match_c1.group(1))
                tm2 = float(match_c2.group(1))
                return tm1, tm2
            else:
                print(f"Attenzione: TM-score non trovato nell'output per {path_pdb_id_1} e {path_pdb_id_2}.")
                return None

        except subprocess.CalledProcessError as e:
            print(f"Errore nell'esecuzione di TM-align: {e.stderr}")
            return None
        except FileNotFoundError:
            print("Errore: Eseguibile TMalign non trovato. Verifica il percorso!")
            return None
