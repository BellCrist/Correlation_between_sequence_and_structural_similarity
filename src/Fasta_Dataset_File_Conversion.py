from Bio import SeqIO
import pandas as pd
import re


class FastaDatasetConversion:
    @staticmethod
    def fasta_to_csv_conversion():
        # File FASTA di input
        fasta_file = "../Storage/Files/wood_pearson dataset.fasta"

        records = []

        for record in SeqIO.parse(fasta_file, "fasta"):

            header = record.description

            # ============================
            # Estrazione PDB ID
            # ============================
            pdb_match = re.match(r"([A-Za-z0-9]{4})_\d+", header)

            if pdb_match:
                pdb_id = pdb_match.group(1).upper()
            else:
                pdb_id = ""

            # ============================
            # Estrazione catena
            # ============================
            chain_match = re.search(r"Chain\s+([A-Za-z0-9])", header)

            if chain_match:
                chain = chain_match.group(1)
            else:
                chain = ""

            # ============================
            # Suddivisione dell'header
            # ============================
            parts = header.split("|")

            protein = parts[2].strip() if len(parts) > 2 else ""
            organism = parts[3].strip() if len(parts) > 3 else ""

            records.append({

                "PDB_ID": pdb_id,
                "Chain": chain,
                "Protein": protein,
                "Organism": organism,
                "Sequence": str(record.seq)

            })

        # Creazione DataFrame
        df = pd.DataFrame(records)
        # Salvataggio CSV
        df.to_csv("wood_pearson_dataset.csv", index=False)
        print(df.head())
        print(f"\nNumero di proteine: {len(df)}")
