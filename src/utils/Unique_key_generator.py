def unique_key_generator(pdb_1, pdb_2):
    # 'sorted' ordina alfabeticamente i due ID
    # tuple() la converte in tupla (necessaria come chiave per i dizionari)
    return tuple(sorted([str(pdb_1).strip().upper(), str(pdb_2).strip().upper()]))
