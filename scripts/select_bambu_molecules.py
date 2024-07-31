#!/usr/bin/env bash

from rdkit import Chem
from rdkit.Chem import AllChem
import pandas as pd
import numpy as np
import requests
import glob
import os

for subset in ['fda', 'natural_compounds']:
    #os.system(f'rm -r data/docking/ligands/{subset}')
    os.system(f'mkdir -p data/docking/ligands/{subset}/pdb')
    os.system(f'mkdir -p data/docking/ligands/{subset}/pdbqt')
    zinc_ids = []
    df_data = pd.DataFrame()
    for mol in Chem.SDMolSupplier(f'data/screening/{subset}.sdf'):
        zinc_ids.append(mol.GetPropsAsDict()['zinc_id'])
    df_data['zinc_id'] = zinc_ids
    for csv_file in glob.glob(f'data/screening/bambu/{subset}/*.csv'):
        pubchem_id = csv_file.split('/')[4].split(".")[0]
        df_data[pubchem_id] = pd.read_csv(csv_file)['activity_probability']
    df_data = df_data[['zinc_id', '1239', '1276', '1285', '588852', '602250', '624040', '1468', '1347398', '1347395', 'B3DB']]
    df_data = df_data.query('B3DB >= 0.75')
    df_data.to_csv(f'data/screening/bambu.{subset}.csv', index=False)
    df_data.to_excel(f'data/screening/bambu.{subset}.xlsx', index=False)
    
    df_data['score'] = (df_data['B3DB'] + 1 - df_data['1239'] + (df_data['1276'] + df_data['1285']) / 2 + (df_data['588852'] + df_data['602250'] + df_data['624040']) / 3 + df_data['1468'] + (df_data['1347398'] + df_data['1347395']) / 2) / 6
    df_data = df_data.sort_values(by='score', ascending=False)
    
    for r, row in df_data.head(1000).iterrows():
        print(row.zinc_id)
        df_molecule_data = pd.read_csv(f'https://zinc15.docking.org/substances/{row.zinc_id}.csv')
        mol = Chem.MolFromSmiles(df_molecule_data['smiles'][0])
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol)
        try:
            AllChem.MMFFOptimizeMolecule(mol)
        except:
            pass
        with open(f'data/docking/ligands/{subset}/pdb/{row.zinc_id}.pdb', 'w') as writer:
            writer.write(Chem.MolToPDBBlock(mol))
        

