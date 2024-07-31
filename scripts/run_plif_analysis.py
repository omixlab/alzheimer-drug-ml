import xml.etree.ElementTree as ET
from biopandas.pdb import PandasPdb
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
import pandas as pd
import pickle
import glob
import os
import shlex
import json

target_blacklist = ['2ACE', '1ACL', '1BVN', '1BRU', '6QAA', '6XBY', '3FP3', '3MGN', '1EVE', '4F5S', '1DX6', '1VOT']
target_list = ['2UZS', '1Z0Q', '5I3V', '2J5D', '4TPK', '2W4O', '2I1M', '3ANR', '4ACD', '1E3G', '3BKL', '2Z5X', '4A79', '5K4I', '5ES1', '2MZ7', '2FV5']
target_list = [v.lower() for v in target_list]
def format_residue(residue):
    residue, position, chain = residue.split(':')
    return f"{residue.lower().capitalize()}{position} ({chain.upper()})"

residues_data = {}

for subset in ['fda', 'natural_compounds']:
    print(subset)
    plif_data = []
    for plip_xml in tqdm(glob.glob(f'data/docking/results/{subset}/*/*/plip/report.xml')):
        _,_,_,_,target,zinc_id,_,_ = plip_xml.split('/')
        if target not in target_list:
            continue
        site = []
        if not os.path.isfile(f'data/docking/sites/{target}.pdb'):
            continue
        for r, row in PandasPdb().read_pdb(f'data/docking/sites/{target}.pdb').df['ATOM'][['chain_id', 'residue_name', 'residue_number']].drop_duplicates().iterrows():
            site.append(f'{row.residue_name}:{row.residue_number}:{row.chain_id}')
        tree = ET.parse(plip_xml)
        root = tree.getroot()
        residues = []
        for interaction in root.iterfind('bindingsite/interactions/*/*'):
            residue = {}
            for attribute in interaction.getchildren():
                if attribute.tag == 'resnr':
                    residue['residue_number'] = attribute.text
                if attribute.tag == 'restype':
                    residue['residue_name'] = attribute.text
                if attribute.tag == 'reschain':
                    residue['chain_id'] = attribute.text
            residues.append(f"{residue['residue_name']}:{residue['residue_number']}:{residue['chain_id']}")
        residues = list(set(residues))
        residues_data[(subset, target, zinc_id)] = residues
        plif_score = len(set(site).intersection(set(residues))) / len(set(site))
        plif_data.append({'target': target, 'zinc_id': zinc_id, 'plif_score': plif_score})
    df = pd.DataFrame(plif_data)
    df = df.pivot(index='zinc_id', columns=['target'], values=['plif_score'])
    df = df.reset_index()
    df = df.fillna(0)
    df['mean plif score'] = df[df.columns[1::]].mean(axis=1)
    df = df.sort_values(by='mean plif score', ascending=False)
    df.to_csv(f'data/docking/{subset}.csv')

    targets = list(df.drop(['zinc_id', 'mean plif score'], axis=1).columns)
    targets = [target[1] for target in targets]

    rows = []

    for r, row in df.iterrows():
        row_data = row.to_dict()
        row_data['high-affinity-targets'] = []
        row_data['high-affinity-interactions'] = {}
        row_data['high-affinity-targets-LE'] = {}
        row_data['high-affinity-plif-values'] = {}
        for target in target_list:
            for column in df.columns:
                if column[0] == 'plif_score' and column[1] in target_list:
                    row_data['high-affinity-plif-values'][column[1]] = row[column]
        for target in targets:
            zinc_id  = row.zinc_id[0]
            zinc_mol = AllChem.MolFromPDBFile(f'data/docking/ligands/{subset}/pdb/{zinc_id}.pdb') 
            heavy_count = Descriptors.HeavyAtomCount(zinc_mol)
            try:
                vina_result = float(shlex.split(open(f'data/docking/results/{subset}/{target}/{zinc_id}/result_ligand_1.pdbqt').readlines()[0])[3])
                if vina_result < (-7.0) and (vina_result / heavy_count) < -0.3:
                    row_data['high-affinity-targets'].append([target, vina_result])
                    row_data['high-affinity-interactions'][target] = residues_data[(subset, target, zinc_id)]
                    row_data['high-affinity-interactions'][target] = ', '.join(list(map(format_residue, row_data['high-affinity-interactions'][target])))
                    row_data['high-affinity-targets-LE'][target] = (vina_result / heavy_count)
            except:
                vina_result = None
        row_data['targets'] = []
        for target in targets:
            if target in [t[0] for t in row_data['high-affinity-targets']] and target in row_data['high-affinity-interactions'] and row_data['high-affinity-interactions'][target]:
                row_data['targets'].append(target)
        
        row_data['high-affinity-targets'] = json.dumps(dict(row_data['high-affinity-targets']))
        row_data['high-affinity-targets-LE'] = json.dumps(dict(row_data['high-affinity-targets-LE']))
        row_data['high-affinity-plif-values'] = json.dumps(dict([(key, value) for (key, value) in row_data['high-affinity-plif-values'].items() if key in row_data['targets']]))
        rows.append(row_data)
    df_with_targets = pd.DataFrame(rows)
    df_with_targets = df_with_targets[(df_with_targets['high-affinity-targets'] != {}) & (df_with_targets['high-affinity-interactions'] != {})]
    
    for column in df_with_targets.columns:
        if column[0] == 'plif_score':
            del df_with_targets[column]

    with open(f'data/docking/{subset}.with_binding_energies.pandas', 'wb') as writer:
        writer.write(pickle.dumps(df_with_targets))
    
    df_with_targets.to_csv(f'data/docking/{subset}.with_binding_energies.csv')
    df_with_targets.to_excel(f'data/docking/{subset}.with_binding_energies.xlsx')

    with open(f'data/docking/{subset}.with_binding_energies.top_100.pandas', 'wb') as writer:
        writer.write(pickle.dumps(df_with_targets.head(100)))
    
    df_with_targets.head(100).to_csv(f'data/docking/{subset}.top_100.with_binding_energies.csv')
    df_with_targets.head(100).to_excel(f'data/docking/{subset}.top_100.with_binding_energies.xlsx')

    rows = []
    
    for r, row in df_with_targets.iterrows():
        for target in row['targets']:
            new_row = {}
            new_row['zinc_id'] = row[("zinc_id", '')]
            new_row['target'] = target
            new_row['affinity'] = json.loads(row['high-affinity-targets'])[target]
            new_row['ligand_efficacy'] = json.loads(row['high-affinity-targets-LE'])[target]
            new_row['plif'] = json.loads(row['high-affinity-plif-values'])[target]
            new_row['residues'] = row['high-affinity-interactions'][target]
            print(new_row)
            rows.append(new_row)
    
    df_with_targets_lines = pd.DataFrame(rows)

    with open(f'data/docking/lines_{subset}.with_binding_energies.pandas', 'wb') as writer:
        writer.write(pickle.dumps(df_with_targets_lines))
    
    df_with_targets_lines.to_csv(f'data/docking/lines_{subset}.with_binding_energies.csv')
    df_with_targets_lines.to_excel(f'data/docking/lines_{subset}.with_binding_energies.xlsx')

    with open(f'data/docking/lines_{subset}.with_binding_energies.top_100.pandas', 'wb') as writer:
        writer.write(pickle.dumps(df_with_targets_lines.head(100)))
    
    df_with_targets_lines.head(100).to_csv(f'data/docking/lines_{subset}.top_100.with_binding_energies.csv')
    df_with_targets_lines.head(100).to_excel(f'data/docking/lines_{subset}.top_100.with_binding_energies.xlsx')
