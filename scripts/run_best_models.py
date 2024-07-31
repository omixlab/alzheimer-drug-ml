#!/usr/bin/env python

import multiprocessing as mp
import pandas as pd
import os

subsets = ['fda', 'natural_compounds']

commands = []

for subset in subsets:
    for r, row in pd.read_csv('data/validation/best_models.csv').iterrows():
        model_path        = os.path.join('data', 'train', str(row.assay_id), row.descriptors, row.algorithm, str(row.assay_id) + '_model.pickle')
        preprocessor_path = os.path.join('data', 'preprocessing', str(row.assay_id), row.descriptors, str(row.assay_id) + '_preprocessor.pickle')
        csv_output        = os.path.join('data', 'screening', 'bambu', subset, str(row.assay_id) + '.csv')
        if os.path.isfile(csv_output):
            continue
        commands.append(f'bambu-predict --input data/screening/{subset}.sdf --model {model_path} --preprocessor {preprocessor_path} --output {csv_output}')

pool = mp.Pool(10)

for _ in pool.map(os.system, commands):
    _
