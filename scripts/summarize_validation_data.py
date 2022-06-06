
from pathlib import Path
import pandas as pd
import json
import os
import glob
import warnings
warnings.simplefilter('ignore')

df = pd.DataFrame(columns=['assay_id', 'feature_set', 'algorithm', 'accuracy', 'recall', 'precision', 'f1', 'roc_auc'])
df_best = pd.DataFrame(columns=['assay_id', 'feature_set', 'algorithm', 'accuracy', 'recall', 'precision', 'f1', 'roc_auc', 'model_path', 'processor_path'])

for json_file in glob.glob('data/validation/*/*/*/validation'):

    print(json_file)
    
    path = Path(json_file)
    assay_id    = path.parts[2]
    feature_set = path.parts[3]
    algorithm   = path.parts[4]

    validation_results = json.loads(open(json_file).read())

    accuracy  = round(float(validation_results['raw_scores']['accuracy'][0])*100, 2)
    recall    = round(float(validation_results['raw_scores']['recall'][0])*100, 2)
    precision = round(float(validation_results['raw_scores']['precision'][0])*100, 2)
    f1        = round(float(validation_results['raw_scores']['f1'][0])*100, 2)
    roc_auc   = round(float(validation_results['raw_scores']['roc_auc'][0])*100, 2)

    df = df.append({
        'assay_id': assay_id,
        'feature_set': feature_set,
        'algorithm': algorithm,
        'accuracy': accuracy,
        'recall': recall,
        'precision': precision,
        'f1': f1,
        'roc_auc': roc_auc
    }, ignore_index=True)

df.to_csv('data/validation/summary.csv', index=False)

for assay_id in df.assay_id.unique():
    row_data = df.query('assay_id == @assay_id').sort_values('roc_auc', ascending=False).iloc[0]
    row_data_dict = dict(row_data)
    row_data['model_path'] = f'data/train/{row_data.assay_id}/{row_data.feature_set}/{row_data.algorithm}/{row_data.assay_id}_model.pickle'
    row_data['processor_path'] = f'data/preprocessing/{row_data.assay_id}/{row_data.feature_set}/{row_data.assay_id}_preprocessor.pickle'
    df_best = df_best.append(row_data, ignore_index=True)

df_best.to_csv('data/validation/best_models.csv', index=False)
