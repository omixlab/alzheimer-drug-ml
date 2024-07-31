import pandas as pd
import numpy as np
import json
import glob

rows = []

for json_file in glob.glob('data/validation/*/*/*/validation.json'):
    validation = json.loads(open(json_file).read())
    _,_,assay_id, descriptor, algorithm, filename = json_file.split('/')
    print(assay_id)
    row = {
            'assay_id': assay_id,
            'descriptors': descriptor,
            'algorithm': algorithm
    }
    for score in ['accuracy', 'recall', 'precision', 'f1', 'roc_auc']:
        default = f'%.2f'%(validation['raw_scores'][score][0])
        y_randomization_mean = np.mean(validation['raw_scores'][score][1::]) * 100
        y_randomization_stdv = np.std(validation['raw_scores'][score][1::])
        y_randomization_str  = f'%.2f +/- %.2f'%(y_randomization_mean, y_randomization_stdv)
        row[f'{score} (default)'] = default
        row[f'{score} (y-randomization)'] = y_randomization_str
    rows.append(row)

df_summary = pd.DataFrame(rows)

df_best_models = pd.DataFrame()

for assay_id in df_summary['assay_id'].unique():
    df_best_models = pd.concat([df_best_models, df_summary.query('assay_id == @assay_id').sort_values('roc_auc (default)', ascending=False).head(1)])

df_summary.to_csv('data/validation/summary.csv', index=False)
df_best_models.to_csv('data/validation/best_models.csv', index=False)

df_summary.to_excel('data/validation/summary.xlsx')
df_best_models.to_excel('data/validation/best_models.xlsx')