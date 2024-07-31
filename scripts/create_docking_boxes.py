from biopandas.pdb import PandasPdb
import glob
import os

os.system('mkdir -p data/docking/boxes/')

for pdb_file in glob.glob('data/docking/sites/*.pdb'):

    pdb_df = PandasPdb().read_pdb(pdb_file)

    x_pos = round(pdb_df.df['ATOM'].x_coord.mean(), 2)
    y_pos = round(pdb_df.df['ATOM'].y_coord.mean(), 2)
    z_pos = round(pdb_df.df['ATOM'].z_coord.mean(), 2)

    x_size = round(abs(pdb_df.df['ATOM'].x_coord.max() - pdb_df.df['ATOM'].x_coord.min()) + 10, 2)
    y_size = round(abs(pdb_df.df['ATOM'].y_coord.max() - pdb_df.df['ATOM'].y_coord.min()) + 10, 2)
    z_size = round(abs(pdb_df.df['ATOM'].z_coord.max() - pdb_df.df['ATOM'].z_coord.min()) + 10, 2)
    
    with open(f'data/docking/boxes/{os.path.basename(pdb_file)}.config', 'w') as writer:
        writer.write(f'center_x = {x_pos}\n')
        writer.write(f'center_y = {y_pos}\n')
        writer.write(f'center_z = {z_pos}\n')
        writer.write(f'size_x = {x_size}\n')
        writer.write(f'size_y = {y_size}\n')
        writer.write(f'size_z = {z_size}\n')