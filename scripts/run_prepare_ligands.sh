#!/usr/bin/env bash

eval "$(conda shell.bash hook)"
conda activate alzheimer-drug-prepare-docking

cd data/docking/ligands/

for subset in fda natural_compounds; do
    cd $subset/pdb
    for pdb_file in *.pdb; do
        pdb_id=$(basename $pdb_file | cut -d'.' -f1)
        prepare_ligand4.py -l $pdb_file -o ../pdbqt/$pdb_id.pdbqt -A bonds_hydrogens
    done
    cd ../../

done