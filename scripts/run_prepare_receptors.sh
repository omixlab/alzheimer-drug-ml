#!/usr/bin/env bash

eval "$(conda shell.bash hook)"

cd data/docking/
rm -r target/pdbqt/
mkdir -p target/pdbqt/

for pdb_file in target/pdb/*.pdb; do
    conda activate alzheimer-drug-ml
    pdbfixer $pdb_file
    conda activate alzheimer-drug-prepare-docking
    pdb_id=$(basename $pdb_file | cut -d'.' -f1)
    cat $pdb_file | grep ATOM > $pdb_file\_fix
    mv $pdb_file\_fix $pdb_file
    echo $pdb_id
    prepare_receptor4.py -r $pdb_file -o target/pdbqt/$pdb_id.pdbqt -A bonds_hydrogens > /dev/null
done
