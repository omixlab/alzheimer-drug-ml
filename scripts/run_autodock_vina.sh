#!/usr/bin/env bash

eval "$(conda shell.bash hook)"
conda activate alzheimer-drug-vina

rm -r data/docking/results/

for subset in fda natural_compounds; do
    for receptor_pdbqt in data/docking/target/pdbqt/*.pdbqt; do
        pdb_id=$(basename $receptor_pdbqt | cut -d'.' -f1)
        config_file=data/docking/boxes/$pdb_id.pdb.config
        for ligand_pdbqt in data/docking/ligands/$subset/pdbqt/*.pdbqt; do
            zinc_id=$(basename $ligand_pdbqt | cut -d'.' -f1)
            result_directory=data/docking/results/$subset/$pdb_id/$zinc_id
            mkdir -p $result_directory
            mkdir -p $result_directory/plip
            score=$(vina \
                --receptor $receptor_pdbqt \
                --ligand $ligand_pdbqt \
                --config $config_file \
                --out $result_directory/result.pdbqt \
                --cpu 16 \
                --exhaustiveness 20 > /dev/null
            cat $result_directory/result.pdbqt | grep 'REMARK VINA RESULT:' | head -1 | python -c "import shlex; print(shlex.split(input())[3])" \
            )
            vina_split --input $result_directory/result.pdbqt
            obabel -ipdbqt $result_directory/result_ligand_1.pdbqt -opdb -O $result_directory/result.pdb
            cat data/docking/target/pdb/$pdb_id.pdb $result_directory/result.pdb > $result_directory/complex.pdb
            plip -f $result_directory/complex.pdb -o $result_directory/plip/ --txt --xml
            echo -e "$pdb_id\t$zinc_id\t$score"
        done 
    done
done