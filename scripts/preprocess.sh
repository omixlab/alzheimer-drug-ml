while read PUBCHEM_BIOASSAY_ID; do
 
        for FEATURE_SET in descriptors morgan-1024 morgan-2048 mol2vec; do
 
                mkdir -p data/preprocessing/$PUBCHEM_BIOASSAY_ID/$FEATURE_SET
                cd data/preprocessing/$PUBCHEM_BIOASSAY_ID/$FEATURE_SET
 
                bambu-preprocess --input ../../../raw/$PUBCHEM_BIOASSAY_ID.csv --output $PUBCHEM_BIOASSAY_ID\_preprocess.csv --output-preprocessor $PUBCHEM_BIOASSAY_ID\_preprocessor.pickle --feature-type $FEATURE_SET --mol2vec-model-path ../../../../mol2vec.pickle --train-test-split 0.75 --undersample
                cd ../../../..
 
        done
 
done < data/raw/pubchem_ids
