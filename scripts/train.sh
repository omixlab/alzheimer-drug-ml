while read PUBCHEM_BIOASSAY_ID; do
 
        for FEATURE_SET in descriptors morgan-1024 morgan-2048 mol2vec; do
            for ALGORITHM in rf extra_tree gradient_boosting logistic_regression; do
                mkdir -p data/train/$PUBCHEM_BIOASSAY_ID/$FEATURE_SET/$ALGORITHM
                cd data/train/$PUBCHEM_BIOASSAY_ID/$FEATURE_SET/$ALGORITHM 
                bambu-train \
                    --input-train ../../../../preprocessing/$PUBCHEM_BIOASSAY_ID/$FEATURE_SET/$PUBCHEM_BIOASSAY_ID\_preprocess_train.csv \
                    --output $PUBCHEM_BIOASSAY_ID\_model.pickle \
                    --model-history \
                    --time-budget 3600 \
                    --estimators $ALGORITHM \
                    --threads 16
 
                cd ../../../../../ 
                done
 
        done
 
done < ../../../pubchem_ids