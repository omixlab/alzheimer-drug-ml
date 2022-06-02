while read PUBCHEM_BIOASSAY_ID; do
 
        for FEATURE_SET in descriptors morgan-1024 morgan-2048 mol2vec; do
            for ALGORITHM in rf extra_tree gradient_boosting logistic_regression; do
                mkdir -p data/train/$PUBCHEM_BIOASSAY_ID/$FEATURE_SET/$ALGORITHM
                cd data/train/$PUBCHEM_BIOASSAY_ID/$FEATURE_SET/$ALGORITHM 

		if [ -f $PUBCHEM_BIOASSAY_ID\_model.pickle ]; then
			echo "skipping $PUBCHEM_BIOASSAY_ID $FEATURE_SET $ALGORITHM"
			cd ../../../../../
			continue
		fi

                bambu-train \
                    --input-train ../../../../preprocessing/$PUBCHEM_BIOASSAY_ID/$FEATURE_SET/$PUBCHEM_BIOASSAY_ID\_preprocess_train.csv \
                    --output $PUBCHEM_BIOASSAY_ID\_model.pickle \
                    --time-budget 3600 \
                    --estimators $ALGORITHM \
                    --threads 1
 
                cd ../../../../../ 
                done
 
        done
 
done < data/raw/pubchem_ids
