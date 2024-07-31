while read PUBCHEM_BIOASSAY_ID; do
 
        for FEATURE_SET in descriptors morgan-1024 morgan-2048 mol2vec; do
            for ALGORITHM in rf extra_tree gradient_boosting logistic_regression; do
                mkdir -p data/validation/$PUBCHEM_BIOASSAY_ID/$FEATURE_SET/$ALGORITHM
                cd data/validation/$PUBCHEM_BIOASSAY_ID/$FEATURE_SET/$ALGORITHM 

		#if [ -f $PUBCHEM_BIOASSAY_ID\_model.pickle ]; then
		#	echo "skipping $PUBCHEM_BIOASSAY_ID $FEATURE_SET $ALGORITHM"
		#	cd ../../../../../
		#	continue
		#fi

                bambu-validate  \
                    --input-train ../../../../preprocessing/$PUBCHEM_BIOASSAY_ID/$FEATURE_SET/$PUBCHEM_BIOASSAY_ID\_preprocess_train.csv \
                    --input-test ../../../../preprocessing/$PUBCHEM_BIOASSAY_ID/$FEATURE_SET/$PUBCHEM_BIOASSAY_ID\_preprocess_test.csv \
                    --model ../../../../train/$PUBCHEM_BIOASSAY_ID/$FEATURE_SET/$ALGORITHM/$PUBCHEM_BIOASSAY_ID\_model.pickle \
		    --output validation.json \
		    --randomizations 100
 
                cd ../../../../../ 
                done
 
        done
 
done < data/raw/pubchem_ids
