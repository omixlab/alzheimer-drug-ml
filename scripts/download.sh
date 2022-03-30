cd data/raw
while read pubchem_id; do 
    bambu-download \
	--pubchem-assay-id $pubchem_id \
	--pubchem-InchI-chunksize 100 \
	--output $pubchem_id.csv
done < pubchem_ids