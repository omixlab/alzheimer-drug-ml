cd data/raw

cat pubchem_ids | xargs -P 8 -i \
    bambu-download \
	--pubchem-assay-id {} \
	--pubchem-InchI-chunksize 100 \
	--output {}.csv
