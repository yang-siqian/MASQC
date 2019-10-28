table=species_table.txt
for Species in $(cat ${table} | cut -f1 )
do
    echo ${Species}
    python3 threshold_analysis.py -s ${Species}
done
