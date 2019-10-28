table=species.txt
for species in $(cat ${table} | cut -d ' ' -f1)
do
    n=$(grep ${species} ${table} | cut -d ' ' -f2)
    python3 generate_csv.py -s ${species} -n ${n}
done