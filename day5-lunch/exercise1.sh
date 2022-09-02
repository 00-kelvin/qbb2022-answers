cut -f 5,6 -d "," aau1043_dnm.csv | sort | uniq -c | grep mother | cut -f 1 -d "," | sort -k2,2 > maternally_inherited.txt 

cut -f 5,6 -d "," aau1043_dnm.csv | sort | uniq -c | grep father | cut -f 1 -d "," | sort -k2,2 > paternally_inherited.txt 

join -1 2 -2 2 paternally_inherited.txt maternally_inherited.txt >> father_mother_joined.txt

echo "Proband_id father_count mother_count Father_age Mother_age" > wrangled.txt
join <(sort -k 1,1 father_mother_joined.txt) <(sed 's/,/ /g' aau1043_parental_age.csv | sort -k 1,1) >> wrangled.txt