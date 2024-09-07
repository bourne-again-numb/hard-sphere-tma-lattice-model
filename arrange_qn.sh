#/bin/bash

#total number of directories starting with temp_*
nreplicas=$(ls -d1 temp_*/ | wc -l)
let nreplicas-=1

rm -f ascending_qn.dat
for temp in $(seq -w 00 ${nreplicas});do
    CURRENT_REPLICA=temp_${temp}
    tail -n 1 ${CURRENT_REPLICA}/qn.c_mcsteps.csv >> ascending_qn_backup1.dat
done

#cat -n ascending_qn_backup1.dat > ascending_qn_backup2.dat
# awk '{
#        if(NR == 1) {
#            shift = $1
#        }

#        print ($1 - shift) "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7
# }' ascending_qn_backup2.dat > ascending_qn.dat

#total number of directories starting with temp_*
nlines=$(ls -d1 temp_*/ | wc -l)
file_temp="ascending_qn_backup3.dat"
data="exchange_info.csv"

tail -n ${nlines} $data > $file_temp

paste ${file_temp} ascending_qn_backup1.dat > ascending_qn.txt

sort -k2 -n ascending_qn.txt > ascending_qn.dat


rm ascending_qn_backup* ascending_qn.txt







#echo $nlines

# file_temp="sorted_exchange_info_temp.csv"
# file="sorted_exchange_info.csv"
# data="exchange_info.csv"

# tail -n ${nlines} $data > $file_temp

# sort -k2 -n $file_temp > $file
# rm $file_temp