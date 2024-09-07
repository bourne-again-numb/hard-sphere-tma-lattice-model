#/bin/bash

#total number of directories starting with temp_*
nlines=$(ls -d1 temp_*/ | wc -l)

#echo $nlines

file_temp="sorted_exchange_info_temp.csv"
file="sorted_exchange_info.csv"
data="exchange_info.csv"

tail -n ${nlines} $data > $file_temp

sort -k2 -n $file_temp > $file
rm $file_temp