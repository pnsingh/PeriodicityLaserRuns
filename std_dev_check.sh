#!/bin/bash

filename=SPE_8.txt

i=$(wc -l < SPE_8.txt)
no_of_rows=$((i/3-1))

#echo $no_of_rows

end=$((no_of_rows*2))
start=$((end-9))

#echo $end
#echo $start

temp=$(head -n$end SPE_8.txt|tail -n10)

mean=0
for i in `seq 1 10`
do
	temp=$(head -n$start SPE_8.txt|tail -n1)
	mean=$((mean+temp))
	start=$((start+1))
done

#echo $temp