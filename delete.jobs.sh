#/bin/bash

irun=""
frun=""

for RUNS in $(seq ${irun} ${frun})
do
	qdel ${RUNS}
done	
