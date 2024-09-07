#/bin/bash

VERSION='v6.2'
RUN='run'
MASTER=$(pwd)/
INPUT="input.dat"
#JOB_QUEUE="pbs.script"
JOB_QUEUE="lsf.script"

#since the dissociation constant for the reaction is very high
#reaction: SN + OH- <=> SI + H2O   Kd = 1.7*10^6
#there would be equal of the amounts of SI and sda
#in this script we don't distingiush between SI and sda 

#---------variables-------
i_snsda="-0.10";  	d_snsda="-0.10";	f_snsda=${i_snsda}
i_teos="5";     	d_teos="5";     	f_teos=${i_teos}
i_sda="1";     	        d_sda="4";   	        f_sda=${i_sda}
#i_h2o="9500";  	d_h2o="1";   	        f_h2o=${i_h2o}

manual_concen_flag="1"
#here we disregard the concentration of water and just enter TEOS and sda molecule number.
#the amount of ionic silica produces is the same as sda, cause we assume that the dissociation
#constant of the above reaction is too high

i_t="0.15";    	d_t="0.15";    	f_t=${i_t}
i_x="4";     	d_x="4";       	f_x=${i_x}
i_y="4";     	d_y="1";       	f_y=${i_y}
i_z="4";     	d_z="1";       	f_z=${i_z}

sda_size="1"

nprocs="4"
let nprocs-=1

runs="1"

ftemp="0.30"
jump_flag="0"

ipen3="0.00";     dpen3="0.10";     fpen3=${ipen3}
ipen4="0.00";     dpen4="0.10";     fpen4=${ipen4}

nsweeps="10000000";nprint="50000";nexchange="50000"
tstarf="0.6d0";m0="0.3";m1="0.35";m2="0.5"
start_config="0";time_limit="6.0"

#testing variables
nsweeps="200000";nprint="50";nexchange="50"
time_limit="0.014"
#-------------------------

# make the directory structure
for pen3 in $(seq ${ipen3} ${dpen3} ${fpen3});do #3 ring penalty

    for pen4 in $(seq ${ipen4} ${dpen4} ${fpen4});do #4 ring penalty

	for snsda in $(seq ${i_snsda} ${d_snsda} ${f_snsda});do #SNSDA interaction
	    
	    for nsda in $(seq ${i_sda} ${d_sda} ${f_sda});do #SDA concentration
		
		for nteos in $(seq ${i_teos} ${d_teos} ${f_teos});do #TEOS concentration
		    
		    for t in $(seq ${i_t} ${d_t} ${f_t});do #temperature
			
			for cux in $(seq ${i_x} ${d_x} ${f_x});do #x length
			    for cuy in $(seq ${i_y} ${d_y} ${f_y});do #y length
				for cuz in $(seq ${i_z} ${d_z} ${f_z});do #z length
				    
				    CURRENT_VERSION=${VERSION}_eSNSDA${snsda}_nTEOS${nteos}_nSDAOH${nsda}_3RP${pen3}_4RP${pen4}_itstar${t}_ftstar${ftemp}_${cux}x${cuy}x${cuz}
				    mkdir -p ${CURRENT_VERSION}
				    
			        #enter the current version directory
				    cd ${CURRENT_VERSION}
				    
			        # once inside the current version directiory do the following
			        # copy the xcutable,*.mod,pbs file from the master directory
			        # create the input.dat 
			        # create the pbs.submit file
				    cp -f ${MASTER}xcutable .
				    cp -f ${MASTER}*.mod .
				    cp -f ${MASTER}${JOB_QUEUE} ./${JOB_QUEUE}
				    #cp -f ${MASTER}*.gnu .
				    cp -f ${MASTER}analysis.sh .
				    cp -f ${MASTER}sort_exchange_info.sh .
				    
			            # -------------input.dat file creation---------------
				    echo "${cux},${cuy},${cuz}" > ${INPUT}
				    echo "${t}d0,${sda_size}" >> ${INPUT}
				    echo "${nsweeps},${nprint},${nexchange}" >> ${INPUT}
				    echo "${nteos},${nsda},${nsda},${manual_concen_flag}" >> ${INPUT}
				    echo "${pen3}d0,${pen4}d0" >> ${INPUT}
				    echo "${snsda}d0,-1.00d0" >> ${INPUT}
				    echo "${tstarf},${m0},${m1},${m2}" >> ${INPUT}
				    echo "${start_config}","${time_limit}" >> ${INPUT}
				    echo "${ftemp}","${nexchange}","${jump_flag}" >> ${INPUT}
			            # ---------------------------------------------------
				    
				    for j in $(seq 1 ${runs});do #3 ring penalty
					
					CURRENT_RUN=${RUN}_${j}
					
					mkdir -p ${CURRENT_RUN}
					
					cp -f ${MASTER}xcutable ${CURRENT_RUN}/
					cp -f ${MASTER}${CURRENT_VERSION}/input.dat ${CURRENT_RUN}/
					cp -f ${MASTER}*.mod ${CURRENT_RUN}/
					cp -f ${MASTER}${JOB_QUEUE} ${CURRENT_RUN}/
					#cp -f ${MASTER}*.gnu ${CURRENT_RUN}/
					cp -f ${MASTER}analysis.sh ${CURRENT_RUN}/
					cp -f ${MASTER}sort_exchange_info.sh ${CURRENT_RUN}/
					
					cd ${CURRENT_RUN}
					
			                # -------------edit the job.submit file--------------
					sed -i s/abcdef/"${CURRENT_VERSION}_${j}"/g ${MASTER}${CURRENT_VERSION}/${CURRENT_RUN}/${JOB_QUEUE}
			                # ---------------------------------------------------			    
					for temp in $(seq -w 00 ${nprocs});do
					    CURRENT_TEMP=temp_${temp}
					    mkdir -p ${CURRENT_TEMP}
					    cp analysis.sh input.dat ${CURRENT_TEMP}
					done
				        #submit the job
                                        #qsub ${JOB_QUEUE}
				        #bsub < ${JOB_QUEUE}
					let nprocs+=1; ~/bin/mpich2/bin/mpirun -np ${nprocs} ./xcutable
				        #./xcutable
					
                                        # exit the run directory
		       		        #cd ${MASTER}/${CURRENT_VERSION}
					
					cd ${MASTER}${CURRENT_VERSION}
					
				    done
				    
				    cd ${MASTER}
				    
				done #z length
			    done #y length
			done #z length
		    done #temperature
		done #TEOS number
	    done #SDA number
	done # SN-SDA interaction
    done
done
