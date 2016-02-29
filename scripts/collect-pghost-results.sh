#!/bin/bash

BASEDIR=/scratch/daint/biddisco/egpgv/scripts

NPERNODE=1
FILTER=0

for NODES in 32 64 128 256 512 1024
do
  for REDUCTION in 1
  do
      for GHOST in 0 0.01 
      do
        JOB_NAME=$(printf 'pghosts-%04d-%06d-%s' ${NODES} ${REDUCTION} ${GHOST})
        JOB_STRING=$(echo $JOB_NAME | tr '-' ' ')
        cd $JOB_NAME
        if [ -e timing-log.txt ] ; then
          printf "$JOB_STRING " 
          MEM=$(grep "Total memory" timing-log.txt | awk '{print $(NF)}' | tr -d '\n')

          if [ -z "$MEM" ]; then
            MEM=" 0.0 "
          fi

          PPF=$(grep "Execute vtkParticlePartitionFilter" timing-log.txt | grep -v head | awk '{print $(NF-1)}' | sort -nr | head -n1)
          D3=$(grep "Execute vtkDistributedDataFilter" timing-log.txt | grep -v head |awk '{print $(NF-1)}' | sort -nr | head -n1)

          if [ -z "$PPF" ]; then
            PPF=$(grep "Execute vtkDistributedDataFilter" timing-log.txt | grep -v head |awk '{print $(NF-1)}' | sort -nr | head -n1)
          fi
          if [ -z "$PPF" ]; then
            PPF="0.0"
          fi

          GP=$(grep "Num Ghost Points " timing-log.txt | awk '{print $(NF-2)}' | sed 's/L,//g' )

          printf "mem $MEM PPF $PPF GP $GP \n"
        fi
        cd $BASEDIR

      done
  done
done
