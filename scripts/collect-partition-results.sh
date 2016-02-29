#!/bin/bash

BASEDIR=/scratch/daint/biddisco/egpgv/scripts

for NODES in 32 64 128 256 512 1024
do
  for NPERNODE in 1
  do
    for REDUCTION in 1 2 4 8 16 32 64 128 256 512 1024 
    do
      for FILTER in 0 1
      do
        JOB_NAME=$(printf 'partition-%04d-%d-%06d-%d' ${NODES} ${NPERNODE} ${REDUCTION} ${FILTER})
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

          printf "mem $MEM PPF $PPF \n"
        fi
        cd $BASEDIR

      done
    done
  done
done
