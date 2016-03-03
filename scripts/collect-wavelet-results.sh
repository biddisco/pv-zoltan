#!/bin/bash

BASEDIR=$(pwd)

for NODES in 32 64 128 256 512 1024 2048 4096
do
  for NPERNODE in 1
  do
    for RESOLUTION in 100 200
    do
      for FILTER in 0 1
      do
        JOB_NAME=$(printf 'wavelet-%04d-%03d-%d' ${NODES} ${RESOLUTION} ${FILTER})
        JOB_STRING=$(echo $JOB_NAME | tr '-' ' ')
        cd $JOB_NAME
        if [ -e timing-log.txt ] ; then
          printf "$JOB_STRING "
          MEM=$(grep "numcells" timing-log.txt | awk '{print $(NF)}' | tr -d '\n')

          if [ -z "$MEM" ]; then
            MEM=" 0.0 "
          fi

          PPF=$(grep "Execute vtkMeshPartitionFilter" timing-log.txt | grep -v head | awk '{print $(NF-1)}' | sort -nr | head -n1)
          D3=$(grep "Execute vtkDistributedDataFilter" timing-log.txt | grep -v head |awk '{print $(NF-1)}' | sort -nr | head -n1)

          STILL=$(grep "STILL" timing-log.txt | awk '{print $(NF-1)}' | sort -nr | head -n1)

          if [ -z "$PPF" ]; then
            PPF=$D3
          fi
          if [ -z "$PPF" ]; then
            PPF="0.0"
          fi

          printf "cells $MEM PPF $PPF still $STILL\n"
        fi
        cd $BASEDIR

      done
    done
  done
done
