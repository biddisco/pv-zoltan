#!/bin/bash

# This function writes a slurm script. 
# We can call it with different parameter 
# settings to create different experiments

function write_script
{
TASKS=$[$NPERNODE * $NODES]
FILE_PREFIX=$(printf '%d' ${REDUCTION})
JOB_NAME=$(printf 'partition-%04d-%d-%06d-%d' ${NODES} ${NPERNODE} ${REDUCTION} ${FILTER})
DIR_NAME=$(printf '%s' ${JOB_NAME})

if [ -f $DIR_NAME/timing-log.txt ] ; then
  echo "$DIR_NAME/timing-log.txt already exists, skipping..."
  return 0
else
  echo "Creating job $DIR_NAME"
fi

mkdir -p $DIR_NAME

cat << _EOF_ > ${DIR_NAME}/slurm-exp.bash
#!/bin/bash
#SBATCH --job-name=${JOB_NAME}
#SBATCH --output=slurm.out
#SBATCH --error=slurm.err
#SBATCH --partition=${QUEUE}
#SBATCH --nodes=${NODES}
#SBATCH --ntasks-per-node=${NPERNODE}
#SBATCH --distribution=cyclic
#SBATCH --time=01:00:00
#SBATCH --exclusive
#SBATCH --constraint=startx

export DISPLAY=:0.0
#export MPICH_MAX_THREAD_SAFETY=multiple

module load boost/1.56.0
module load cudatoolkit
#module load viz

# boost
export LD_LIBRARY_PATH=/apps/daint/boost/1.56.0/gnu_482/lib:$LD_LIBRARY_PATH
# hdf5
export LD_LIBRARY_PATH=/users/biddisco/apps/daint/hdf5_1_8_cmake/lib:$LD_LIBRARY_PATH
# nvidia GL
export LD_LIBRARY_PATH=/opt/cray/nvidia/default/lib64:$LD_LIBRARY_PATH
# paraview
export LD_LIBRARY_PATH=/scratch/daint/biddisco/egpgv/:$LD_LIBRARY_PATH

# plugins
export LD_LIBRARY_PATH=/scratch/daint/biddisco/egpgv:$LD_LIBRARY_PATH
export PV_PLUGIN_PATH=/scratch/daint/biddisco/egpgv

export OMP_NUM_THREADS=1

aprun -n ${TASKS} -N ${NPERNODE} /scratch/daint/biddisco/egpgv/pvbatch --use-offscreen-rendering --disable-xdisplay-test /scratch/daint/biddisco/egpgv/scripts/benchmark-particle-partition-filter.py -g 0 -f $FILTER -p /scratch/daint/biddisco/data/sphflow/resampled/${FILE_PREFIX}_dambreak21.h5part > timing-log.txt

_EOF_

chmod 775 ${DIR_NAME}/slurm-exp.bash

echo "cd ${DIR_NAME}; sbatch slurm-exp.bash; cd \$BASEDIR" >> run_jobs.bash

}

# get the path to this generate script, works for most cases
pushd `dirname $0` > /dev/null
BASEDIR=`pwd`
popd > /dev/null
echo "Generating jobs using base directory $BASEDIR"

# Create another script to submit all generated jobs to the scheduler
echo "#!/bin/bash" > run_jobs.bash
echo "BASEDIR=$BASEDIR" >> run_jobs.bash
chmod 775 run_jobs.bash

# set some vars which are fixed in this test
QUEUE=viz
NODES=32
NPERNODE=1

# Loop through all the parameter combinations generating jobs for each
for NODES in 32 64 128 256 512 1024
do
  for NPERNODE in 1 2 4
  do
    for REDUCTION in 1 2 4 8 16 32 64 128 256 512 1024 2048 4096
    do
      for FILTER in 0 1
      do 
        write_script
      done
    done
  done
done

