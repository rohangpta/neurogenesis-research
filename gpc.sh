#!/bin/bash

#SBATCH -J simulationDataGet
#SBATCH -D ./
#SBATCH --mail-user=grohan@seas.upenn.edu
#SBATCH --mail-type=ALL
#SBATCH -o /home/grohan/out/simulationCode."%j".out
#SBATCH -e /home/grohan/out/simulationCode."%j".err
#SBATCH --time=200:0:0
#SBATCH --partition=highmem
#SBATCH --mem=24gb

. /etc/profile.d/modules.sh
module load matlab

matlab -nodisplay -nodesktop -nosplash -r "simulationCode;exit"