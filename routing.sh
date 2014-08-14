#!/bin/bash
# Job name:
#SBATCH --jog-name=routing
#
# Partition:
#SBATCH --partition=lr1
#
# QoS:
#SBATCH --qos=lr_normal
#
# Account:
#SBATCH --account=ac_chenglyu
#
# Processors:
#SBATCH --ntasks=1
#
# Wall clock limit:
#SBATCH --time=24:0:00
#
# Mail type:
#SBATCH --mail-type=all
#
# Mail user:
#SBATCH --mail-user=chenglu@berkeley.edu

## Run command
module load matlab/R2011b
matlab -nodesktop -nodisplay < wedge_main.m &> output.out &
