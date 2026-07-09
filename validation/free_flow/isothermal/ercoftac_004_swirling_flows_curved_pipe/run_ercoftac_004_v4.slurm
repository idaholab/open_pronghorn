#!/bin/bash

# Submit this script with: sbatch thefilename

#SBATCH --time=0:30:00   # walltime

#SBATCH --ntasks-per-node=112   # number of processor cores (i.e. tasks)

#SBATCH --nodes=1   # number of nodes

#SBATCH --wckey=neams   # Project Code

#SBATCH --job-name="ercoftac_004_v4"

#SBATCH --mail-user=alberic.seydoux@inl.gov   # email address

#SBATCH --mail-type=BEGIN,END,FAIL

# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE

module load use.moose moose-dev-openmpi/2026.06.16
moose-dev-shell

cd /scratch/seydap/projects/open_pronghorn/validation/free_flow/isothermal/ercoftac_004_swirling_flows_curved_pipe

mpiexec -n 112 moose-dev-exec ../../../../open_pronghorn-opt \
    -i curvedpipe_noswirl_v4.i 2>&1 | tee curvedpipe_noswirl_v4.log
