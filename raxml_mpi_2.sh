#!/usr/bin/env bash
#SBATCH -J raxml-bs
#SBATCH -N 6
#SBATCH --ntasks-per-node=1
#SBATCH -c 96
#SBATCH -t 24:00:00
#SBATCH -o bs.%j.log

module load apptainer
export OMPI_MCA_btl=self,tcp,vader
export OMPI_MCA_btl_openib_allow_ib=0

apptainer exec \
  -B "$PWD:$PWD" \
  --pwd "$PWD" \
  /home/isaac12/sif/raxml.sif raxml-ng-mpi \
    --bootstrap \
    --msa T1_parse.raxml.rba \
    --prefix T1_bs_6nodes \
    --seed 123 \
    --threads 96 \
    --workers auto{96} \
    --bs-trees 100
