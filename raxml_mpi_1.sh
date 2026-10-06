#!/usr/bin/env bash
#SBATCH -J raxml-search
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH -c 96
#SBATCH -t 12:00:00
#SBATCH -o search.%j.log

set -euo pipefail

apptainer exec \
  -B "$PWD:$PWD" \
  --pwd "$PWD" \
  /home/isaac12/sif/raxml.sif raxml-ng \
    --search \
    --msa T1_parse.raxml.rba \
    --prefix T1_search \
    --seed 123 \
    --threads auto{96} \
    --workers auto{10} \
    --tree pars{5},rand{5}
