#!/bin/bash
#SBATCH --job-name=af3-data
#SBATCH --time=72:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=220G
#SBATCH --array=0-65%10
#SBATCH --output=slurm-%A_%a.out
#SBATCH --error=slurm-%A_%a.err
module load StdEnv/2023 hmmer/3.4 rdkit/2024.03.5 python/3.12

DOWNLOAD_DIR=$SCRATCH/alphafold/dbs
INPUT_DIR=$SCRATCH/alphafold/AGO_RNA_input
FINAL_OUT=$SCRATCH/alphafold/AGO_RNA_output
AF3_ENV=$HOME/af3-venv
AF3_REPO=$HOME/alphafold3

# Activate the AlphaFold 3 environment
source $AF3_ENV/bin/activate

# Set JAX flags to disable Triton GEMM and specify the compilation cache directory
export XLA_FLAGS="--xla_gpu_enable_triton_gemm=false"
export JAX_COMPILATION_CACHE_DIR=$SLURM_TMPDIR/jax_cache
mkdir -p $FINAL_OUT $JAX_COMPILATION_CACHE_DIR

# Find all JSON files in the input directory and store them in an array
mapfile -t JSONS < <(find $INPUT_DIR -maxdepth 1 -type f -name "*.json" | sort)

JSON_FILE=${JSONS[$SLURM_ARRAY_TASK_ID]}
BASENAME=$(basename "$JSON_FILE" .json)
TASK_OUT=$SLURM_TMPDIR/out_${BASENAME}

# Run Alphafold3 without inference and just build the MSA and template features
python run_alphafold.py \
  --db_dir=$DOWNLOAD_DIR \
  --json_path=$JSON_FILE \
  --output_dir=$TASK_OUT \
  --jax_compilation_cache_dir=$JAX_COMPILATION_CACHE_DIR \
  --nhmmer_n_cpu=2 \
  --jackhmmer_n_cpu=2 \
  --norun_inference

mkdir -p $FINAL_OUT/$BASENAME
cp -vr $TASK_OUT/. $FINAL_OUT/$BASENAME
