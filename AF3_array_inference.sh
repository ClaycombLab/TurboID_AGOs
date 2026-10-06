#!/bin/bash
#SBATCH --job-name=af3-infer
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH --gpus=h100
#SBATCH --output=slurm-%A_%a.out
#SBATCH --error=slurm-%A_%a.err
#SBATCH --array=0-65%5

module load StdEnv/2023 hmmer/3.4 rdkit/2024.03.5 python/3.12 cuda/12.2 cudnn/9.2

AF3_ENV=$HOME/af3-venv
AF3_REPO=$HOME/alphafold3
MODEL_DIR=$HOME/models
CPU_STAGE_DIR=$SCRATCH/alphafold/AGO_RNA_output
FINAL_OUT=$SCRATCH/alphafold/AGO_RNA_inference
MANIFEST=$CPU_STAGE_DIR/data_manifest.txt

# This is an adapted version of the inference script from DRAC
# The main modification is the use of manifest for batch array
# And the use of a fixed python env

# Build manifest
if [[ ! -f $MANIFEST ]]; then
  find $CPU_STAGE_DIR -type f -name '*_data.json' | sort > $MANIFEST
fi

mapfile -t DATA_JSONS < $MANIFEST
N=${#DATA_JSONS[@]}


DATA_JSON=${DATA_JSONS[$SLURM_ARRAY_TASK_ID]}

# Load the environment and varaibles acording to DRAC
source $AF3_ENV/bin/activate
export XLA_FLAGS="--xla_gpu_enable_triton_gemm=false"
export XLA_PYTHON_CLIENT_PREALLOCATE=true
export XLA_CLIENT_MEM_FRACTION=0.95
export TMPDIR="$SLURM_TMPDIR"
export TEMP="$SLURM_TMPDIR"
export TMP="$SLURM_TMPDIR"
export JAX_COMPILATION_CACHE_DIR="$SLURM_TMPDIR/jax_cache"

mkdir -p $JAX_COMPILATION_CACHE_DIR $FINAL_OUT

BASENAME=$(basename "$DATA_JSON" _data.json)
TASK_OUT=$SLURM_TMPDIR/out_${BASENAME}


python run_alphafold.py \
  --json_path=$DATA_JSON \
  --model_dir=$MODEL_DIR \
  --output_dir=$TASK_OUT \
  --jax_compilation_cache_dir=$JAX_COMPILATION_CACHE_DIR \
  --norun_data_pipeline \
  --run_inference


mkdir -p $FINAL_OUT/$BASENAME

if command -v rsync >/dev/null 2>&1; then
  rsync -a $TASK_OUT/ $FINAL_OUT/$BASENAME/
else
  cp -a $TASK_OUT/. $FINAL_OUT/$BASENAME/
fi

