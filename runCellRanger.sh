#!/bin/bash
#SBATCH --job-name=cellranger_run
#SBATCH --output=CellRanger_run_array_%A_%a.out
#SBATCH --error=CellRanger_run_array_%A_%a.err
#SBATCH --time=46:00:00
#SBATCH --mem=80G
#SBATCH --cpus-per-task=16
#SBATCH --partition=amilan
#SBATCH --account=amc-general
#SBATCH --qos=long
#SBATCH --array=0-7%3  # Limit to 3 concurrent jobs
#SBATCH --mail-user=connor.eck@cuanschutz.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# Load modules
module load anaconda
conda activate cellranger_update_env

# Add Cell Ranger to PATH
export PATH=/scratch/alpine/ceck@xsede.org/Kuhn_Data/Sarah_Data/CellRanger_Update/cellranger-9.0.1:$PATH

# Confirm Cell Ranger is accessible
which cellranger
cellranger --version

# Sample names
SAMPLES=(
  863_spleen
  864_spleen  
  865_spleen
  863_epithelium 
  864_epithelium
  865_epithelium
  Pooled_JLN
  Pooled_spleen
)

# Get sample name for current array index
SAMPLE_NAME=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

# Set ID and CSV config
ID="${SAMPLE_NAME}_output"
CSV_FILE="/scratch/alpine/ceck@xsede.org/Kuhn_Data/Sarah_Data/CellRanger_Update/CONFIG_FILES/${SAMPLE_NAME}_config.csv"

# Create and move to isolated working directory
BASE_WORK_DIR="/scratch/alpine/ceck@xsede.org/Kuhn_Data/Sarah_Data/CellRanger_Update/CellRanger_All_Samples/SC5P_PE_V3_Chem"
WORK_DIR="${BASE_WORK_DIR}/${SAMPLE_NAME}"
mkdir -p "$WORK_DIR"
cd "$WORK_DIR"

# Ignore inode checks
export MRO_DISK_SPACE_CHECK=disable

# Run Cell Ranger multi
echo "Running Cell Ranger multi for $SAMPLE_NAME in $WORK_DIR"
cellranger multi \
  --id="$ID" \
  --csv="$CSV_FILE" \
  --localcores=16 \
  --localmem=80


