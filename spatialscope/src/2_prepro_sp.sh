#!/bin/bash
#SBATCH --job-name=sp_prepro # Nom du job
#SBATCH --nodes=1 # 1 seul noeud
#SBATCH --ntasks=1 # 1 seule tâche
#SBATCH --mem=60gb # mémoire
#SBATCH --array=1-1  # Runs 7 jobs with different parameters
#SBATCH --time=04:00:00 # Temps limite en hrs:min:sec
#SBATCH --account=dev # Utilisation du compte d’accès dev
#SBATCH --partition=batch  # Utilisation de la partition batch
#SBATCH --output=sp_prepro_%j.out # Standard output
#SBATCH --error=sp_prepro_%j.err # Standard error log

hostname
# SLIDES=("17P02529" "18P06762" "18P08140" "17P04394" "18P06593" "18P03122" "18P02831")
SLIDE=("18P06593")
# SLIDE=${SLIDES[$SLURM_ARRAY_TASK_ID - 1]}
echo "Running sp data preprocessing job $SLURM_ARRAY_TASK_ID with parameter: $SLIDE"
apptainer exec /mnt/beegfs/common/containers/singularity/dev/spatialscope/spatialscope.sif \
/opt/conda/envs/SpatialScope/bin/python \
/mnt/beegfs/home/asobkow1/persistent/spatialscope/2_prepro_sp.py \
--Slide_Path /mnt/beegfs/home/asobkow1/persistent/data/LUSC_v2/$SLIDE \
--Analysis_level 2
