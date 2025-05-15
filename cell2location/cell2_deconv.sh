#!/bin/bash
#SBATCH --job-name=lusc_deconv # Nom du job
#SBATCH --nodes=1 # 1 seul noeud
#SBATCH --ntasks=1 # 1 seule tâche
#SBATCH --mem=60gb # mémoire
#SBATCH --array=1-7  # Runs 7 jobs with different parameters
#SBATCH --time=03:00:00 # Temps limite en hrs:min:sec
#SBATCH --account=dev_gpu # Utilisation du compte d’accès dev_gpu
#SBATCH --partition=dev_gpu  # Utilisation de la partition dev_gpu
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:1
#SBATCH --output=/mnt/beegfs/home/asobkow1/persistent/logs/lusc_deconv_%J_%j.out # Standard output
#SBATCH --error=/mnt/beegfs/home/asobkow1/persistent/logs/lusc_deconv_%J_%j.err # Standard error log

hostname
PARAMS=("17P02529" "18P06762" "18P08140" "17P04394" "18P06593" "18P03122" "18P02831")
PARAM=${PARAMS[$SLURM_ARRAY_TASK_ID - 1]}
echo "Running job $SLURM_ARRAY_TASK_ID on GPU with parameter: $PARAM"
apptainer exec --nv /mnt/beegfs/common/containers/singularity/dev/cell2location/cell2location.sif python cell2_deconv.py --param $PARAM
