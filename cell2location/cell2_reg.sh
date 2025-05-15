#!/bin/bash
#SBATCH --job-name=lusc_deconv # Nom du job
#SBATCH --nodes=1 # 1 seul noeud
#SBATCH --ntasks=1 # 1 seule tâche
#SBATCH --mem=60gb # mémoire
#SBATCH --time=02:00:00 # Temps limite en hrs:min:sec
#SBATCH --account=dev_gpu # Utilisation du compte d’accès dev_gpu
#SBATCH --partition=dev_gpu  # Utilisation de la partition dev_gpu
#SBATCH --gres=gpu:1
#SBATCH --output=lusc_deconv_%j.out # Standard output
#SBATCH --error=lusc_deconv_%j.err # Standard error log

hostname
apptainer exec --nv /mnt/beegfs/common/containers/singularity/dev/cell2location/cell2location.sif python cell2_reg.py

