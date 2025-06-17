#!/bin/bash
#SBATCH --job-name=single_cell_ref # Nom du job
#SBATCH --nodes=1 # 1 seul noeud
#SBATCH --ntasks=1 # 1 seule tâche
#SBATCH --mem=60gb # mémoire
#SBATCH --time=60:00:00 # Temps limite en hrs:min:sec
#SBATCH --account=dev_gpu # Utilisation du compte d’accès dev_gpu
#SBATCH --partition=batch_gpu  # Utilisation de la partition dev_gpu
#SBATCH --gres=gpu:1
#SBATCH --output=single_cell_ref_%j.out # Standard output
#SBATCH --error=single_cell_ref_%j.err # Standard error log

hostname
nvidia-smi
apptainer exec --nv /mnt/beegfs/common/containers/singularity/dev/spatialscope/spatialscope.sif \
/opt/conda/envs/SpatialScope/bin/python \
/mnt/beegfs/home/asobkow1/persistent/spatialscope/SpatialScope/src/Train_scRef.py \
--ckpt_path /mnt/beegfs/home/asobkow1/persistent/spatialscope/results/Ckpts_scRefs_epochs_1000_normalised \
--scRef /mnt/beegfs/home/asobkow1/persistent/data/Sc_data/Sc_cellxgene_normalised.h5ad \
--cell_class_column sample_cell_type \
--epoch 1000 \
--gpus 03
