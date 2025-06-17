#!/bin/bash
#SBATCH --job-name=cell_iden # Nom du job
#SBATCH --nodes=1 # 1 seul noeud
#SBATCH --ntasks=1 # 1 seule tâche
#SBATCH --cpus-per-task=10
#SBATCH --mem=150gb # mémoire
#SBATCH --array=1-7  # Runs 7 jobs with different parameters
#SBATCH --time=10:00:00 # Temps limite en hrs:min:sec
#SBATCH --account=dev # Utilisation du compte d’accès dev
#SBATCH --partition=batch # Utilisation de la partition dev
#SBATCH --output=cell_iden_%j.out # Standard output
#SBATCH --error=cell_iden_%j.err # Standard error log

hostname
SLIDES=("17P02529" "18P06762" "18P08140" "17P04394" "18P06593" "18P03122" "18P02831")
SLIDE=${SLIDES[$SLURM_ARRAY_TASK_ID - 1]}
apptainer exec --nv /mnt/beegfs/common/containers/singularity/dev/spatialscope/spatialscope.sif \
/opt/conda/envs/SpatialScope/bin/python \
/mnt/beegfs/home/asobkow1/persistent/spatialscope/SpatialScope/src/Cell_Type_Identification.py \
--tissue $SLIDE \
--out_dir /mnt/beegfs/home/asobkow1/persistent/spatialscope/results/deconv_normalised \
--ST_Data /mnt/beegfs/home/asobkow1/persistent/spatialscope/results/deconv_normalised/$SLIDE/sp_adata_ns.h5ad \
--SC_Data /mnt/beegfs/home/asobkow1/persistent/data/Sc_data/Sc_cellxgene_normalised.h5ad \
--cell_class_column sample_cell_type
