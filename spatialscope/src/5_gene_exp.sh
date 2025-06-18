#!/bin/bash
#SBATCH --job-name=gene_exp
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=200gb
#SBATCH --array=1-4  # Runs 4 jobs with different parameters
#SBATCH --time=100:00:00
#SBATCH --account=dev_gpu # Utilisation du compte d’accès dev_gpu
#SBATCH --partition=batch_gpu  # Utilisation de la partition dev_gpu
#SBATCH --gres=gpu:1
#SBATCH --output=gene_exp_%j.out
#SBATCH --error=gene_exp_%j.err

hostname
SLIDES=("18P06762" "17P04394" "18P02831" "18P03122")
# SLIDES=("17P02529" "18P08140" "18P06593")
SLIDE=${SLIDES[$SLURM_ARRAY_TASK_ID - 1]}
echo "Running sp data preprocessing job $SLURM_ARRAY_TASK_ID with parameter: $SLIDE"

H5AD_FILE="/mnt/beegfs/home/asobkow1/persistent/spatialscope/results/deconv_normalised/${SLIDE}/sp_adata.h5ad"
END_INDEX=$(apptainer exec --nv /mnt/beegfs/common/containers/singularity/dev/spatialscope/spatialscope.sif \
    /opt/conda/envs/SpatialScope/bin/python - <<END
import scanpy as sc
adata = sc.read_h5ad("$H5AD_FILE")
result = adata.obs.shape[0] - (adata.obs["cell_count"] == 0).sum()
print(result)
END
)

CHUNK_SIZE=200
spot_ranges=()
start=0
while [ $start -lt $END_INDEX ]; do
    end=$((start + CHUNK_SIZE))
    if [ $end -gt $END_INDEX ]; then
        end=$END_INDEX
    fi
    spot_ranges+=("${start},${end}")
    start=$end
done

echo "Found spot ranges:"
for r in "${spot_ranges[@]}"; do
    echo "  $r"
done

for range in "${spot_ranges[@]}"; do
    echo "Processing spot_range $range"
    apptainer exec --nv /mnt/beegfs/common/containers/singularity/dev/spatialscope/spatialscope.sif \
    /opt/conda/envs/SpatialScope/bin/python \
    /mnt/beegfs/home/asobkow1/persistent/spatialscope/SpatialScope/src/Decomposition.py \
    --tissue $SLIDE \
    --out_dir /mnt/beegfs/home/asobkow1/persistent/spatialscope/results/deconv_normalised \
    --SC_Data /mnt/beegfs/home/asobkow1/persistent/data/Sc_data/Sc_cellxgene_normalised.h5ad \
    --cell_class_column sample_cell_type \
    --ckpt_path /mnt/beegfs/home/asobkow1/persistent/spatialscope/results/Ckpts_scRefs_epochs_1000_normalised/model_01000.pt \
    --spot_range $range \
    --gpu 1,2,3
done
