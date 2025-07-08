#!/bin/bash
#SBATCH --job-name=hypomap_hires_cibersortx_job
#SBATCH --nodes=2
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=32
#SBATCH --mem=200G
#SBATCH --time=48:00:00
#SBATCH --gpus=a100:2
#SBATCH --gpus-per-node a100:1
#SBATCH --constraint=interconnect_hdr,chip_manufacturer_intel


apptainer exec -B /project/jgeorg4/georgelab/Prakrit/Deconvolution_Analysis_December_2024/input:/src/data \
  -B /home/psubba/Prakrit_georgelab_shortcut/Deconvolution_Analysis_December_2024/output/Hypomap_output:/src/outdir \
  /project/jgeorg4/georgelab/Prakrit/Deconvolution_Analysis_December_2024/hires_latest.sif \
  /src/CIBERSORTxHiRes --username psubba@clemson.edu --token 6ac6a4fa284c94e1485fa6439f932f8a \
  --mixture /src/data/CIBERSORTx_Mixtures_Adjusted.txt \
  --sigmatrix /src/data/Hypomap_custom_signature_matrix.txt \
  --QN FALSE \
  
  