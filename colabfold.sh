#!/bin/sh
### General options
### –- specify queue --
#BSUB -q gpua100

### -- set the job Name --
#BSUB -J testjob

### -- ask for number of cores (default: 1) --
#BSUB -n 4

### -- Select the resources: 1 gpu in exclusive process mode --
#BSUB -gpu "num=1:mode=exclusive_process"

### -- set walltime limit: hh:mm --  maximum 24 hours for GPU-queues right now
#BSUB -W 24:00

# request system-memory per core
#BSUB -R "rusage[mem=4.2GB]"

### -- specify that the cores must be on the same host -- 
#BSUB -R "span[hosts=1]"

### -- set the email address --
# please uncomment the following line and put in your e-mail address,
# if you want to receive e-mail notifications on a non-default address
##BSUB -u s204116@student.dtu.dk

### -- send notification at start --
#BSUB -B

### -- send notification at completion--
#BSUB -N

### -- Specify the output and error file. %J is the job-id --
### -- -o and -e mean append, -oo and -eo mean overwrite --
#BSUB -o /work3/s204116/Projects/BO/error_messages/%J.out
#BSUB -e /work3/s204116/Projects/BO/error_messages/%J.err
# -- end of LSF options --

module load cuda/12.1

source /dtu/projects/RFdiffusion/setup.sh
module load colabfold

cd /work3/s204116/Projects/BO/af2pipeline

input_path="Inputs/bc_max.fasta"
output_colab_path="outputs/max_test_csntx"
output_metrics_path="outputs/max_test_csntx.csv"

colabfold_batch "${input_path}" "${output_colab_path}"
source /work3/s204116/anaconda3/bin/activate
conda activate BO
python extract_AF_metrics.py --colab_input "${input_path}" --colab_output_folder "${output_colab_path}" --output_csv "${output_metrics_path}"


