#!/bin/bash
#SBATCH -A secretome
#SBATCH --mem 4G
#SBATCH -o slurm.%j.out
#SBATCH --array=0-6      ## we have 6 fasta of 300-seqs in each set


module load python/3.7
conda create -y -n outcyte
conda activate outcyte


executable_dir="/projects/secretome/Local_Install_OutCyte" #> where run_outcyte.py script is called from

workdir="/projects/secretome/"
cd ${workdir}

fasta_dir=${workdir}/Protein_Sequences/my_gene_set
## "PepSeq_[....].fasta"

results_dir=${workdir}/Results_Outcyte/my_gene_set
mkdir -p ${results_dir}

## TO RUN:python run_outcyte.py /path/to/your/fasta/file /outcyte-sp/or/outcyte-ups/or/outcyte

##A/ when fasta contain only 100 seq:
pepseqs=(${fasta_dir}/*.fasta)
python ${executable_dir}/run_outcyte_edMK.py ${pepseqs[$SLURM_ARRAY_TASK_ID]} outcyte ${results_dir}


## results saved in "results/" folder, by default
## Note: run_outcyte.py --> modified to _edMK, by adding a custom results_folder in user-input-arguments
