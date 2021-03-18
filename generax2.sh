#!/bin/bash
#SBATCH --job-name=generax2
#SBATCH --partition=compute
#SBATCH --time=3-0
#SBATCH --mem=150G
#SBATCH --mail-user=jigyasa.arora@oist.jp
#SBATCH --mail-type=BEGIN,FAIL,END

#SBATCH --output=generax2_%A-%a.out
#SBATCH -o out_generax2.%j
#SBATCH -e err_generax2.%j

##SBATCH --array 0-4
##num=$(printf "%04d" $SLURM_ARRAY_TASK_ID) #add two zeros infront

#load modules-

#files=(uniq-aligned-*fasta)
##echo "list: " ${files[${SLURM_ARRAY_TASK_ID}]} # this generates a list of $sf files
#file1=${files[${SLURM_ARRAY_TASK_ID}]} #it reads each index at a time
#file2=${file1/all.protein.aligned.fas/all.fas}

mpiexec -np 1 generax -f families.method2.txt -s hosttree-treponemacog0099_generax.nwk -p cog0099_method2 
