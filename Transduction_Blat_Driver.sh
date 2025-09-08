#!/bin/bash
#SBATCH --mail-user=blacksmi@umich.edu
#SBATCH --mail-type=FAIL,ARRAY_TASKS
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=12G
#SBATCH --time=1:00:00
#SBATCH --account=jmkidd99
#SBATCH --partition=standard
#SBATCH --output=logs/%x-%A_%a.out.log
#SBATCH --export=ALL
#SBATCH --array=1-12
cd Revised_All_Canines_List_2025_08_06
run-by-id-log.pl Transduction_Detection/Blat_CMDs.txt ./logs/Transduction_Detection_Blat.log $SLURM_ARRAY_TASK_ID
