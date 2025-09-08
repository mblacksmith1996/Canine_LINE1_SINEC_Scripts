#!/bin/bash
#SBATCH --job-name=Detect_Hallmarks
#SBATCH --mail-user=blacksmi@umich.edu
#SBATCH --mail-type=FAIL,ARRAY_TASKS
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1  
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=6G 
#SBATCH --time=40:00:00
#SBATCH --account=jmkidd1
#SBATCH --partition=standard
#SBATCH --output=logs/%x-%A_%a.out.log
#SBATCH --export=ALL

python Detect_Hallmarks_Revamp.py --sample_file ./Prep_Files/Canine_Comparison_Revised_Canine_List.txt --out_path_name Revised_All_Canines_List_2025_03_25/ --RM_path RM_info_revised/
