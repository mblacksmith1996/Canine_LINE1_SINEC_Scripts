import os
import subprocess


directory = "/nfs/turbo/umms-jmkidd/matt-projects/Former_KiddLabScratch_Data/inter-genome_comparisons/Aligning_Using_2.26/Revised_All_Canines_List_2025_08_06"

processed_filtered_files = []
for item in os.listdir(directory):
    if item.endswith("processed_filtered.bed"):
        processed_filtered_files.append(item)
print(processed_filtered_files,len(processed_filtered_files))

cmds = []
for file in processed_filtered_files:
    cmd = f"python scale_up_testing_poly_a_all_genomes.py --parse_file {directory}/{file} --output_dir {directory}/Relaxed_Calling\n"
    print(cmd)
    cmds.append(cmd)
with open(f"{directory}/Relaxed_Calling/poly_a_cmds.txt",'wt') as outfile:
    for cmd in cmds:
        outfile.write(cmd)

#make driver file
job_count = len(cmds)
header = f"#!/bin/bash\n\
#SBATCH --mail-user=blacksmi@umich.edu\n\
#SBATCH --mail-type=FAIL,ARRAY_TASKS\n\
#SBATCH --ntasks=1\n\
#SBATCH --cpus-per-task=1\n\
#SBATCH --ntasks-per-node=1\n\
#SBATCH --mem-per-cpu=8G\n\
#SBATCH --time=2:00:00\n\
#SBATCH --account=jmkidd99\n\
#SBATCH --partition=standard\n\
#SBATCH --output=logs/%x-%A_%a.out.log\n\
#SBATCH --export=ALL\n\
#SBATCH --array=1-{job_count}\n\n"

if not os.path.exists(f"{directory}/Relaxed_Calling/logs"):
    os.mkdir(f"{directory}/Relaxed_Calling/logs")
with open(f"{directory}/Relaxed_Calling/poly_a_cmds.sh",'wt') as out:
    out.write(header)
    out.write(f"run-by-id-log.pl {directory}/Relaxed_Calling/poly_a_cmds.txt {directory}/Relaxed_Calling/logs/poly_a_cmds.log $SLURM_ARRAY_TASK_ID")
subprocess.run(f"chmod ugo+x {directory}/Relaxed_Calling/poly_a_cmds.sh", shell=True)
subprocess.run(f"sbatch {directory}/Relaxed_Calling/poly_a_cmds.sh",shell=True)
