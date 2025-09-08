import subprocess
import os
import argparse

parser = argparse.ArgumentParser(prog="Create_Detect_Hallmarks_CMDs.py",description="Creates commands to run Detect_Hallmarks_Revamp script in parallel. Use fullpaths to input files")
parser.add_argument("--sample_file",dest="sample_file",required=True,help="Formatted input file where each line contains (in a tab delimited list) a sample name, the reference path, a RepeatMasker path, the path to reference Gaps, and a path to segmental duplications in the reference")
parser.add_argument("--out_path_name",dest="aln_dir",required=True,help="The directory that the commands file will be placed in, this will be the same path that the alignments files are located in")
parser.add_argument("--RM_path",dest="RM_dir",required=True,help="The directory that contains RepeatMasker files. Currently, this is a subdirectory of the directory the script is launched from")
args = parser.parse_args()

os.chdir(args.aln_dir)
aln_dir_content = os.listdir(".")

#set variables
path_to_detect_hallmarks = "/nfs/turbo/umms-jmkidd/matt-projects/Former_KiddLabScratch_Data/inter-genome_comparisons/Aligning_Using_2.26/Detect_Hallmarks_Revamp.py"
count_of_jobs = 0


with open("Hallmarks_CMDs.txt",'wt') as outfile:
    for file_name in aln_dir_content:
        if "filled_intersect_with" not in file_name or not file_name.endswith(".bed") or "processed" in file_name:
            continue
        #elif "LINE" in file_name:#TODO temporary since we're not there yet.
        #    continue 
        cmd = f"python {path_to_detect_hallmarks} --sample_file {args.sample_file} --out_path_name {args.aln_dir} --RM_path {args.RM_dir} --file_name {file_name}"
        print(cmd)
        cmd = cmd + "\n"
        outfile.write(cmd)
        count_of_jobs+=1
        
header = f"#!/bin/bash\n\
#SBATCH --mail-user=blacksmi@umich.edu\n\
#SBATCH --mail-type=FAIL,ARRAY_TASKS\n\
#SBATCH --ntasks=1\n\
#SBATCH --cpus-per-task=1\n\
#SBATCH --ntasks-per-node=1\n\
#SBATCH --mem-per-cpu=8G\n\
#SBATCH --time=20:00:00\n\
#SBATCH --account=jmkidd99\n\
#SBATCH --partition=standard\n\
#SBATCH --output=logs/%x-%A_%a.out.log\n\
#SBATCH --export=ALL\n\
#SBATCH --array=1-{count_of_jobs}\n\n"

if not os.path.exists(f"logs"):
    os.mkdir(f"logs")
with open(f"Detect_Hallmarks_Launcher.sh",'wt') as out:
    out.write(header)
    out.write("date;\n")
    out.write(f"run-by-id-log.pl Hallmarks_CMDs.txt logs/Hallmarks_CMDs.log $SLURM_ARRAY_TASK_ID\n")
    out.write("date;\n")
subprocess.run(f"chmod ugo+x Detect_Hallmarks_Launcher.sh", shell=True)
subprocess.run(f"sbatch Detect_Hallmarks_Launcher.sh",shell=True)

