import os
import subprocess
import argparse
import sys

def get_paftools_path():
    minimap2_path = subprocess.Popen("which minimap2",shell=True,stdout=subprocess.PIPE)
    minimap2_path,err = minimap2_path.communicate()
    minimap2_path = minimap2_path.decode()
    paf_path = minimap2_path[0:-9] + "paftools.js"
    return paf_path
def check_path_validity(input_file_list):
    print(input_file_list)
    for file in input_file_list:
        if file == "N/A":
            continue
        else:
            if not os.path.isfile(file):
                sys.exit(f"ERROR: {file} does not exist!!")

parser = argparse.ArgumentParser(prog="Performing_Alignments.py",description="Produces a list of minimap2 alignment commands.")
parser.add_argument("--sample_file",dest="sample_file",required=True,help="Formatted input file where each line contains (in a tab delimited list) a sample name, the reference path, a RepeatMasker path, the path to reference Gaps, and a path to segmental duplications in the reference")
parser.add_argument("--reference_sample",dest = "aln_ref",required=True,help="Sample to which all others will be aligned. Must be present in the sample file")
parser.add_argument("--out_path_name",dest="out_dir",required=True,help="The directory that contents will be placed in. Currently, this is a subdirectory of the directory the script is launched from")
args = parser.parse_args()

#create the genome_dict
genome_dict = {}
with open(args.sample_file,'rt') as infile:
    for line in infile:
        print(line)
        if line.startswith("#"):
            continue
        split_line = line.split()
        check_path_validity(split_line[1::])
        genome_dict[split_line[0]] = split_line[1]
print(genome_dict)

paftools_path = get_paftools_path()
print(paftools_path)

current_dir = os.getcwd()
aln_path = f"{current_dir}/{args.out_dir}"
if not os.path.exists(aln_path):
    os.mkdir(aln_path)
os.chdir(aln_path)

with open(f"{aln_path}/Aln_Commands.txt",'wt') as out:
    #Generate Minimap2 commands
    for genome in genome_dict.keys():
        if genome == args.aln_ref:
            continue
        cmd1 = f"minimap2 -c -x asm5 --cs {genome_dict[args.aln_ref]} {genome_dict[genome]} | sort -k6,6 -k8,8n > {genome}_against_{args.aln_ref}.paf && k8 {paftools_path} call {genome}_against_{args.aln_ref}.paf > {genome}_against_{args.aln_ref}.paf-vars \n" # && {paftools_path} call -f {genome_dict[args.aln_ref]} -s {args.aln_ref} {genome}_against_{args.aln_ref}.paf > {genome}_against_{args.aln_ref}.vcf \n"
        print(cmd1)
        out.write(cmd1)
job_count = len(genome_dict.keys())-1
header = f"#!/bin/bash\n\
#SBATCH --mail-user=blacksmi@umich.edu\n\
#SBATCH --mail-type=FAIL,ARRAY_TASKS\n\
#SBATCH --ntasks=1\n\
#SBATCH --cpus-per-task=3\n\
#SBATCH --ntasks-per-node=1\n\
#SBATCH --mem-per-cpu=8G\n\
#SBATCH --time=10:00:00\n\
#SBATCH --account=jmkidd1\n\
#SBATCH --partition=standard\n\
#SBATCH --output=logs/%x-%A_%a.out.log\n\
#SBATCH --export=ALL\n\
#SBATCH --array=1-{job_count}\n\n"

if not os.path.exists(f"logs"):
    os.mkdir(f"logs")
with open(f"Aln_Commands_Driver.sh",'wt') as out:
    out.write(header)
    out.write(f"run-by-id-log.pl Aln_Commands.txt logs/Aln_Commands.log $SLURM_ARRAY_TASK_ID")
subprocess.run(f"chmod ugo+x Aln_Commands_Driver.sh", shell=True)
subprocess.run(f"sbatch Aln_Commands_Driver.sh",shell=True)

