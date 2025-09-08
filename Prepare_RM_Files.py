import os
import sys
import subprocess
import argparse

parser = argparse.ArgumentParser(prog="Performing_Alignments.py",description="Produces a list of minimap2 alignment commands.")
parser.add_argument("--sample_file",dest="sample_file",required=True,help="Formatted input file where each line contains (in a tab delimited list) a sample name, the reference path, a RepeatMasker path, the path to reference Gaps, and a path to segmental duplications in the reference")
parser.add_argument("--cmd_descriptor",dest = "cmd_descriptor",required=True,help="string that will be present in the commands file to prevent seperate command sets from overwriting each other")
parser.add_argument("--RM_path",dest="RM_path",required=True,help="THe directory that contains RepeatMasker files. Currently, this is a subdirectory of the directory the script is launched from")
parser.add_argument("--species_name",dest="species",default="dog",required=False,help="The species being processed. Currently only supports 'dog' and 'human'")
args = parser.parse_args()


#create RM_dict
RM_dict = {} 
with open(args.sample_file) as infile:
    for line in infile:
        if not line.startswith("#"):
            line = line.split()
            RM_dict[line[0]] = line[2]
print(RM_dict)

#Extract relevent portions from RM.bed files
current_dir = os.getcwd()
RM_path = f"{current_dir}/{args.RM_path}"
if not os.path.exists(RM_path):
    os.mkdir(RM_path)
os.chdir(RM_path)

cmd_file = f"{RM_path}/RM_Commands_{args.cmd_descriptor}.sh"
print(cmd_file)
additional_filters_dict = {}
if args.species.lower() == "dog":
    TEs = ["LINE/L1","SINEC"]
    additional_filters_dict["LINE/L1"] = "| grep -v 'HAL'"
    additional_filters_dict["SINEC"] = ""
elif args.species.lower() == "human":
    TEs = ["LINE/L1","SINE/Alu"]
    additional_filters_dict["LINE/L1"] = ""
    additional_filters_dict["SINE/Alu"] = ""
else:
    sys.exit("Invalid species name. Currently supported species are dog and human") #if this gets longer we can replace this step with a dictionary 
with open(cmd_file,'wt') as out:
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
#SBATCH --export=ALL\n\n"
    out.write(header)
    for TE in TEs:
        for key in RM_dict.keys():
            cmd = f'cat {RM_dict[key]} | grep "{TE}" {additional_filters_dict[TE]} > {RM_path}/{key}_{TE.split("/")[0]}_RM.out.bed\n'
            print(cmd)
            out.write(cmd)

if not os.path.exists(f"{RM_path}/logs"):
    os.mkdir(f"{RM_path}/logs")
subprocess.run(f"chmod ugo+x {cmd_file}", shell=True)
subprocess.run(f"sbatch {cmd_file}",shell=True)