import os
import subprocess
import sys
import argparse

parser = argparse.ArgumentParser(prog="Produce_Miroepeats_Image.py",description="Creates an annotated Miropeats Image")
parser.add_argument("--genome1_path",dest="gen1_path",required=True,help="Path to fasta file containing sequence 1 (filled). Currently must come from a reference assembly.")
parser.add_argument("--genome2_path",dest = "gen2_path",required=True,help="Path to fasta file containing sequence 2 (empty). Currently must come from a reference assembly")
parser.add_argument("--genome1_coords",dest="gen1_coords",required=True,help="Genomic location in sample 1 (filled). Assumes 1 indexed")
parser.add_argument("--genome2_coords",dest="gen2_coords",required=True,help="Genomic location in sample 2 (empty). Assemes 1 indexed")
parser.add_argument("--reverse",dest="rev",required=False,default="False",help="If true, rev comp the empty site")
args = parser.parse_args()
   

#populate relevant variables into locus_info
locus_info = {}
locus_info["coords1"] = args.gen1_coords
locus_info["coords2"] = args.gen2_coords
locus_info["gen1_path"] = args.gen1_path
locus_info["gen2_path"] = args.gen2_path

#convert coords to lists for later
locus_info["coords1_list"] = args.gen1_coords.replace(":"," ").replace("-"," ").split()
locus_info["coords2_list"] = args.gen2_coords.replace(":"," ").replace("-"," ").split()
locus_info["flank_dist"] = max((int(locus_info['coords1_list'][2])-int(locus_info['coords1_list'][1])+1)*2,1000)
#create the coords which include flanking distances
relevant_keys = ["coords1_list","coords2_list"]
for key in relevant_keys:
    locus_info[key][1] = int(locus_info[key][1])
    locus_info[key][2] = int(locus_info[key][2])
    
    locus_info[f"{key}_with_flank"] = locus_info[key][0],locus_info[key][1]-locus_info["flank_dist"],locus_info[key][2]+locus_info["flank_dist"]
    locus_info[f"{key}_string"] = f'{locus_info[key][0]}:{locus_info[key][1]-locus_info["flank_dist"]}-{locus_info[key][2]+locus_info["flank_dist"]}'

print(locus_info)

extract1_prefix = f"{locus_info['coords1_list_string']}_extract1_miropeats.fa"
extract2_prefix = f"{locus_info['coords1_list_string']}_extract2_miropeats.fa"
if args.rev == "True" or args.rev == True:
    modifier = "-i"
else:
    modifier = ""
extract1 = f"samtools faidx {locus_info['gen1_path']} {locus_info['coords1_list_string']} > {extract1_prefix}"
extract2 = f"samtools faidx {modifier} {locus_info['gen2_path']} {locus_info['coords2_list_string']} > {extract2_prefix}"
subprocess.run(extract1,shell=True)
subprocess.run(extract2,shell=True)

rm1 = f"RepeatMasker -species dog {extract1_prefix} {extract2_prefix}"
subprocess.run(rm1,shell=True)

command = f'miropeats -onlyinter -s 200 -seq {extract2_prefix} {extract1_prefix} | tee > 200.mrout'
subprocess.run(command,shell=True)
#sys.exit()
subprocess.run('ps2pdf threshold200.ps', shell=True)

#sys.exit()
command = f'annotate-miropeats-2seqs.py --miroin threshold200.ps --siteID test \
--topRM {extract2_prefix}.out --bottomRM {extract1_prefix}.out --topName {locus_info["coords2_list_string"]} \
--bottomBreak {locus_info["flank_dist"]},{locus_info["flank_dist"]+int(locus_info["coords1_list"][2])-int(locus_info["coords1_list"][1])} --topBreak {locus_info["flank_dist"]},{locus_info["flank_dist"] +int(locus_info["coords2_list"][2])-int(locus_info["coords2_list"][1])}'
subprocess.run(command,shell=True)
command = f'ps2pdf -dEPSCrop threshold200.ps.annotated.ps Annotated_Image_{args.gen1_coords}_{args.gen2_coords}.pdf'
subprocess.run(command,shell=True)