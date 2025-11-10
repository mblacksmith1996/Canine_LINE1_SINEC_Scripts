import os
import subprocess
import sys
import argparse
import shutil

parser = argparse.ArgumentParser(prog="Process_Locus_for_Diagram.py",description="Facilitates Easier Processing of Comparing a Filled and Empty Site. Requires manual curation")
parser.add_argument("--Processed_File",dest="Processed",required=True,help="Path to processed file containing information about the locus being investigated")
parser.add_argument("--Locus_File",dest="Loci",required=True,help="File containing loci to be processed")
parser.add_argument("--outdir",dest="outdir",required=False,default=".",help="Directory to place output files")
args = parser.parse_args()
#This script has been designed to undergo manual curation. There are probably edge cases that will not be good candidates for this script such as loci with massive TSDs


def check_for_modules():
    print("Get module version information")
    required_modules = ['python','samtools','RepeatMasker','bedtools','miropeats']
    
    for module in required_modules:
        if shutil.which(module) == "None":
            sys.exit(f"Missing {module}")
        print(module,shutil.which(module))
            
def extract_from_Locus_file(input_file):
    infile = open(input_file,'rt') 
    list_of_variants = []
    for line in infile.readlines():
        line = line.split()
        list_of_variants.append([line[0],line[1]])
    print(list_of_variants)
    infile.close()
    return list_of_variants
    
def extract_information_from_Processed_File(input_file,list_of_variants):
    infile = open(input_file,'rt')
    list_of_info = []
    for line in infile.readlines():
        line = line.split()
        if line[0] == "#chr":
            header = line
        else:
            if [line[0],line[1]] in list_of_variants:
                dict_of_line = {}
                #print(list_of_variants)
                for i in range(len(line)):
                    dict_of_line[header[i]] = line[i]
                list_of_info.append(dict_of_line)
                #print(dict_of_line)
                
    if len(list_of_info) != len(list_of_variants):
        print("Num found in processed file", len(list_of_info), "\nNum variants to process", len(list_of_variants))
        sys.exit("At least one variant not found in the processed file!!!!")
    infile.close()
    
    return list_of_info
    
def convert_fasta_to_all_caps(input_name):
    infile = open(input_name,'rt')
    outfile = open(f"{input_name}.upper",'wt')
    
    for line in infile.readlines():
        if line.startswith(">"):
            outfile.write(line)
        else:
            outfile.write(line.upper().rstrip())
    outfile.write("\n")
    
    infile.close()
    outfile.close()
    
def mask_LINE1_segments(input_file):
    if input_file.endswith("upper"):
        input_file = input_file[0:-6]
        
    RM_file = open(f"{input_file}.out",'rt')
    list_of_mask_coords = []
    for line in RM_file.readlines()[3::]:
        #print(line)
        line = line.split()
        if not line[10] == "LINE/L1":
            continue
        list_of_mask_coords.append([line[4],str(int(line[5])-1),line[6]])
    #print(list_of_mask_coords)
    if len(list_of_mask_coords) == 0:
        sys.exit("No LINE-1 segments. Shouldn't be possible")
    
    outname = f"{input_file}_coords_for_softmask"
    outfile = open(outname,'wt')
    for coords in list_of_mask_coords:
        out_line = "\t".join(coords) + "\n"
        outfile.write(out_line)
    outfile.close()
    
    cmd = f"bedtools maskfasta -fi {input_file}.upper -bed {outname} -fo {input_file}.softmasked -soft"
    subprocess.run(cmd,shell=True)
    #process the fasta.out from RM
  
def print_file_contents(infile):
    x = open(infile,'rt')
    for line in x.readlines():
        print(line.rstrip())
    x.close()
    #print('here')
#check for modules
check_for_modules()
#sys.exit()
#process the variants
list_of_variants = extract_from_Locus_file(args.Loci)

list_of_dicts = extract_information_from_Processed_File(args.Processed,list_of_variants)

for dictionary in list_of_dicts:
    #print(dictionary)
    path_to_filled = dictionary['Ref']
    path_to_empty = dictionary['Query']
    
    #define how many extra bp to extract
    #maybe start with TSD + 15 if there is a TSD of 10+ 
    #otherwise 25
    
    extra_dist = 15
    if len(dictionary['Age_TSD_Seq_right']) >= 10:
        filled_end = int(dictionary['right_TSD_filled'].split('\\')[-1])
        filled_start = int(dictionary['left_TSD_filled'].split('\\')[0])
        
        filled_end_with_extra = filled_end + extra_dist
        filled_start_with_extra = filled_start - extra_dist
        print(filled_start,filled_end)
        
        mod = len(dictionary['Age_TSD_Seq_right']) + extra_dist -1
        #sys.exit()
    else:
        filled_start = int(dictionary['Refined_Filled_left_bound'])
        filled_end = int(dictionary['Refined_Filled_right_bound'])
        
        filled_end_with_extra = filled_end + extra_dist +9
        filled_start_with_extra = filled_start - extra_dist -9
        mod = extra_dist +9
    empty_start = min(int(dictionary['Refined_Empty_left_bound']),int(dictionary['Refined_Empty_right_bound'])) 
    empty_end = max(int(dictionary['Refined_Empty_left_bound']),int(dictionary['Refined_Empty_right_bound']))
    
    empty_start_with_mod = empty_start - mod
    empty_end_with_mod = empty_end + mod
    
    print(dictionary['#chr'],filled_start_with_extra,filled_end_with_extra,empty_start_with_mod,empty_end_with_mod)
    
    
    file_prefix = f"{dictionary['#chr']}-{dictionary['Filled_start']}"
    
    
    if int(dictionary['Refined_Empty_left_bound']) < int(dictionary['Refined_Empty_right_bound']):
        mod = ""
        reverse_status = "False"
    else:
        mod = "-i"
        reverse_status = "True"
    #extract with samtools
    extract_cmd = f"samtools faidx {path_to_filled} {dictionary['#chr']}:{filled_start_with_extra}-{filled_end_with_extra} > {args.outdir}/{file_prefix}_filled_extract.fa"
    subprocess.run(extract_cmd,shell=True)
    
    extract_cmd = f"samtools faidx {mod} {path_to_empty} {dictionary['#chr']}:{empty_start_with_mod}-{empty_end_with_mod} > {args.outdir}/{file_prefix}_empty_extract.fa"
    subprocess.run(extract_cmd,shell=True)
    
    #use RM to get the mask of the filled site
    RM_cmd = f"RepeatMasker --species dog {args.outdir}/{file_prefix}_filled_extract.fa"
    subprocess.run(RM_cmd,shell=True)
    
    
    #Next, lets try to mask only the LINE-1 segment rather than all of the sequence.
    #To do so I will convert the entire sequence to all caps.
    convert_fasta_to_all_caps(f"{args.outdir}/{file_prefix}_filled_extract.fa")
    convert_fasta_to_all_caps(f"{args.outdir}/{file_prefix}_empty_extract.fa")
    
    #mask TE segments 
    mask_LINE1_segments(f"{args.outdir}/{file_prefix}_filled_extract.fa.upper")
    
    #prep miropeats command
    miro_path = "/nfs/turbo/umms-jmkidd/matt-projects/Former_KiddLabScratch_Data/inter-genome_comparisons/Aligning_Using_2.26/Produce_Miropeats_Image.py"
    cmd = f"python {miro_path} --genome1_path {path_to_filled} --genome2_path {path_to_empty} --genome1_coords {dictionary['#chr']}:{filled_start}-{filled_end} --genome2_coords {dictionary['#chr']}:{empty_start}-{empty_end} --reverse {reverse_status}"
    subprocess.run(cmd,shell=True)


    print("Prepping Variant Final Output")
    for key in dictionary.keys():
        if key == "poly_a_list":
            continue
        print(key, dictionary[key])
    print_file_contents(f"{args.outdir}/{file_prefix}_filled_extract.fa.softmasked")
    print_file_contents(f"{args.outdir}/{file_prefix}_empty_extract.fa.upper")
    print_file_contents(f"{args.outdir}/{file_prefix}_filled_extract.fa.out")

    
    