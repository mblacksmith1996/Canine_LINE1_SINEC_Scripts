import sys
import os
import subprocess
#sys.path.append("/home/blacksmi/links/kidd-lab-scratch/matt-projects/CLAINE_pull/CLAINE/")
sys.path.append("/home/blacksmi/links/kidd-lab-umms/matt-projects/Hallmark_Detection_python_Scripts_Dev/Hallmark_Detection_Python_Scripts")
import Detecting_Hallmarks_Functions as Hallmarks
import argparse

#parse input args
parser = argparse.ArgumentParser(prog="Detect_Hallmarks_Revamp.py",description="Detects Target site duplications, Poly(A) tails and .")
parser.add_argument("--sample_file",dest="sample_file",required=True,help="Formatted input file where each line contains (in a tab delimited list) a sample name, the reference path, a RepeatMasker path, the path to reference Gaps, and a path to segmental duplications in the reference")
parser.add_argument("--out_path_name",dest="aln_dir",required=True,help="The directory that contents will be placed in. Currently, this is a subdirectory of the directory the script is launched from")
parser.add_argument("--RM_path",dest="RM_dir",required=True,help="The directory that contains RepeatMasker files. Currently, this is a subdirectory of the directory the script is launched from")
parser.add_argument("--file_name",dest="file_name",required=False,help="An optional file path if it desired that the script only run on a single file")
args = parser.parse_args()
#Todo. Within the post processing utils get rid of all references to mCan and Mischka directly
#create dictionary of reference sequences
genome_dict = {}
with open(args.sample_file) as infile:
    for line in infile.readlines()[1::]:
        line = line.split()
        genome_dict[line[0]] = line[1]
os.chdir(args.aln_dir)


print(os.getcwd())


required_keys = [] #this will be used to ensure that every locus has the same keys.
header = [] #this is used to dynamically create the header. 

if args.file_name != None:
    file_names = [args.file_name]
else:
    file_names = []
    dir_contents = os.listdir(".")
    for file_name in dir_contents:
        if "filled_intersect_with" not in file_name or not file_name.endswith(".bed") or "processed" in file_name:
            continue
        else:
            file_names.append(file_name)
print("File(s) to be processed, as well as number of files being processed.",file_names,len(file_names))
for file_name in file_names:
    if "filled_intersect_with" not in file_name or not file_name.endswith(".bed") or "processed" in file_name:
        continue
    split_file_name = file_name.split("_")
    if split_file_name[0] not in genome_dict.keys():
        continue
    
    if "LINE" in file_name:
        #continue #ToDo. Remove once done with SINE.
        TE = "LINE"
        TE_class = "LINE/L1"
        RM_modifier = "LINE"
    elif "SINE" in file_name:
        #continue #ToDo. Remove once done with LINE
        TE = "SINE"
        TE_class = "SINE/tRNA"
        RM_modifier = "SINEC"
    else:
        sys.exit("File name does not contain name of a valid TE")
    print(file_name)
    
    
    split_name = file_name.split("/")[-1].split("_")
    if "ref" in file_name: #present in reference absent in query
        status = "ref"
        ref = genome_dict[split_name[1]]
        query = genome_dict[split_name[0]]
        RM_name = f"{args.RM_dir}/{split_name[1]}_{RM_modifier}_RM.out.bed"
    elif "query" in file_name: #present in query absent in reference
        status = "query"
        ref = genome_dict[split_name[0]]
        query = genome_dict[split_name[1]]
        RM_name = f"{args.RM_dir}/{split_name[0]}_{RM_modifier}_RM.out.bed"
    else:   
        sys.exit("Improper File Name")
        
    #create dict of chrom lengths
    chr_dir = {}
    with open(f"{ref}.fai",'rt') as infile:
        for line in infile:
            line = line.split()
            chr_dir[line[0]] = int(line[1])
    processed_file = open(f"{file_name[0:-4]}_processed.bed", "wt")
    if header != []:
        processed_file.write(header)

    #prep tmp rm_files
    with open(file_name,'rt') as infile:
        if not os.path.exists(f"./{split_name[0]}_{status}_{TE}"):
            os.mkdir(f"./{split_name[0]}_{status}_{TE}")
        os.chdir(f"./{split_name[0]}_{status}_{TE}")
        
        #prep tmp rm_files
        for i in range(1,39):
            cmd = f'grep -w "^chr{i}" {RM_name} > chr{i}.tmp'
            subprocess.run(cmd,shell=True)
        cmd = f'grep -w "^chrX" {RM_name} > chrX.tmp'
        subprocess.run(cmd,shell=True)
        
        seperate_loci = open("seperate_loci.txt",'wt')
        #z = 0
        for line in infile.readlines():#[0:10]:
            line = line.split()
            locus_info = {}
            locus_info["chr"] = line[0]
            locus_info["Filled_start"] = int(line[1]) # These intro coordinates are 0 based but everything else is 1 based.
            locus_info["Filled_end"] = int(line[2])# These intro coordinates are 0 based but everything else is 1 based.
            locus_info["Empty_start"] = int(line[4])# These intro coordinates are 0 based but everything else is 1 based.
            locus_info["Empty_end"] = int(line[5])# These intro coordinates are 0 based but everything else is 1 based.
            locus_info["Insertion_sequence"] = line[6]
            locus_info["Reference_orientation"] = line[7]
            locus_info["Ref"] = ref #ref refers to the Filled site
            locus_info["Query"] = query #query refers to the Empty site
            #locus_info["chrom_length"] = chr_dir[locus_info["chr"]]
            locus_info["RM_name"] = RM_name
            locus_info["RM_Class"] = TE_class
            locus_info["Interior_TE_dist"] = 50
            locus_info["Prefix"] = file_name
            if TE == "SINE":
                locus_info["age_flank"] = 1000
            elif TE == "LINE":
                locus_info["age_flank"] = 1000
            print(locus_info)

            #get modifier for age command
            if locus_info["Reference_orientation"] == "-":
                locus_info["Age_orientation"] = "-revcom2"
            elif locus_info["Reference_orientation"] == "+":
                locus_info["Age_orientation"] = ""
            else:
                sys.exit("Invalid reference orientation")
            
            #do an age step here
            loop_break = True
            while loop_break:
                extract_file_prefix = f"{locus_info['chr']}-{locus_info['Filled_start']}-{locus_info['Filled_end']}"
                locus_info['age_flank'] = min(locus_info['age_flank'],locus_info['Filled_start']-1,locus_info['Empty_start']-1)
                print(locus_info)
                filled_extract = f"samtools faidx {locus_info['Ref']} {locus_info['chr']}:{max(locus_info['Filled_start']-locus_info['age_flank']+1,1)}-{locus_info['Filled_end']+locus_info['age_flank']} > {extract_file_prefix}_filled.fa"
                empty_extract = f"samtools faidx {locus_info['Query']} {locus_info['chr']}:{max(locus_info['Empty_start']-locus_info['age_flank']+1,1)}-{locus_info['Empty_end']+locus_info['age_flank']} > {extract_file_prefix}_empty.fa" 
                print(filled_extract)
                print(empty_extract)
                subprocess.run(filled_extract,shell=True)
                subprocess.run(empty_extract,shell=True)
                age_align_cmd = f"age_align {extract_file_prefix}_filled.fa {extract_file_prefix}_empty.fa {locus_info['Age_orientation']} > {extract_file_prefix}_age.txt"
                print(age_align_cmd)
                subprocess.run(age_align_cmd,shell=True)
                age_out = Hallmarks.get_info_from_subprocess(f"cat {extract_file_prefix}_age.txt")
                print(age_out,locus_info['age_flank'])

                re_run_dist = 150
                Hallmarks.process_age(age_out,locus_info)
                print(locus_info)
                
                if not locus_info['Refined_Filled_start'] == "N/A":#this will allow the loop to iterate even if no alignment is found at the current alignment length
                #check if either locus is far from the expected locus. 
                    filled_start_comp = abs(locus_info['Filled_start']+1-locus_info["Refined_Filled_left_bound"])
                    filled_end_comp = abs(locus_info['Filled_end']+1-locus_info["Refined_Filled_right_bound"])
                    if locus_info['Empty_start'] == locus_info['Empty_end']:
                        print()
                        empty_start_comp = locus_info['Empty_start']-locus_info["Refined_Empty_left_bound"]
                        empty_end_comp = locus_info['Empty_end']-locus_info["Refined_Empty_right_bound"]
                    else:
                        sys.exit("need to get this working")#This is not relevant and does not arise in my canine dataset.
                    print(filled_start_comp,filled_end_comp,empty_start_comp,empty_end_comp)
                if locus_info['Refined_Filled_start'] == "N/A" or max(filled_start_comp, filled_end_comp, empty_start_comp, empty_end_comp) >= re_run_dist:
                    if locus_info['age_flank'] >= 250:
                        locus_info['age_flank'] = int(round(locus_info['age_flank'] / 2,0))
                    else:
                        #locus_info['No_predicted_locus'] = True
                        loop_break = False

                else:
                    loop_break = False          
            
            if locus_info['Refined_Filled_start'] == "N/A":
                seperate_output = [locus_info['chr'],locus_info['Filled_start'],locus_info['Filled_end'],locus_info['Empty_start'],locus_info['Empty_end'],ref,query,TE,locus_info["Age_Insertion_Found"]]
                seperate_output = [str(x) for x in seperate_output]
                seperate_output = "\t".join(seperate_output) + "\n"
                seperate_loci.write(seperate_output)
                continue
                
            #Determine if there is less than 30bp between the TSDs. If there is, I will skip the locus.
            if locus_info['right_TSD_filled'][0] != "N/A":
                left = int(locus_info['left_TSD_filled'][1])+1 
                right = int(locus_info['right_TSD_filled'][0])-1
                print(f"Right and left are {right} and {left}")

            else:
                left = int(locus_info['Refined_Filled_start'])+1
                right = int(locus_info['Refined_Filled_end'])-1
                print(f"Right and left are {right} and {left}")            
                
            if locus_info['Refined_Filled_right_bound']-locus_info['Refined_Filled_left_bound']+1 < 50 or right-left+1 < 30:
                if locus_info['Refined_Filled_right_bound']-locus_info['Refined_Filled_left_bound']+1 < 50:
                    fail_reason = "Less_than_50_bp_in_sequence"
                elif right-left+1 < 30:
                    fail_reason = "Less_than_30_bp_between_TSDs"
                else:
                    sys.exit("Invalid stop reason")
                seperate_output = [locus_info['chr'],locus_info['Filled_start'],locus_info['Filled_end'],locus_info['Empty_start'],locus_info['Empty_end'],ref,query,TE,fail_reason]
                seperate_output = [str(x) for x in seperate_output]
                seperate_output = "\t".join(seperate_output) + "\n"
                seperate_loci.write(seperate_output) 
                continue
            
                
            #mask the locus with RepeatMasker
            if TE == "LINE":
                extract_seq = f"samtools faidx {locus_info['Ref']} {locus_info['chr']}:{left}-{right} > {extract_file_prefix}_internal_extract.fa"
                subprocess.run(extract_seq,shell=True)
                print(extract_seq)
                
                #remove existing RM if exists"
                if os.path.exists(f"{extract_file_prefix}_internal_extract.fa.out"):
                    rm_cmd = f"rm {extract_file_prefix}_internal_extract.fa.*"
                    subprocess.run(rm_cmd,shell=True)
                mask_cmd = f"RepeatMasker --species dog {extract_file_prefix}_internal_extract.fa"
                subprocess.run(mask_cmd,shell=True)
                print(mask_cmd)
                
                with open(f"{extract_file_prefix}_internal_extract.fa.out",'rt') as RM_in:
                    RM_coords = []
                    for line in RM_in.readlines()[3::]:
                        line = line.split()
                        print(line)
                        if line[10] != TE_class:
                            #print(line[10])
                            #sys.exit()
                            continue
                        RM_coords.append([int(line[5]),int(line[6])])
                if len(RM_coords) == 0:
                    total_length_of_TE_sequence = 0
                    print("No L1 sequence Identified at this locus")
                    #sys.exit()
                    Multiple = False
                elif len(RM_coords) == 1:
                    total_length_of_TE_sequence = RM_coords[0][1] - RM_coords[0][0] + 1
                    print(f"One LINE-1 at this locus", RM_coords,total_length_of_TE_sequence)
                    Multiple = False
                    #sys.exit()
                else:#this implicitly assumes that a single nucleotide cannot be part of 3 repeats. To my knowledge this isn't possible.
                    total_overlap = 0
                    total_length_of_TE_sequence = 0
                    for i in range(len(RM_coords)-1):
                        if RM_coords[i+1][0] <= RM_coords[i][1]:
                            total_overlap = total_overlap+ (RM_coords[i+1][0]-RM_coords[i][1])+1
                        total_length_of_TE_sequence += RM_coords[i][1]- RM_coords[i][0] +1
                    total_length_of_TE_sequence += RM_coords[-1][1]- RM_coords[-1][0] +1

                    #print(total_overlap,RM_coords,total_length_of_TE_sequence)
                    total_length_of_TE_sequence -= total_overlap
                    print(total_overlap,RM_coords,total_length_of_TE_sequence)
                    print("Multiple LINE-1 sequences present at this locus")
                    Multiple = True
                    #sys.exit()
                if not total_length_of_TE_sequence > (right - left +1)*.7:
                   #remove because the remaining locus does not resolve as 70% LINE-1 as called by RepeatMasker
                    seperate_output = [locus_info['chr'],locus_info['Filled_start'],locus_info['Filled_end'],locus_info['Empty_start'],locus_info['Empty_end'],ref,query,TE,"Less_than_70%_TE_sequence_by_RM"]
                    seperate_output = [str(x) for x in seperate_output]
                    seperate_output = "\t".join(seperate_output) + "\n"
                    seperate_loci.write(seperate_output)
                    print("Too short",total_length_of_TE_sequence,(right-left+1)*.7)
                    #sys.exit()
                    continue
                print("Sufficient LINE-1 content",total_length_of_TE_sequence,(right-left+1)*.7)

                
            
                
            if locus_info['right_TSD_filled'][0] != "N/A":
                left = int(locus_info['left_TSD_filled'][1]) #set for bed
                right = int(locus_info['right_TSD_filled'][0])-1
            else:
                left = int(locus_info['Refined_Filled_start'])-1 #set for bed
                right = int(locus_info['Refined_Filled_end'])
            
            
            #Lets try extracting from the genomewide RM file (Do not ru-run)          
            with open("one_locus.bed",'wt') as one_locus:
                one_locus.write(f"{locus_info['chr']}\t{left}\t{right}\n")
                print(f"{locus_info['chr']}\t{left}\t{right}\n")
            intersect_cmd = f"bedtools intersect -wa -a {locus_info['chr']}.tmp -b one_locus.bed"
            rm_out = Hallmarks.get_info_from_subprocess(intersect_cmd)#Originally wanted to do wao to not lose information. Don't think it is necessary here.
            print(rm_out)
            rm_out = rm_out.split("\n")[0:-1]
            rm_out = [x.split() for x in rm_out]
            print(rm_out,len(rm_out))
            #handle the overlap
            TE_extends = False
            if len(rm_out) != 0:
                if int(rm_out[0][1]) + 20  < left: #TE extends beyond the left flank of the variant (before TSD)
                    TE_extends = True
                if int(rm_out[-1][2]) - 20 > right: #TE extends beyond the right flank of the variant (before TSD)
                    TE_extends = True

                #correct the RM loci as if -wa was not used
                if int(rm_out[0][1]) < left:
                    rm_out[0][1] = left
                if int(rm_out[-1][2]) > right:
                    rm_out[-1][2] = right
                print(rm_out)
            locus_info["TE_extends"] = TE_extends
            print(TE_extends)
            #subfamily analysis
            orientations = []
            subfamilies = []
            if len(rm_out) != 0:
                for rm in rm_out:
                    if rm == []:
                        continue
                    print(rm)
                    if int(rm[2])-int(rm[1]) > 20:
                        orientations.append(rm[8])
                        subfamilies.append(rm[3])
            locus_info['RM_subfamilies'] = subfamilies 
            if len(set(orientations)) > 1:
                print(orientations)
                locus_info['RM_orientation'] = "Multiple"
                locus_info["poly_a_list"] = [['N/A']]
            elif len(set(orientations)) == 0:
                locus_info['RM_orientation'] = "None"
                locus_info["poly_a_list"] = [['N/A']]
            else:
                if orientations[0] == "+":
                    locus_info['RM_orientation'] = "Forward"
                    locus_info["poly_a_list"] = [['N/A']]
                if orientations[0] == "C":
                    locus_info['RM_orientation'] = "Reverse"
                Hallmarks.homopolymer_extraction(locus_info)
                
                
            if "N/A" not in locus_info['left_TSD_filled'] and len(locus_info['Age_TSD_Seq_left']) >= 5:
                if locus_info['RM_orientation'] == "Forward":
                    reformatted_TSD = f"{locus_info['chr']}:{locus_info['left_TSD_filled'][0]}-{locus_info['left_TSD_filled'][1]}"
                    endonuclease_list = [reformatted_TSD,"","",""]
                    locus_info["age_endo"] = Hallmarks.identify_endo_site(locus_info["RM_orientation"],endonuclease_list,locus_info["Ref"])
                elif locus_info["RM_orientation"] == "Reverse":
                    reformatted_TSD = f"{locus_info['chr']}:{locus_info['right_TSD_filled'][0]}-{locus_info['right_TSD_filled'][1]}"
                    endonuclease_list = [locus_info['chr'],"",reformatted_TSD,""]
                    locus_info["age_endo"] = Hallmarks.identify_endo_site(locus_info["RM_orientation"],endonuclease_list,locus_info["Ref"])
                else:
                    locus_info["age_endo"] = "N/A"
            else:
                locus_info['age_endo'] = "N/A"   

            print(locus_info)
            for key in locus_info.keys():
                #print(key)
                if type(locus_info[key]) == str and  " " in locus_info[key]:
                    locus_info[key] = (locus_info[key].split())
                if type(locus_info[key]) == int or type(locus_info[key]) == float:
                    locus_info[key] = str(locus_info[key])
                elif type(locus_info[key]) == list:
                    locus_info[key] = [str(x) for x in locus_info[key]]
                    locus_info[key] = "\\".join(locus_info[key])
                else:
                    locus_info[key] = str(locus_info[key])
            
            print(required_keys)
            if required_keys == []:
                required_keys = list(locus_info.keys())
            else:
                for key in required_keys:
                    if key not in locus_info.keys():
                        sys.exit(f"key {key} missing from keys list")
                for key in locus_info.keys():
                    if key not in required_keys:
                        sys.exit(f"key {key} present in output but absent in required keys")
            
            output = []
            starting_keys = ['chr','Filled_start','Filled_end','Empty_start','Empty_end']
            if header == []:
                for key in starting_keys:
                    header.append(key)
            for key in starting_keys:
                output.append(locus_info[key])
                del locus_info[key]
            
            removable_keys = ['Insertion_sequence']
            for key in removable_keys:
                del locus_info[key]
                
            for key in sorted(locus_info.keys()):
                if " " in locus_info[key]:
                    locus_info[key] = locus_info[key].replace(" ","")
                if len(locus_info[key]) == 0:
                    locus_info[key] = "N/A"
                output.append(locus_info[key])
            if type(header) ==  list:
                for key in sorted(locus_info.keys()):
                    header.append(key)
                header = "#" + "\t".join(header) + "\n"
                processed_file.write(header)
            output = "\t".join(output) + "\n"
            print(output)
            
            processed_file.write(output)
            processed_file.flush()
            #z+=1
            
    #clean up tmp RM files
    subprocess.run("rm chr*.tmp",shell=True)
    seperate_loci.close()
    os.chdir("..")
