import argparse
import os 
import pandas as pd
import subprocess
import sys
sys.path.append("/home/blacksmi/links/kidd-lab-umms/matt-projects/Hallmark_Detection_python_Scripts_Dev/Hallmark_Detection_Python_Scripts")
import Detecting_Hallmarks_Functions as Hallmarks


#parse input args
parser = argparse.ArgumentParser(prog="Transduction_Detection.py",description="Identifies and extracts transductions from LINE-1 data")
parser.add_argument("--sample_file",dest="sample_file",required=True,help="Formatted input file where each line contains (in a tab delimited list) a sample name, the reference path, a RepeatMasker path, the path to reference Gaps, and a path to segmental duplications in the reference")
parser.add_argument("--Data_Dir",dest="aln_dir",required=True,help="The directory that contents will be placed in. Currently, this is a subdirectory of the directory the script is launched from")
parser.add_argument("--RM_path",dest="RM_dir",required=True,help="The directory that contains RepeatMasker files. Currently, this is a subdirectory of the directory the script is launched from")
args = parser.parse_args()


full_path_RM = os.path.realpath(args.RM_dir)
sdust_path = "/home/blacksmi/links/kidd-lab-umms/jmkidd-projects/ballard-assemblies/run-minigraphpan/analysis/progs/sdust/sdust"

Reference_Path_Dict = {} #extract the paths to the full references
with open(args.sample_file,'rt') as infile:
    for line in infile:
        if line.startswith("#"):
            continue
        line = line.split()
        Reference_Path_Dict[line[0]] = line[1]
print(Reference_Path_Dict)

os.chdir(args.aln_dir)
dir_contents = os.listdir(".")
if not os.path.exists("./Transduction_Detection"):
    os.mkdir("./Transduction_Detection")
blat_out = open("Transduction_Detection/Blat_CMDs.txt",'wt')
for file in dir_contents:
    if "LINEs_processed_filtered" not in file:
        continue
    #if "Mischka" not in file or "query" not in file:#Todo remove
    #    continue
    print(file)
    split_file = file.split("_")
    query = split_file[0]
    ref = split_file[1]
    state = split_file[3]

    coords_list = []
    with open(file,'rt') as infile: #create dataframe containing all information about the locus
        columns = []
        rows = []
        lines = []
        for line in infile.readlines():
            if line.startswith("#"):
                #print(line[1::])
                for element in enumerate(line[1::].split()):
                    rows.append(element[0])
                    columns.append(element[1])
            else:
                lines.append(line.split())
    df = pd.DataFrame(lines,columns=list(columns))
    pandas_columns = len(df.columns.tolist())
    print(list(df))
    print(df.shape)
    
    for row in df.itertuples(): #iterate through dataframe. If the TSD is 10bp or longer do not include the TSD in the insertion contents. If shorter, include left TSD in the contents.
        if len(row.Age_TSD_Seq_left) >=10:
            left_bound = int(row.left_TSD_filled.split("\\")[1])+1
            right_bound = int(row.right_TSD_filled.split("\\")[0])-1
            
        else:
            left_bound = int(row.Refined_Filled_start)
            right_bound = int(row.Refined_Filled_end)
            
        #extract the mCanlor coordinates for later
        ref_left = int(row.Refined_Empty_left_bound) #this is affected by orientation of each comparison relative to mCanlor. So ref_left can be larger or smaller than ref right. Handled below by subtracting 1 from the larger value. 1 is subtracted because it is the first base after the insertion and bed is half closed. 1 is not subtracted from left bound because it is the base before the insertion starts.
        ref_right = int(row.Refined_Empty_right_bound)
            
        coords = [row.chr,left_bound-1,right_bound,row.chr,min(ref_left,ref_right),max(ref_left,ref_right)-1] #-1 to put it in bed coordinates
    
        coords_list.append(coords)
    
    with open(f"Transduction_Detection/{query}_{ref}_{state}_coord_extract.bed",'wt') as out: #write to output in bedformat
        for coord in coords_list:
            coord = [str(x) for x in coord]
            coord = "\t".join(coord) + "\n"
            out.write(coord)
    if state == "query":
        sample_name = query
    else:
        sample_name = ref
    
    cmd = f"bedtools intersect -a Transduction_Detection/{query}_{ref}_{state}_coord_extract.bed -b {full_path_RM}/{sample_name}_LINE_RM.out.bed -wao > Transduction_Detection/{query}_{ref}_{state}_coord_intersect.bed"
    subprocess.run(cmd,shell=True)
    
    
    rm_start_index = 7
    rm_end_index = 8
    orientation_index = -3

    #parse the bedtools intersect file.
    sam_name = f"Transduction_Detection/{query}_{ref}_{state}_fasta_extract.fa"
    sam_out = open(sam_name,'wt')
    with open(f"Transduction_Detection/{query}_{ref}_{state}_coord_intersect.bed",'rt') as infile:
        transduction_coordinates_list = []
        intersect_contents = infile.readlines()
        #print(len(intersect_contents))
        current_content = []
        for i in range(len(intersect_contents)):
            if "Null" in intersect_contents[i]:
                print(intersect_contents[i])
                sys.exit()
            if current_content == []:
                current_content.append(intersect_contents[i].split())
            if i != len(intersect_contents)-1 and intersect_contents[i+1].split()[1] == current_content[0][1]:
                current_content.append(intersect_contents[i+1].split())
                continue
            else:
                #print(current_content)
                
                #identify the longest segment. This will be where the transductions are evaluated from.
                longest_segment_length = -1
                longest_segment_index = ""
                for i in range(len(current_content)):
                    if int(current_content[i][-1]) > longest_segment_length:
                        longest_segment_length = int(current_content[i][-1])
                        longest_segment_index = i
                #print(longest_segment_length,longest_segment_index)
                
                tmp_orientation_index = orientation_index
                while True:
                    #print(tmp_orientation_index,type(tmp_orientation_index),longest_segment_index,current_content)
                    if current_content[longest_segment_index][tmp_orientation_index] == "+":
                        #print('Forward')
                        orientation_mod = ""
                        #print(current_content[longest_segment_index])
                        transduction_coordinates = [current_content[longest_segment_index][0],int(current_content[longest_segment_index][rm_end_index]),int(current_content[longest_segment_index][2])]
                        #sys.exit()
                        break
                    elif current_content[longest_segment_index][tmp_orientation_index] == "C":
                        #print("Reverse")
                        orientation_mod = "-i"
                        #print(current_content[longest_segment_index])
                        transduction_coordinates = [current_content[longest_segment_index][0],int(current_content[longest_segment_index][1])+1,int(current_content[longest_segment_index][rm_start_index])] #is this and the above correct
                        break
                    else:
                        tmp_orientation_index -=1
                        #print(current_content[longest_segment_index])
                        #sys.exit("Invalid Orientation")
                        if tmp_orientation_index == -5:#sometimes there can be more columns. But there should be no other column that is + or C
                            print(current_content[longest_segment_index])
                            sys.exit("Invalid Orientation")
                    
                if transduction_coordinates[2]-transduction_coordinates[1] >= 24:    
                    #print(transduction_coordinates)
                    transduction_coordinates = [str(x) for x in transduction_coordinates]
                    transduction_coordinates = (f"{transduction_coordinates[0]}:{transduction_coordinates[1]}-{transduction_coordinates[2]}")
                    
                    
                    samtools_cmd = f"samtools faidx {orientation_mod} {Reference_Path_Dict[sample_name]} {transduction_coordinates}"
                    #print(samtools_cmd)
                    samtools_out  = Hallmarks.get_info_from_subprocess(samtools_cmd)
                    
                    #check if the segment is primarily degenrate A sequence. 
                    seq = ""
                    for row in samtools_out.split("\n")[1::]:  
                        seq = seq + row.rstrip()
                    #print(samtools_out.split("\n")[0],seq)
                    #a_count = 0
                    #for nuc in seq:
                    #    if nuc == "a" or nuc == "A":
                    #        a_count +=1
                    #if a_count/len(seq)*100 > 70:
                    #    pass
                    #else:
                    samtools_name = samtools_out.split("\n")[0] + " " + ",".join(current_content[longest_segment_index][0:6]) +"," + query + "," + state + "\n"
                    seq = seq + "\n"
                    sam_out.write(samtools_name)
                    sam_out.write(seq)
                        
                    
                    
                    
                current_content = []#reset counter since the next line (if there is one) is at a different locus
        #print(transduction_coordinates_list,len(transduction_coordinates_list))
        
        #for
    sam_out.close()
    
    #RepeatMask the fasta file
    #if not os.path.exists(f"{sam_name}.out"):
    RM_cmd = f"RepeatMasker --species dog {sam_name}"
    subprocess.run(RM_cmd,shell=True)
    
    #Process the fasta.out file
    #if not os.path.exists(f"{sam_name}.out.bed"):
    out_fa_bed = open(f"{sam_name}.out.bed",'wt')
    with open(f"{sam_name}.out",'rt') as infile:
        for line in infile.readlines()[3::]:
            print(line)
            line = line.split()
            coords = f"{line[4]}\t{int(line[5])-1}\t{line[6]}\n"
            out_fa_bed.write(coords)
    out_fa_bed.close()
    
    
    #if not os.path.exists(f"{sam_name}.sdust.bed"):
    sdust_cmd = f"{sdust_path} {sam_name} > {sam_name}.sdust.bed"
    print(sdust_cmd)
    subprocess.run(sdust_cmd,shell=True)

    
    #create merged file of the sdust and RM loci
    cmd = f"cat {sam_name}.out.bed {sam_name}.sdust.bed | bedtools sort -i - | bedtools merge -i - > {sam_name}.RM.sdust.bed.merge"
    subprocess.run(cmd, shell=True)
    
    #mask the fasta file
    cmd = f"bedtools maskfasta -fi {sam_name} -bed {sam_name}.RM.sdust.bed.merge -fo {sam_name}.masked"
    subprocess.run(cmd,shell=True)
    
    #Remove loci which do not have at least 25bp of unmasked sequence
    def get_next_sequence(readlines,i):
        if i==len(readlines):
            return "","", ""
        if not readlines[i].startswith(">"):
            sys.exit("sequence name does not start with '>'")
        seq_name = readlines[i].rstrip()
        i+=1
        seq = ""
        while True:
            seq = seq + readlines[i].rstrip()
            if i+1 == len(readlines):
                return seq_name, seq,i+1
            if readlines[i+1].startswith(">"):
                return seq_name, seq, i+1
            i = i+1
            
            
    loci_to_retain = []
    out_filtered_mask = open(f"{sam_name}_filtered.masked",'wt')
    with open(f"{sam_name}.masked",'rt') as infile:
        readlines = infile.readlines()
        index=0
        sequence = None
        while sequence != "":
            name,sequence,index = get_next_sequence(readlines,index)
            unmasked = 0
            for nuc in sequence:
                if nuc != "N" and nuc != "n":
                    unmasked+=1
            if unmasked >= 25:
                #print(sequence,unmasked)
                out_filtered_mask.write(f"{name}\n{sequence}\n")
                loci_to_retain.append(name)
    out_filtered_mask.close()
    #sys.exit()
    out_filtered_sam = open(f"{sam_name}.filtered",'wt')
    with open(sam_name,'rt') as infile:
        readlines = infile.readlines()
        index=0
        sequence = None
        while sequence != "":
            name,sequence,index = get_next_sequence(readlines,index)
            if name == "":
                break
            split_name= name.split()[0]
            if split_name in loci_to_retain:
                #print(name,sequence)
               out_filtered_sam.write(f"{name}\n{sequence}\n")
    out_filtered_sam.close()
    #sys.exit()
    
    blat_command = f"blat -out=pslx -minIdentity=95 {Reference_Path_Dict['mCanlor']} {sam_name}.filtered  {sam_name[0:-3]}.filtered.pslx\n"
    blat_out.write(blat_command)
    print(blat_command)
    #break
blat_out.close()

#steps.
#extract variants and get their transductions.
#repeatmask and sdusst mask the sequences.
#Blat the unmasked
#filter

