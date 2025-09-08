import pandas as pd
import numpy as np
import Python_Functions_Post_Processing.Blacksmith_Post_Processing_Utils  as PFPP
import sys
import pkg_resources

def process_file(input_file):
    infile = open(input_file,'rt')
    columns = []
    rows = []
    lines = []
    for line in infile.readlines():#uses all variants in the dataset. This includes "RM fail" 
        #print(line)
        if line.startswith("#"):
            #print(line[1::])
            for element in enumerate(line[1::].split()):
                rows.append(element[0])
                columns.append(element[1])
        else:
            lines.append(line.split())
        df = pd.DataFrame(lines,columns=list(columns))
        pandas_columns = len(df.columns.tolist())
    print(df.shape)
    infile.close()
    return df
    

loaded_modules = PFPP.automatically_extract_modules(globals().copy(),sys.modules.keys())
for module in loaded_modules:
    if module == "builtins":
        continue
    try:   
        print(f"Version of {module} is {pkg_resources.get_distribution(module).version}")
    except:
        print(f"No version information for {module}")

def identify_if_deletion_in_age(path):
    infile = open(path,'rt')
    file_content = infile.readlines()
    
    for i in range(len(file_content)):
        if file_content[i].startswith("EXCISED REGION(S)"):
            EXCISED_file1 = file_content[i+1].replace("]","").replace("[","").replace(","," ").split()
            EXCISED_file2 = file_content[i+2].replace("]","").replace("[","").replace(","," ").split()
            #print(EXCISED_file1,EXCISED_file2)
            #print(EXCISED_file2[3])
            if EXCISED_file2[3] == "0":
                continue
            #print(int(EXCISED_file2[3]),path)
            return int(EXCISED_file2[3])
            #print(EXCISED_file1,EXCISED_file2)
            #if EXCISED_file2[3] != "18":
            infile.close()
            break
            #for content in file_content[i+3::]:
            #    print(content.rstrip())
            #print(path)
    #for line in infile:
    #    print(line)
    return "N/A"
    infile.close()

input_dir = "Revised_All_Canines_List_2025_08_06"
query_genome = "Mischka"
ref_genome = "mCanlor"
state = "query"
TE = "LINE"
input_file = f"{input_dir}/{query_genome}_{ref_genome}_SV_{state}_filled_intersect_with_{TE}s_processed_filtered.bed"
#print(input_file)

df = process_file(input_file)
dict_of_locus_counts = {}
dict_of_deletion_lengths = {}
for i in range(0,11):
    dict_of_locus_counts[i] = 0
    dict_of_deletion_lengths[i] = []
subdir = f"{query_genome}_{state}_{TE}"
loci_with_deletions = []
for row in df.itertuples():
    deletion_len = identify_if_deletion_in_age(f"{input_dir}/{subdir}/{row.chr}-{row.Filled_start}-{row.Filled_end}_age.txt")
    #print(deletion_len)
    
    if "N/A" in row.Age_TSD_Seq_right:
        TSD = 0
        dict_of_locus_counts[TSD] +=1
    elif len(row.Age_TSD_Seq_right) < 10:
        TSD = len(row.Age_TSD_Seq_right)
        dict_of_locus_counts[TSD] +=1
    else:
        TSD = 10
        dict_of_locus_counts[TSD] +=1
    if deletion_len != "N/A":# and deletion_len != 1:
        dict_of_deletion_lengths[TSD].append(deletion_len)
        loci_with_deletions.append([row.chr,row.Filled_start])
print(dict_of_deletion_lengths)

print("TSD_length\tNumber_of_LINEs\tNumber_with_deletion\tavg_deletion_length\tminimum_deletion_length\tmaximum_deletion_length\tmedian_deletion_length")
for i in range(0,11):
    if len(dict_of_deletion_lengths[i]) != 0:
        avg = sum(dict_of_deletion_lengths[i])/len(dict_of_deletion_lengths[i])
        minimum = min(dict_of_deletion_lengths[i])
        maximum = max(dict_of_deletion_lengths[i])
        median = np.median(dict_of_deletion_lengths[i])
    else:
        avg = "N/A"
        minimum = "N/A"
        maximum = "N/A"
        median = "N/A"
    if i == 10:
        mod = str(i) + "+"
    else:
        mod = i
    print(f"{mod}\t{dict_of_locus_counts[i]}\t{len(dict_of_deletion_lengths[i])}\t{avg}\t{minimum}\t{maximum}\t{median}")
    

#create a file of loci with deletion
Age_out_path = "Age_Identify_target_site_deletion"
with open(f"{Age_out_path}/{query_genome}_{ref_genome}_{TE}.txt",'wt') as out:
    for deletion in loci_with_deletions:
        out_line = f"{deletion[0]}\t{deletion[1]}\n"
        out.write(out_line)
        
#run the command.
#cmd = f"run-job.py --cmd 'Process_Locus_for_Diagram.py --Processed_File {input_file} --Locus_File {Age_out_path}/{query_genome}\t{ref_genome}_{TE}.txt --outdir {Age_out_path}'"
