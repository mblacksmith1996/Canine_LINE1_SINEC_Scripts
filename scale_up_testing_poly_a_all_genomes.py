import pandas as pd
import subprocess
import argparse
import os

#parse input args
parser = argparse.ArgumentParser(prog="scale_up_testing_poly_a_all_genomes.py",description="Impliment the relaex filtering on all genomes")
parser.add_argument("--output_dir",dest="out_directory",required=True,help="Path to directory where output will be placed")
parser.add_argument("--parse_file",dest="input_file",required=True,help="Path to input file")
args = parser.parse_args()


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
i = 0

def samtools_extract(path,chrom, left_bound, right_bound,mod):
    cmd = f"samtools faidx {path} {mod} {chrom}:{left_bound}-{right_bound}"
    #print(cmd)
    #subprocess.run(cmd,shell=True)
    proc = subprocess.run(cmd,shell=True, encoding='utf-8', stdout=subprocess.PIPE)
    seq = ""
    for line in proc.stdout.split('\n'):
        if line.startswith(">"):
            header = line
        else:
            seq = seq + line
    seq =  seq.lower()
    return header, seq

def scanning(sequence):
    scan_length=15
    
    max_length = len(sequence)+1
    max_a = 0
    for i in range(0,max_length-scan_length):
        #print(sequence[i:i+scan_length])
        current_a_percent = sequence[i:i+scan_length].count("a")
        #print(current_a_percent)
        if max_a < current_a_percent:
            max_a = current_a_percent
    return max_a
def process_sequence(sequence):
    max_a = scanning(sequence)
    #print(max_a)
    
    
    return max_a
def process_df(df):
    dict_of_loci = {}
    for row in df.itertuples():    
        #get genome to extract from
        if state == "query":
            samtools_genome = row.Ref
        else:
            samtools_genome = row.Query
        dist = 30
        if len(row.Age_TSD_Seq_right) >= 5:
            TSD_found = True
            TSD_len = len(row.Age_TSD_Seq_right)
            if row.RM_orientation == "Forward":
                left_bound= int(row.left_TSD_filled.split("\\")[1])+1
                #print(f"left bound of variant is {left_bound}")
                TSD_coords = row.right_TSD_filled.split("\\")
                poly_coords_left = int(TSD_coords[0])-dist
                poly_coords_right = int(TSD_coords[0])-1
                if poly_coords_left < left_bound:
                    poly_coords_left = left_bound
                head,sequence = samtools_extract(samtools_genome,row.chr,poly_coords_left,poly_coords_right, "")

                #print(head,sequence)
                #max_a = process_sequence(sequence)

            elif row.RM_orientation == "Reverse":
                right_bound = int(row.right_TSD_filled.split("\\")[0])-1
                TSD_coords = row.left_TSD_filled.split("\\")
                poly_coords_left = int(TSD_coords[1])+1
                poly_coords_right = int(TSD_coords[1])+dist
                if poly_coords_right > right_bound:
                    poly_coords_right = right_bound
                head,sequence = samtools_extract(samtools_genome,row.chr,poly_coords_left,poly_coords_right, "-i")
                #print(head,sequence)
            else:
                continue
        else:
            TSD_found = False
            TSD_len = len(row.Age_TSD_Seq_right)
            if row.RM_orientation == "Forward":
                left_bound= int(row.Refined_Filled_start)
                #TSD_coords = row.right_TSD_filled.split("\\")
                poly_coords_left = int(row.Refined_Filled_end)-(dist-1)
                poly_coords_right = int(row.Refined_Filled_end)
                if poly_coords_left < left_bound:
                    poly_coords_left = left_bound
                head,sequence = samtools_extract(samtools_genome,row.chr,poly_coords_left,poly_coords_right, "")

                #print(head,sequence)
                #max_a = process_sequence(sequence)

            elif row.RM_orientation == "Reverse":
                right_bound = int(row.Refined_Filled_end)
                TSD_coords = row.left_TSD_filled.split("\\")
                poly_coords_left = int(row.Refined_Filled_start)
                poly_coords_right = int(row.Refined_Filled_start)+(dist-1)
                if poly_coords_right > right_bound:
                    poly_coords_right = right_bound
                head,sequence = samtools_extract(samtools_genome,row.chr,poly_coords_left,poly_coords_right, "-i")
                #print(head,sequence)
            else:
                continue
        max_a = process_sequence(sequence)
        if "N" in row.Adjacent_Poly:
            poly_len = 0
        else:
            split_poly = row.Adjacent_Poly.split("\\")
            poly_len = int(split_poly[1])-int(split_poly[0])+1
        dict_of_locus = {'seq':sequence, 'max':max_a, 'TSD_found':TSD_found,'TSD_len':TSD_len, "Adjacent_Poly_length":poly_len,"RM_orientation":row.RM_orientation}
        dict_of_loci[f"{row.chr}_{row.Filled_start}_{row.Filled_end}"] = dict_of_locus
    #print(dict_of_loci)
    multiple_RM_orientations = len(df)-len(dict_of_loci.keys())
    return dict_of_loci,multiple_RM_orientations

total = 0
long_poly = 0
def restrictive_processing(dictionary):
    Both = 0
    TSD = 0
    Poly = 0
    Neither = 0
    Neither_dict = {}
    for key in dictionary.keys():
        subdict = dictionary[key]
        if subdict['TSD_len'] >= 10 and subdict["Adjacent_Poly_length"] >= 10:
            Both +=1
        elif subdict['TSD_len'] >= 10:# and subdict['Adjacent_Poly_length'] != 0:
            TSD +=1
        elif subdict["Adjacent_Poly_length"] >= 10:
            Poly +=1
        else:
            Neither +=1
            Neither_dict[key] = subdict
    return_dict = {"Both":Both,'TSD':TSD,'Poly':Poly,"Neither":Neither}#,"Neither_Dict":Neither_dict}
    return return_dict
def permissive_processing(dictionary):
    Both = 0
    TSD = 0
    Poly = 0
    Neither = 0
    Neither_dict = {}
    for key in dictionary.keys():
        subdict = dictionary[key]
        if subdict['TSD_len'] >= 7 and (subdict["Adjacent_Poly_length"] >= 10 or subdict['max'] >= 10):
            Both +=1
        elif subdict['TSD_len'] >= 7:# and subdict['Adjacent_Poly_length'] != 0:
            TSD +=1
            #print(key,subdict)
        elif subdict["Adjacent_Poly_length"] >= 10 or subdict['max'] >= 10:
            Poly +=1
        else:
            Neither +=1
            Neither_dict[key] = subdict
    return_dict = {"Both":Both,'TSD':TSD,'Poly':Poly,"Neither":Neither}#,"Neither_Dict":Neither_dict}
    return return_dict
#directory = "/nfs/turbo/umms-jmkidd/matt-projects/Former_KiddLabScratch_Data/inter-genome_comparisons/Aligning_Using_2.26/Revised_All_Canines_List_2025_08_06"

if not os.path.exists(args.out_directory):
    os.mkdir(args.out_directory)

filled_canine = "Mischka"
empty_sample = "mCanlor"
state = "query"

df = process_file(args.input_file)
print(df.columns)
dict_of_loci,multi_orientation = process_df(df)


#filter dictionary 
restrictive_dict = restrictive_processing(dict_of_loci)
restrictive_dict['RM_fail'] = multi_orientation
print(restrictive_dict,sum(restrictive_dict.values()))
permissive_dict = permissive_processing(dict_of_loci)
permissive_dict['RM_fail'] = multi_orientation
print(permissive_dict,sum(permissive_dict.values()))


#output counts
header = ["#Sample_name"]
for key in permissive_dict.keys():
    header.append(key)
header = "\t".join(header)+"\n"

file_name = args.input_file.split("/")[-1].split("_") 
sample_name = file_name[0] + "_" + file_name[3] + "_" + file_name[7]
permissive_list = [sample_name+"_permissive"]
restrictive_list = [sample_name+"_restrictive"]

for key in permissive_dict:
    permissive_list.append(str(permissive_dict[key]))
permissive_list="\t".join(permissive_list)+"\n"

for key in restrictive_dict:
    restrictive_list.append(str(restrictive_dict[key]))
restrictive_list="\t".join(restrictive_list)+"\n"

print(header,restrictive_list,permissive_list)

with open(f"{args.out_directory}/{sample_name}_counts.txt",'wt') as outfile:
    outfile.write(header)
    outfile.write(restrictive_list)
    outfile.write(permissive_list)

#also output the info held in the dictionary
dict_header = list(dict_of_loci[list(dict_of_loci.keys())[0]].keys())#get the list of keys from the first value of the dict of loci dictionary
print(dict_header)
dict_header.insert(0,"locus")
print(dict_header)
dict_header = "\t".join(dict_header)+"\n"
print(dict_header)
with open(f"{args.out_directory}/{sample_name}_dict_of_loci.txt",'wt') as out:
    out.write(dict_header)
    for key in dict_of_loci.keys():
        vals = [str(x) for x in dict_of_loci[key].values()]
        print(vals)
        combined_val = key+"\t"+"\t".join(vals)+"\n"
        out.write(combined_val)

