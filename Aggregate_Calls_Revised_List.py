#prep the upset plot
#needs bedtools
import os
import subprocess
import sys

processed_dir = "/home/blacksmi/links/kidd-lab-umms/matt-projects/Former_KiddLabScratch_Data/inter-genome_comparisons/Aligning_Using_2.26/Revised_All_Canines_List_2025_08_06"
canines = ["Clu","Nala","Sandy","Tasha","Zoey","Mischka"]
reference_genome = 'mCanlor'

cmd = ""
os.chdir(processed_dir)

subdir_of_interest = "Locus_Sharing"
if not os.path.exists(subdir_of_interest):
    os.mkdir(subdir_of_interest)
    
def acceptable_chroms():
    a_chroms = []
    for i in range(1,39):
        a_chroms.append(f"chr{i}")
    a_chroms.append("chrX")
    print(a_chroms)
    return a_chroms
chrom_list = acceptable_chroms()
#get mCanlor chrom lengths.
output = []
with open("/home/blacksmi/links/kidd-lab-umms/from-locker/genomes/mCanLor1.2/ref/mCanLor1.2.fa.fai",'rt') as infile:
    for line in infile.readlines():
        line = line.split()
        if line[0] in chrom_list:
            output.append(f"{line[0]}\t0\t{line[1]}\n")
#print(output)
with open(f"{subdir_of_interest}/mCanlor_chroms.bed",'wt') as outfile:
    for line in output:
        outfile.write(line)


subprocess.run("module list",shell=True) #bedtools is the important one here
list_of_genomes = ["Mischka","Nala","Tasha","Zoey","Sandy","Clu"]

def bedtools_analysis(dir_of_interest,input_list,ignore_list,subdir_of_interest=None):
    revised_list = []
    if ignore_list == ['None']:
        outfile = open(f"{subdir_of_interest}/Intersected_Filtered_Retro_Loci.bed",'wt')
    for item in input_list:
        if item not in ignore_list:
            revised_list.append(item)
    cmd = f"bedtools intersect -a {dir_of_interest}/{revised_list[0]}_filtered_retro.bed -b {dir_of_interest}/{revised_list[1]}_filtered_retro.bed"
    if len(revised_list) >2:
        for genome in revised_list[2::]:
            cmd += f" | bedtools intersect -a - -b {dir_of_interest}/{genome}_filtered_retro.bed "
    print(revised_list)
    print(cmd)
    process =subprocess.Popen(cmd,shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    stdout, stderr = process.communicate()
    stdout = stdout.decode()
    total = 0
    for line in stdout.split("\n"):
        if ignore_list == ['None'] and len(line.split()) > 0:
            outfile.write(f"{line}\n")
        line = line.split()
        if line!=[]:
            #print(line)
            total += int(line[2]) - int(line[1])

    print(total)
    if ignore_list == ['None']:
        print('outfile written')
        outfile.close()

        
    return total
dict_of_overlap = {}
dict_of_overlap['None'] = bedtools_analysis(".",list_of_genomes,['None'],subdir_of_interest)

cmd = f"bedtools subtract -a {subdir_of_interest}/mCanlor_chroms.bed -b {subdir_of_interest}/Intersected_Filtered_Retro_Loci.bed > {subdir_of_interest}/Inverse_Intersected_Filtered_Loci.bed"
subprocess.run(cmd,shell=True)


#create files with reduced columns and place all variants in mCanlor coordinates
contents = os.listdir(".")
for content in contents:
    calls = []
    if reference_genome not in content or not content.endswith(f"INEs_processed_filtered.bed"):
        continue
    split_content = content.split("_")
    if split_content[0] not in canines:
        continue
        
    def extract_coordinates(line):
        line = line.split()
        if line[19] == "N/A":
            empty_start = min(int(line[21]),int(line[20]))
            empty_end = max(int(line[21]),int(line[20]))-1
        else:
            empty_start = int(line[22])-1
            empty_end = int(line[19])
        filled_start = int(line[26])-1
        filled_end = int(line[23])
        chrom = line[0]   
        return chrom, str(empty_start), str(empty_end), str(filled_start), str(filled_end)
    with open(f"{content}",'rt') as infile:
        if "ref" in content:
            for line in infile.readlines():
                if line.startswith("#"):
                    continue
                
                #print(line.split()[19:27])
                #sys.exit()
                chrom, empty_start, empty_end, filled_start, filled_end = extract_coordinates(line)
                
                new_line = "\t".join([chrom,filled_start,filled_end,f"{chrom}:{empty_start}-{empty_end}",line[6]])
                new_line = f"{new_line}\tRef_Filled\t{content.split('_')[1]}/{content.split('_')[0]}\n"
                calls.append(new_line)
                #print(new_line)
        elif "query" in content:
            for line in infile.readlines():
                if line.startswith("#"):
                    continue
                chrom, empty_start, empty_end, filled_start, filled_end = extract_coordinates(line)
                new_line = "\t".join([chrom,empty_start,empty_end,f"{chrom}:{filled_start}-{filled_end}",line[6]])
                new_line = f"{new_line}\tQuery_Filled\t{content.split('_')[0]}\n"
                #print(new_line)
                calls.append(new_line)
    with open(f"{subdir_of_interest}/{content[0:-4]}_truncated.bed",'wt') as outfile:
        for call in calls:
            outfile.write(call)

subdir_contents = os.listdir(subdir_of_interest)
for file in subdir_contents:
    if not file.endswith("INEs_processed_filtered_truncated.bed"):
        continue
    cmd = f"bedtools intersect -v -a {subdir_of_interest}/{file} -b {subdir_of_interest}/Inverse_Intersected_Filtered_Loci.bed > {subdir_of_interest}/{file[0:-4]}_filtered_retro.bed"
    subprocess.run(cmd,shell=True)

       
#concatenate relevant files
subdir_contents = os.listdir(subdir_of_interest)
LINE_command = "cat "
SINE_command = "cat "
for file in subdir_contents:
    if not file.endswith("filtered_retro.bed"):
        continue
    elif "LINEs" in file:
        LINE_command = f"{LINE_command} {subdir_of_interest}/{file} "
    elif "SINEs" in file:
        SINE_command = f"{SINE_command} {subdir_of_interest}/{file} "
LINE_command = f"{LINE_command} > {subdir_of_interest}/aggregated_hits_LINEs.bed"
SINE_command = f"{SINE_command} > {subdir_of_interest}/aggregated_hits_SINEs.bed"
subprocess.run(LINE_command,shell=True)
subprocess.run(SINE_command,shell=True)
print(LINE_command,SINE_command)

TEs = ["LINEs","SINEs"]
for TE in TEs:
    cmd = f"sort -k 1,1 -k 2,2n -k 3,3n {subdir_of_interest}/aggregated_hits_{TE}.bed > {subdir_of_interest}/aggregated_hits_{TE}.sorted.bed && bedtools merge -i {subdir_of_interest}/aggregated_hits_{TE}.sorted.bed -d 100 -c 7,4 -o collapse,collapse > {subdir_of_interest}/aggregated_hits_{TE}_meged.bed"
    subprocess.run(cmd,shell=True)
