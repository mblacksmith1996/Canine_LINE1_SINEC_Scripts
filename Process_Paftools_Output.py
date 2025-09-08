import os
import subprocess
import sys
import argparse

parser = argparse.ArgumentParser(prog="Performing_Alignments.py",description="Produces a list of minimap2 alignment commands.")
parser.add_argument("--sample_file",dest="sample_file",required=True,help="Formatted input file where each line contains (in a tab delimited list) a sample name, the reference path, a RepeatMasker path, the path to reference Gaps, and a path to segmental duplications in the reference")
parser.add_argument("--reference_sample",dest = "aln_ref",required=True,help="Sample to which all others will be aligned. Must be present in the sample file")
parser.add_argument("--out_path_name",dest="out_dir",required=True,help="The directory that contents will be placed in. Currently, this is a subdirectory of the directory the script is launched from")
parser.add_argument("--RM_path",dest="RM_path",required=True,help="THe directory that contains RepeatMasker files. Currently, this is a subdirectory of the directory the script is launched from")
parser.add_argument("--species_name",dest="species",default="dog",required=False,help="The species being processed. Currently only supports 'dog' and 'human'")
args = parser.parse_args()

def extract_input_information(sample_file):
    canines = []
    files_dict = {}
    with open(sample_file,'rt') as infile:
        for line in infile:
            if line.startswith("#"):
                continue
            split_line = line.split()
            files_dict[split_line[0]] = split_line[1::]
            canines.append(split_line[0])
    return files_dict,canines
    
def get_TE_RMs(RM_path):
    TEs = ["LINE","SINE"]
    LINE_1_dict = {}
    SINE_dict = {}


    #generate RM dicts
    for file in os.listdir(RM_path):
        split_file = file.split("_")
        if split_file[0] in canines and split_file[1] == "LINE":
            LINE_1_dict[split_file[0]] = f"{RM_path}/{file}"
        elif split_file[0] in canines and "SINE" in split_file[1]:
            SINE_dict[split_file[0]] = f"{RM_path}/{file}"
    return SINE_dict,LINE_1_dict, TEs
    
def chrom_and_mut_rate(species):#no longer utilizes mutation rate. Not needed in this script.
    #add filter for chromosome
    chrs = []
    if args.species.lower() == "dog":
        chr_count = 38
    else:
        sys.exit("Invalid species name. Currently supported species are dog") #if this gets longer we can replace this step with a dictionary 
    for i in range(1,chr_count+1):
        chrs.append(f"chr{i}")
    chrs.append("chrX")
    #print(chrs)
    return chrs
    
    
def create_gap_and_segdup_filter_file(sample_file):
    combined_dict = {}
    with open(sample_file,'rt') as infile:
        for line in infile:
            if line.startswith("#"):
                continue
    #for key in gaps_dict.keys():
            total = 0
            split_line = line.split()
            
            #new method for TRF, gaps, segdups
            present_files = []
            for i in [3,4,6]:
                if split_line[i] != "N/A":
                    present_files.append(split_line[i])
            if present_files != []:
                with open(f"{split_line[0]}_combined.bed",'wt') as combined_out:
                    for fl in present_files:
                        print(fl)
                        with open(fl,'rt') as filtering_file:
                            for locus in filtering_file:
                                locus = "\t".join(locus.split("\t")[0:3]) + "\n"
                                combined_out.write(locus)
                    #sys.exit()
                if len(present_files) > 1:
                    cmd = f"bedtools sort -i {split_line[0]}_combined.bed | bedtools merge -i - > {split_line[0]}_combined.merge.bed"
                    subprocess.run(cmd,shell=True)
                    Retro_Filter = f"{split_line[0]}_combined.merge.bed"
                else:
                    Retro_Filter = f"{split_line[0]}_combined.bed"
                    #sys.exit()
            else:   
                Retro_Filter = "N/A"

            if split_line[5] != "N/A" and Retro_Filter != "N/A":
                cmd = f"cp {Retro_Filter} {split_line[0]}_SNP_combined.bed"
                subprocess.run(cmd,shell=True)
                with open(f"{split_line[0]}_SNP_combined.bed",'at') as output_file:
                    with open(split_line[5],"rt") as trf_input:
                        for line in trf_input:
                            line = "\t".join(line.split()[0:3]) + "\n"
                            output_file.write(line)
                        
                cmd = f"bedtools sort -i {split_line[0]}_SNP_combined.bed | bedtools merge -i - > {split_line[0]}_SNP_combined.merge.bed"
                print(cmd)
                subprocess.run(cmd,shell=True)
                SNP_Filter = f"{split_line[0]}_SNP_combined.merge.bed"
            elif split_line[5] != "N/A":
                SNP_Filter = split_line[5]
            elif Retro_Filter != "N/A":
                SNP_Filter = Retro_Filter
            else:
                SNP_Filter = "N/A"
            combined_dict[split_line[0]] = [Retro_Filter, SNP_Filter]
    return(combined_dict)
    
    
#START
original_dir = os.getcwd()
#get into working dir
os.chdir(args.out_dir)

files_dict, canines = extract_input_information(args.sample_file)
SINE_dict,LINE_1_dict,TEs = get_TE_RMs(args.RM_path)
chrs = chrom_and_mut_rate(args.species)
combined_dict = create_gap_and_segdup_filter_file(args.sample_file)


contents = os.listdir(".")

for content in contents:
    if content.endswith("paf-vars"):
        split_content = content.split(".")[0].split("_")
        query_sample = split_content[0]
        reference_sample = split_content[2]
        print(query_sample)
        
        #initialize counters
        total_autosome_length = 0
        refined_autosome_length = 0
        unfiltered_SNP_count = 0
        
        infile = open(content,'rt')
        outfile = open(f"{content}.regions.bed",'wt')
        snp_file_name = f"{split_content[0]}_{split_content[2]}_SNPs.unfiltered.txt"
        SV_file_name = f"{split_content[0]}_{split_content[2]}_SVs.unfiltered.txt"
        SNP_out = open(snp_file_name,'wt')
        SV_out = open(SV_file_name,'wt') 

        #extract regions, SNPs, and indels
        for line in infile:
            line = line.split()
            if line[1] not in chrs:
                continue
            if line[0] == "R":
                total_autosome_length += (int(line[3])-int(line[2]))
                refined_line = "\t".join(line[1::])+"\n"
                outfile.write(refined_line)
            else:    
                if line[4] != "1": #suggested by minimap2 manual
                    continue
                if line[1] != line[8]:
                    continue  
                    
                #isolate SNPS
                if len(line[6]) == 1 and len(line[7]) == 1 and line[6] != "-" and line[7] != "-": 
                    line_out = "\t".join(line[1::]) + "\n"
                    SNP_out.write(line_out)
                    unfiltered_SNP_count+=1
                    
                #isolate structural variants 50bp or larger 
                elif len(line[6]) < 50 and len(line[7]) < 50:
                        continue
                elif len(line[6]) >= 50 or len(line[7]) >= 50:
                    line_out = [line[1],line[2],line[3],line[8],line[9],line[10],line[6],line[7],line[11]]
                    line_out = "\t".join(line_out) + "\n"
                    SV_out.write(line_out)
                else:
                    sys.exit("Invalid Locus")
                    
        infile.close()
        outfile.close()
        SNP_out.close()
        SV_out.close()
        
        
        ##filter the aligned regions.
        
        #calculate refined autosome length
        cmd = f"bedtools subtract -a {content}.regions.bed -b {combined_dict[args.aln_ref][1]} > {query_sample}_filtered.bed"
        print(cmd)
        subprocess.run(cmd,shell=True)
        
        #also generate a filtered length using the SV filter. This is used for getting loci that are present in all samples. Not used here.
        cmd = f"bedtools subtract -a {content}.regions.bed -b {combined_dict[args.aln_ref][0]} > {query_sample}_filtered_retro.bed"
        print(cmd)
        subprocess.run(cmd,shell=True)
        
        with open(f"{query_sample}_filtered.bed",'rt') as infile:
            for line in infile:
                line = line.split()
                if line[0] != "chrX":
                    refined_autosome_length += (int(line[2])-int(line[1]))
                
        #Filter SVs
        cmd = f"bedtools window -w 100 -v -a {SV_file_name} -b {combined_dict[args.aln_ref][0]} > {SV_file_name}.reference_filtered"
        print(cmd)
        subprocess.run(cmd,shell=True)
        
        with open(f"{SV_file_name}.reference_filtered",'rt') as infile:
            with open(f"{SV_file_name}.reference_filtered_reordered",'wt') as outfile:
                for line in infile:
                    line = line.split()
                    new_line = [line[3],line[4],line[5]]
                    for item in line:
                        new_line.append(item)
                    out_line = "\t".join(new_line) + "\n"
                    outfile.write(out_line)
        cmd = f"bedtools window -w 100 -v -a {SV_file_name}.reference_filtered_reordered -b {combined_dict[query_sample][0]} > {SV_file_name}.reference_filtered_reordered_query_filtered"
        subprocess.run(cmd,shell=True)
        print(cmd)
        
        #split based on which sample is the filled site
        query_filled = []
        ref_filled = []
        with open(f"{SV_file_name}.reference_filtered_reordered_query_filtered",'rt') as infile:
            for line in infile:
                line = line.split()[3::]
                if line[6] == "-":
                    reformatted_line = [line[3],line[4],line[5],line[0],line[1],line[2],line[7],line[8]]
                    query_filled.append(reformatted_line)
                elif line[7] == "-":
                    reformatted_line = [line[0],line[1],line[2],line[3],line[4],line[5],line[6],line[8]]
                    ref_filled.append(reformatted_line)
        ref_file = f"{query_sample}_{reference_sample}_SV_ref_filled.txt"
        query_file = f"{query_sample}_{reference_sample}_SV_query_filled.txt"
        
        with open(query_file,'wt') as qf:
            for item in query_filled:
                item = "\t".join(item) + "\n"
                qf.write(item)
            
        with open(ref_file,'wt') as wf:
            for item in ref_filled:
                item = "\t".join(item) + "\n"
                wf.write(item)
                

        cmd = f"bedtools intersect -wa -f .7 -a {ref_file} -b {LINE_1_dict[split_content[2]]} | sort -k1,1 -k2,2n > {ref_file[0:-4]}_intersect_with_LINEs.bed"
        print(cmd)
        subprocess.run(cmd,shell=True)
        
        cmd = f"bedtools intersect -wa -f .7 -a {ref_file} -b {SINE_dict[split_content[2]]} | sort -k1,1 -k2,2n > {ref_file[0:-4]}_intersect_with_SINEs.bed"
        print(cmd)
        subprocess.run(cmd,shell=True)

       
        cmd = f"bedtools intersect -wa -f .7 -a {query_file} -b {LINE_1_dict[split_content[0]]} | sort -k1,1 -k2,2n > {query_file[0:-4]}_intersect_with_LINEs.bed"
        print(cmd)
        subprocess.run(cmd,shell=True)    
        
        cmd = f"bedtools intersect -wa -f .7 -a {query_file} -b {SINE_dict[split_content[0]]} | sort -k1,1 -k2,2n > {query_file[0:-4]}_intersect_with_SINEs.bed"
        print(cmd)
        subprocess.run(cmd,shell=True)   
    

        #remove SNPs from duplicated regions, low complexity regions, and gaps
        cmd = f"bedtools intersect -a {snp_file_name} -b {query_sample}_filtered.bed > {split_content[0]}_{split_content[2]}_SNPs.filtered.txt"
        subprocess.run(cmd,shell=True)
        print(cmd)
        final_SNP_count = 0
        with open(f"{split_content[0]}_{split_content[2]}_SNPs.filtered.txt") as filtered_SNP:
            for line in filtered_SNP:
                line = line.split()
                if not line[0] == "chrX":
                    final_SNP_count +=1


        print("Original aligned length is " + str(total_autosome_length))
        print("Refined aligned length is" + str(refined_autosome_length))
        print("Unfiltered SNP count is" + str(unfiltered_SNP_count))
        print(f"Filtered number of SNPs for {query_sample} vs {reference_sample} is {final_SNP_count}")
        
        rate_information_dict = {}
        rate_information_dict["SNP_total"] = str(final_SNP_count)
        rate_information_dict["auto_length"] = str(refined_autosome_length)

        print(rate_information_dict)   
        
        #create output for the rate step
        with open(f"{content[0:-9]}_Rate_Information.txt",'wt') as outfile:
            header_info = []
            content_info = []
            for key in rate_information_dict.keys():
                header_info.append(key)
                content_info.append(rate_information_dict[key])
            header_info = "\t".join(header_info)+"\n"
            content_info = "\t".join(content_info)+"\n"
            outfile.write(header_info)
            outfile.write(content_info)





sys.exit()