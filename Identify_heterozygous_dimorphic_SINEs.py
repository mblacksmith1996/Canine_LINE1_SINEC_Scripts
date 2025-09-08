import os
import sys
import subprocess

canine_list = ["Mischka"]
sv_file = "/home/blacksmi/links/kidd-lab-umms/from-locker/jmkidd/people-projects/peter-projects/dog-sv-longreads/browser-tracks/2025-02-24/SV_12_sample.all.LISTTECH.bed"
out_directory = "Revised_All_Canines_List_2025_08_06"
out_subdir = "Heterozygous_Analysis"

os.chdir(out_directory)
if not os.path.exists(out_subdir):
    os.mkdir(out_subdir)
os.chdir(out_subdir)

#canine_dict = {}
#for canine in canine_list:
#    canine_dict[canine] = []

mischka_out = open("Mischka_SV_del.bed",'wt')

with open(sv_file,'rt') as infile:
    for line in infile:
        #print(line)
        for canine in canine_list:
            if "Mischka" in line and "DEL" in line:
                mischka_out.write(line)

                
mischka_out.close()


#use bedtools to intersect with the Mischka variants.
cmd = "bedtools intersect -wa -a Mischka_SV_del.bed -b ../Mischka_mCanlor_SV_query_filled_intersect_with_SINEs.bed -f .9 -r > Mischka_SV_del_intersected.bed"
subprocess.run(cmd,shell=True)

SR_hets = []
LR_hets = []
PA_hets = []
with open("Mischka_SV_del_intersected.bed",'rt') as infile:
    for line in infile.readlines():
        line = line.split()
        print(line)
        zygosity = line[-1].split(";")
        if line[0] == "chrX":
            continue
        zygosity_out = []
        for zygo in zygosity:
            if zygo.startswith("Mischka"):
                zygosity_out.append(zygo)
        
        LR=0
        SR=0
        PA=0
        for zygo in zygosity_out:
            if "LR" in zygo and "0/1" in zygo:
                LR=1
            elif "SR" in zygo and "0/1" in zygo:
                SR=1
            elif "PA" in zygo and "0/1" in zygo:
                PA=1
                
        LR_hets.append(LR)
        SR_hets.append(SR)
        PA_hets.append(PA)
        #print(line)
        
print(sum(SR_hets),sum(LR_hets),sum(PA_hets))
all_agree = 0
LR_SR_agree=0
LR_PA_agree=0
SR_PA_agree=0
for i in range(len(LR_hets)):
    if SR_hets[i] == LR_hets[i] and SR_hets[i] == PA_hets[i]:
        all_agree+=1
        LR_SR_agree+=1
        LR_PA_agree+=1
        SR_PA_agree+=1
    elif SR_hets[i] == LR_hets[i]:
        LR_SR_agree+=1
    elif SR_hets[i] == PA_hets[i]:
        SR_PA_agree+=1
    elif LR_hets[i] == PA_hets[i]:
        LR_PA_agree+=1
print(all_agree,LR_SR_agree,LR_PA_agree,SR_PA_agree)
#print(SR_hets,LR_hets,PA_hets)

