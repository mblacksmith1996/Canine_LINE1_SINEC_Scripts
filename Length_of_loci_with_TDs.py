import pandas as pd
import matplotlib.pyplot as plt
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
    
    
directory = "/nfs/turbo/umms-jmkidd/matt-projects/Former_KiddLabScratch_Data/inter-genome_comparisons/Aligning_Using_2.26/Revised_All_Canines_List_2025_08_06"
with open(f"{directory}/Transduction_Detection/Transduction_list.txt",'rt') as infile:
    lengths = []
    list_of_variants = []
    for line in infile:
        line = line.split()
        if line[0] != "Mischka":
            continue
        if line[1] == "ref":
            continue
            
        variant = [line[2],line[3],line[4]]#NOTE. Only works for Mischka for now
        if variant not in list_of_variants:
            list_of_variants.append(variant)
        #print(line)
    print(list_of_variants,len(list_of_variants))
    mischka_query = f"{directory}/Mischka_mCanlor_SV_query_filled_intersect_with_LINEs_processed_filtered.bed"
    df = process_file(mischka_query)
        
    for row in df.itertuples():
        for i in range(len(list_of_variants)):
            left_bound_empty = int(row.Refined_Empty_left_bound)
            right_bound_empty = int(row.Refined_Empty_right_bound)
            
            left = min(left_bound_empty,right_bound_empty)
            right = max(left_bound_empty,right_bound_empty)-1
            left,right = str(left),str(right)
            if list_of_variants[i] == [row.chr,left,right]:
                lengths.append(int(row.Refined_Filled_end)-int(row.Refined_Filled_start)+1)
    print(list_of_variants)
    print(lengths)
    
    
    #print(x)

    plt.hist(lengths,bins=30)
    plt.xlabel("Length of Variant")
    plt.ylabel("Number of Variants")
    plt.savefig(f"{directory}/Fig_Save/Length_of_Variants_with_TDs_Mischka_filled.png",dpi=300,bbox_inches="tight")

    filtered_lengths = [50,100,500,1000,2000,3000,4000,5000,6000]
    for filtered_length in filtered_lengths:
        current = 0
        for length in lengths:
            if length >= filtered_length:
                current+=1
        print(f"{current} ({current/len(lengths)*100}% of) of variants are at least {filtered_length} bp long.")