#!/usr/bin/env python3
def automatically_extract_modules(global_vars,sys_modules):
    modules = []
    for value in global_vars.values():
        if "module" not in str(type(value)) or "class" not in str(type(value)):
            continue
        if str(value).split()[1].replace("'","") in sys_modules:
           modules.append(str(value).split()[1].replace("'",""))

    #print(modules)
    return modules
    
def logo_maker(Endo_sites,pd,logomaker,mpl):
    endo_dict = {}
    for i in range(len(Endo_sites[0])):
        endo_dict[i] = {"A":0,"T":0,"C":0,"G":0}
    #print(endo_dict)
    for cut_site in Endo_sites:
        cut_site = cut_site.upper()
        for i in range(len(cut_site)):
            endo_dict[i][cut_site[i]] +=1
            if endo_dict[i][cut_site[i]] == "N":
                sys.exit()
    for key in endo_dict.keys():
        #print(endo_dict[key])
        for key2 in endo_dict[key].keys():
            #print(key2)
            endo_dict[key][key2] = endo_dict[key][key2]/len(Endo_sites)
    endo_dict
    endo_pd = pd.DataFrame(endo_dict)
    endo_pd = endo_pd.T
    print(endo_pd)
    print(len(Endo_sites))
    # create Logo object
    ss_logo = logomaker.Logo(endo_pd,
    width=.8,
    vpad=.05,
    #fade_probabilities=True,
    stack_order='big_on_top',
    color_scheme='classic')
    # style using Logo methods
    ss_logo.style_spines(spines=['left', 'right'], visible=False)
    # style using Axes methods
    ss_logo.ax.set_xticks(range(len(endo_pd)))
    ss_logo.ax.set_xticklabels('%+d'%x for x in [-5, -4, -3, -2, -1, +1, +2])
    ss_logo.ax.set_yticks([0,.2,.4,.6,.8, 1])
    ss_logo.ax.axvline(4.5, color='k', linewidth=1, linestyle=':')
    ss_logo.ax.set_ylabel('probability',fontsize="15")
    ss_logo.ax.set_xlabel("position in motif",fontsize="15")
    #mpl.pyplot.show()
    return ss_logo
    
def filter_processed_files(df, remove_indices_TE_extends_or_no_RM_orientation,processed_file):
    outname = processed_file[0:-4] + "_filtered.bed"
    output_file = open(outname,'wt')

    header = "#" + "\t".join(df.columns.tolist()) + "\n"
    output_file.write(header)

    for i in range(len(df)):
        if i in remove_indices_TE_extends_or_no_RM_orientation:
            continue
        row_info = df.iloc[i].values
        row_info[-1] = f"{row_info[-1][0]}\{row_info[-1][1]}" #hardcoded to handle the Adjacent_Poly column
        row_info = [str(x) for x in row_info]
        row_info = "\t".join(row_info)
        row_info += "\n"
        output_file.write(row_info)

        
    #input_file.close()
    output_file.close()
def create_dictionary_of_LINE1_and_SINEC_dataframes(misc_info_dict,pd,sys):
    """
    Uses pandas to generate dataframes for the hallmark containing output from hallmark testing.
    Removes variants which have a repeatmasker tract which extends 20bp or more beyond the end of the inserted sequence as called by Age.
    
    Parameters
    misc_info_dict['canines']: list 
        names of canines in study. names must correlate with file names
    misc_info_dict['Directory_of_Data']: str
        path to directory containing processed bedfiles
    pd: module
        pandas
        
    Returns
    canine_dict_LINE: dictionary where the first element of the dictionary is the pandas dataframe associated with a file of processed LINE-1 containing variants.
        The second element of the list is the numbe of variants which have multiple or no RM variants afterprocessing.
    canine_dict_SINE: dictionary where the first element of the dictionary is the pandas dataframe associated with a file of processed SINE-1 containing variants.
        The second element of the list is the numbe of variants which have multiple or no RM variants afterprocessing.
    """
    canine_dict_SINE = {}
    canine_dict_LINE = {}

    elements = ["LINE","SINE"]
    for TE in elements:
        for canine in misc_info_dict['canines']:
            print(canine)

            states = ['query','ref']
            for state in states:
                if canine == misc_info_dict['ref_name']:
                    continue
                path = f"{misc_info_dict['Directory_of_Data']}/{canine}_{misc_info_dict['ref_name']}_SV_{state}_filled_intersect_with_{TE}s_processed.bed"
                with open(path,'rt') as infile:
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
                #print(df.shape)

                #filtering options
                continuing = True

                remove_indices = []
                remove_indices_TE_extends_or_no_RM_orientation = []#a seperate list of indices which will remove loci which are low quality and have a RepeatMasker identified TE which extends beyond the boundaries of the element. Suggesting a deletion or rearragnement rather than an insertion
                for i in range(len(df["Age_TSD_Seq_left"])):
                    if TE == "LINE":#Todo. add to methods
                        if df.loc[i]["TE_extends"] == "True" and len(df.loc[i]["Age_TSD_Seq_left"]) < 10:
                            remove_indices_TE_extends_or_no_RM_orientation.append(i)
                            continue
                    if df.loc[i]["RM_orientation"] == "None":
                        remove_indices_TE_extends_or_no_RM_orientation.append(i)
                        print("I AM HERE. HERE IS WHERE I AM")
                        print(df.loc[i])
                        ##sys.exit("Illegal")
                        continue
                    if df.loc[i]['RM_orientation'] == "Multiple":#I think keeping these in is not correct. I will include them with the knowledge that it increases TSD fraction at the expense of Poly(A)
                        #print(df.loc[i])
                        remove_indices.append(i) 
                    #if df.loc[i]['No_predicted_locus'] == True:
                        #remove_indices.append(i)
                        #sys.exit()


                #print(df.shape)

                #handle processing the poly(A)s. Specifically, extract the longest poly(A) within 5bp of the corresponding TSD or end of insertion.
                #i.e. cannot be more than 5bp separating the poly(A) and TSD/end of element
                new_column = []
                for row in df.itertuples():
                    if "N/A" in row.poly_a_list:
                        new_column.append("N/A")
                        continue
                    append_value = []
                    if row.RM_orientation == "Forward":
                        if "N/A" in row.right_TSD_filled:
                            boundary = int(row.Refined_Filled_right_bound)
                        else:
                            boundary = int(row.right_TSD_filled.split("\\")[0])
                        #print(boundary)
                        for poly_a_coords in row.poly_a_list.split("\\")[::-1]:
                            poly_a_list = [int(x) for x in poly_a_coords.replace("]","").replace("[","").split(",")] #prep the coordinates to get into list
                            distance = boundary - poly_a_list[1] -1 #calculate n nucleotides separating tsd and boundary
                            if distance <= 5:
                                if append_value == []:
                                    append_value = poly_a_list
                                else:
                                    if poly_a_list[1]-poly_a_list[0] > append_value[1]- append_value[0]:
                                        append_value = poly_a_list
                                    else:
                                        continue

                    elif row.RM_orientation == "Reverse":
                        if "N/A" in row.right_TSD_filled:
                            boundary = int(row.Refined_Filled_left_bound)
                        else:
                            boundary = int(row.left_TSD_filled.split("\\")[1])
                        for poly_a_coords in row.poly_a_list.split("\\"):
                            poly_a_list = [int(x) for x in poly_a_coords.replace("]","").replace("[","").split(",")] #prep the coordinates to get into list
                            distance =poly_a_list[0] - boundary -1 #calculate n nucleotides separating tsd and boundary
                            if distance <= 5:
                                if append_value == []:
                                    append_value = poly_a_list
                                else:
                                    if poly_a_list[1]-poly_a_list[0] > append_value[1]- append_value[0]:
                                        append_value = poly_a_list
                                    else:
                                        continue
                    else:
                        sys.exit()
                    if append_value != []:
                        new_column.append(append_value)
                    else:
                        new_column.append("N/A")
                #if TE == "LINE":
                df.insert(pandas_columns, "Adjacent_Poly", new_column)
                #else:#Todo. Once I re-run the SINE this won't be necessary
                #    df.insert(pandas_columns-2, "Adjacent_Poly", new_column)
                print(df.shape)
                
                
                
                print(remove_indices_TE_extends_or_no_RM_orientation)
                filter_processed_files(df,list(set(remove_indices_TE_extends_or_no_RM_orientation)),path)
                
                #sys.exit()
                remove_indices = set(remove_indices)
                print("Len remove indices is ", len(remove_indices))
                #remove_indices_TE_extends_or_no_RM_orientation = set(remove_indices_TE_extends_or_no_RM_orientation)
                remove_indices_combined = set(list(remove_indices)+remove_indices_TE_extends_or_no_RM_orientation)
                if len(remove_indices_combined) > 0:
                    for index in sorted(remove_indices_combined, reverse=True):
                        df.drop(index, inplace=True)
                        
                print(df.shape)
                if TE == "LINE":
                    canine_dict_LINE[f"{canine}_{state}"] = [df.copy(),len(remove_indices)]
                    print(canine_dict_LINE[f"{canine}_{state}"][1])
                else:
                    canine_dict_SINE[f"{canine}_{state}"] = [df.copy(),len(remove_indices)]
                    print(canine_dict_SINE[f"{canine}_{state}"][1])
    print(canine_dict_LINE.keys(),canine_dict_SINE.keys())
    #sys.exit()
    return canine_dict_LINE,canine_dict_SINE


def parse_Hallmarks_in_dataframe(pd,df):
    print(df.shape)
    TSD,Poly_A,Both,Neither,Both_Auto,TSD_Auto = 0,0,0,0,0,0
    min_TSD_len = 10
    min_Poly_A_len = 10
    endos = []
    short_endo = []
    for row in df.itertuples():
        if len(row.Age_TSD_Seq_right) >= min_TSD_len and row.Age_TSD_Seq_right != "N/A" and row.Adjacent_Poly != "N/A" and row.Adjacent_Poly[1]-row.Adjacent_Poly[0]+1 >= min_Poly_A_len and row.RM_orientation != "Multiple": #has a TSD and a poly a
            Both +=1
            if row.chr != "chrX":
                Both_Auto +=1
        elif len(row.Age_TSD_Seq_right) >= min_TSD_len and row.Age_TSD_Seq_right != "N/A":# and row.Adjacent_Poly != "N/A":
            TSD +=1
            if row.chr != "chrX":
                TSD_Auto +=1
        elif row.Adjacent_Poly != "N/A" and row.Adjacent_Poly[1]-row.Adjacent_Poly[0]+1 >= min_Poly_A_len:
            Poly_A +=1
        else:
            Neither +=1
        if row.age_endo == "N/A":
            pass
        #elif  len(row.Age_TSD_Seq_left) >= 5 and  len(row.Age_TSD_Seq_left) < min_TSD_len:
        #    short_endo.append(row.age_endo.lower())
        elif len(row.Age_TSD_Seq_left) >= min_TSD_len:
            endos.append(row.age_endo.lower())
    return TSD,Poly_A,Both,Neither,Both_Auto,TSD_Auto,endos

def calculate_TSD_lengths(pd,df):
    TSD_lengths = []
    #true_lengths = []
    for TSD_seq in df['Age_TSD_Seq_right']:
        if TSD_seq == "N/A":
            TSD_lengths.append(0)
        elif len(TSD_seq) > 30:
            TSD_lengths.append(30)
        else:
            TSD_lengths.append(len(TSD_seq))
    return TSD_lengths
    
def calculate_poly_A_lengths(pd,df):
    Poly_lengths = []
    for poly in df['Adjacent_Poly']:
        if poly == "N/A":
            Poly_lengths.append(0)
            continue
        length = poly[1]-poly[0]+1
        if length > 40:
            Poly_lengths.append(40)
        else:
            Poly_lengths.append(length)
    return Poly_lengths
    
def make_histograms(pd,mpl,logomaker,TE_type,canine_dict,misc_info_dict):
    Categories_dict = {}
    Auto_only_dict = {}

    for key in canine_dict.keys():
        #name = misc_info_dict['canines'][i]
        df = canine_dict[key][0]
        split_key = key.split("_")
        #if key[1] == "query":
        for i in range(len(misc_info_dict['canines'])):
            if misc_info_dict['canines'][i] == split_key[0]:
                altname = misc_info_dict["canines_alt_names"][i]
        print(altname)
        TSD,Poly_A,Both,Neither,Both_Auto,TSD_Auto,endos = parse_Hallmarks_in_dataframe(pd,df)
        TSD_lengths = calculate_TSD_lengths(pd,df)
        Poly_lengths = calculate_poly_A_lengths(pd,df)
        print(len(TSD_lengths),len(Poly_lengths))
        #print(f"For the sample {altname}. The average TSD length is {sum(TSD_lengths)/len(TSD_lengths)}. The average poly(A) length is {sum(Poly_lengths)
        #generate histograms
        fig = mpl.pyplot.subplots(dpi=300)
        mpl.pyplot.hist(Poly_lengths,bins=41)
        mpl.pyplot.axvline(x=10, color='red', linestyle='dashed', linewidth=1)
        mpl.pyplot.xticks([0,5,10,15,20,25,30,35,40],labels=["0","5","10","15","20","25","30","35","40+"])
        mpl.pyplot.xlabel("Poly(A) Lengths")
        if TE_type == "LINE":
            mpl.pyplot.ylabel("LINE-1 count")
            mpl.pyplot.ylim(0,150)
        else:
            mpl.pyplot.ylabel("SINEC Count")
            mpl.pyplot.ylim(0,1200)
        #mpl.pyplot.show()
        mpl.pyplot.savefig(f"{misc_info_dict['Subdir_for_plots']}/{TE_type}_{altname}_{split_key[1]}_TSD.png",dpi=300)
        mpl.pyplot.close()
        
        fig = mpl.pyplot.subplots(dpi=300)
        mpl.pyplot.hist(TSD_lengths,bins=31)
        mpl.pyplot.axvline(x=10, color='red', linestyle='dashed', linewidth=1)
        mpl.pyplot.xticks([0,5,10,15,20,25,30],labels=["0","5","10","15","20","25","30+"])
        mpl.pyplot.xlabel("TSD Lengths")
        if TE_type == "LINE":
            mpl.pyplot.ylabel("LINE-1 count")
            mpl.pyplot.ylim(0,400)
        else:
            mpl.pyplot.ylabel("SINEC Count")
            mpl.pyplot.ylim(0,3500)
        #mpl.pyplot.show()
        mpl.pyplot.savefig(f"{misc_info_dict['Subdir_for_plots']}/{TE_type}_{altname}_{split_key[1]}_Poly_A.png",dpi=300)
        mpl.pyplot.close()

        Categories_dict[key] = [Both,TSD,Poly_A,Neither]
        Auto_only_dict[key] = [Both_Auto,TSD_Auto]
        print(TSD)
        print(max(TSD_lengths))

        print(len(endos),endos[0])
        logo = logo_maker(endos,pd,logomaker,mpl)
        mpl.pyplot.savefig(f"{misc_info_dict['Subdir_for_plots']}/{TE_type}_{altname}_{split_key[1]}_EN_CLEAVE.png",dpi=300)
        mpl.pyplot.close()
        #logo_maker(short_endo,pd,logomaker,mpl)
    print(Categories_dict,Auto_only_dict)
    #sys.exit()
    return Categories_dict, Auto_only_dict
    
    
def create_hallmark_stack_barplot(np,mpl,canine_dict,Categories_dict,TE_type,misc_info_dict):
    fractions_out = open(f"{misc_info_dict['Subdir_for_plots']}/Hallmark_fractions_{TE_type}.txt",'wt')
    fractions_out.write("\t".join(["Sample_Name","Both","TSD_only","Poly_A_only","Neither","RM_fail\n"]))
    both,poly_a_only,TSD_only,neither,RM_fail = [],[],[],[],[]
    list_of_fraction_with_both = []
    lists = [both,TSD_only,poly_a_only, neither,RM_fail]
    capital_dogs = []
    states = ['query','ref']
    for state in states:
        for i in range(len(misc_info_dict['canines'])):
            canine = misc_info_dict['canines'][i]
            
            if canine == misc_info_dict['ref_name']:
                continue
            key_name = f"{canine}_{state}"
            print(canine,key_name)
            if len(Categories_dict[key_name]) <= 4:
                Categories_dict[key_name].append(canine_dict[key_name][1])
            if state == "ref" and canine == misc_info_dict['canines'][0]: #canines[0] will be Mischka for the main dataset. 
                capital_dogs.append(misc_info_dict['canines'][-1]) #the altname for the reference is always set as last in the canine names
            elif state == "query":
                capital_dogs.append(misc_info_dict["canines_alt_names"][i])
            if state == "query" or (state == "ref" and canine == misc_info_dict['canines'][0]):
                print('here')
                for j in range(len(Categories_dict[key_name])):
                    lists[j].append(Categories_dict[key_name][j])
                    print(Categories_dict[key_name][j])
            #aggregate and print the data in text format
            bth = str(round(Categories_dict[key_name][0]/sum(Categories_dict[key_name])*100,2))
            tsd_nly = str(round(Categories_dict[key_name][1]/sum(Categories_dict[key_name])*100,2))
            poly_nly = str(round(Categories_dict[key_name][2]/sum(Categories_dict[key_name])*100,2))
            nthr = str(round(Categories_dict[key_name][3]/sum(Categories_dict[key_name])*100,2))
            rm_fail = str(round(Categories_dict[key_name][4]/sum(Categories_dict[key_name])*100,2))
            
            print(capital_dogs,i,misc_info_dict['canines_alt_names'])
            out_info = [f"{capital_dogs[i]}_{state}",bth,tsd_nly,poly_nly,nthr,rm_fail]
            out_info = "\t".join(out_info) + "\n"
            fractions_out.write(out_info)
            
            list_of_fraction_with_both.append(Categories_dict[key_name][0]/sum(Categories_dict[key_name])*100)
    print(lists)
    N = len(misc_info_dict['canines'])
    ind = np.arange(N)   
    width = 0.35 

    bars_1 = np.add(both, TSD_only).tolist()
    bars_2 = np.add(np.add(both, TSD_only), poly_a_only).tolist()
    bars_3 = np.add(np.add(np.add(both, TSD_only), poly_a_only), neither).tolist()

    fig = mpl.pyplot.subplots(figsize =(20, 10),dpi=300)
    p1 = mpl.pyplot.bar(ind, both, width, color = "tab:red")
    p2 = mpl.pyplot.bar(ind, TSD_only, width, bottom = both, color = "tab:blue")
    p3 = mpl.pyplot.bar(ind, poly_a_only, width, bottom = bars_1, color = "tab:olive")
    p4 = mpl.pyplot.bar(ind, neither, width, bottom = bars_2, color = "tab:pink")
    p5 = mpl.pyplot.bar(ind, RM_fail,width, bottom = bars_3, color='tab:grey')

    mpl.pyplot.ylabel(f'{TE_type} Count', size = 20)
    #mpl.pyplot.title('Retrocopy Count by Detected Hallmarks', size = 20)
    
    mpl.pyplot.xticks(ind, capital_dogs, size = 20)
    mpl.pyplot.yticks(size=15)
    mpl.pyplot.legend((p1[0], p2[0], p3[0], p4[0], p5[0]), ('Both', 'TSD only','Poly(A) only', 'Neither','Invalid RM'), bbox_to_anchor=(1.02, 1), fontsize=15)
    
    mpl.pyplot.savefig(f"{misc_info_dict['Subdir_for_plots']}/Stacked_Barplot_{TE_type}.png",dpi=300,bbox_inches='tight')
    mpl.pyplot.close()
    fractions_out.close()
    print(f"Average % both {sum(list_of_fraction_with_both)/len(list_of_fraction_with_both)}")
    
    return
    
    

def process_Rate_document(content,canine,reference_name):
    output_dict = {}
    column_names = content[0].split()
    values = content[1].split()
    for i in range(len(column_names)):
        column_names[i] = column_names[i].replace(canine,"query")
        column_names[i] = column_names[i].replace(reference_name,'ref')
        output_dict[column_names[i]]=values[i]
    return output_dict    
    
def chrX_counter(input_file_name):
    chrX_count = 0
    with open(input_file_name,'rt') as infile:
        for line in infile:
            if line[0:4] == "chrX":
                chrX_count +=1
    return chrX_count
    
def produce_rate_estimates(pd,np,mpl,misc_info_dict,canine_dict,auto_dict,TE_type):
    #This permissive rate estimate includes the RM fail because they're real loci. As a result, I have to use the processed_filter files.
    #The restrictive rate estimate cannot include RM fail because by definition they cannot have a poly(A) and a TSD. Only a TSD
    
    comparison_sample = misc_info_dict['canines'][0]#mischka for all vs mcan comparison
    SNP_rates_dict = {"upper_SNP_bound":7.1e-9, "lower_SNP_bound":2.6e-9, "point_SNP_estimate":4.5e-9} #From E.Koch 2019. PMID: 31297530
    #print(TE_type)    
    dict_of_canines = {}
    for canine in misc_info_dict['canines']:
        
        #Get SNPs, auto_length. All on the autosomes
        if canine != misc_info_dict['ref_name']:
            with open(f"{misc_info_dict['Directory_of_Data']}/{canine}_against_{misc_info_dict['ref_name']}_Rate_Information.txt",'rt') as infile:
                info_dict = process_Rate_document(infile.readlines(),canine,misc_info_dict['ref_name'])
                if canine == comparison_sample:
                    comp_info_dict = info_dict.copy()
        else:
            continue
        #extract permissive variant count.
        states = ["query","ref"]
        for state in states:
            x_count,auto_count = 0,0#handle getting the extra graphs here. Maybe some side by side charts of the point estimate of rates.
            with open(f"{misc_info_dict['Directory_of_Data']}/{canine}_{misc_info_dict['ref_name']}_SV_{state}_filled_intersect_with_{TE_type}s_processed_filtered.bed",'rt') as infile:
                for line in infile:
                    if line.startswith("#"):
                        continue
                    if line.startswith("chrX"):#assumes that only autosomes and x are in file.
                        x_count +=1
                    else:
                        auto_count +=1
                #print(x_count,auto_count)
                info_dict[f"x_count_permissive_{state}_{TE_type}"] = x_count
                info_dict[f"auto_count_permissive_{state}_{TE_type}"] = auto_count
            
        
                        
            info_dict[f"{TE_type}_auto_count_hc_{state}"] = sum(auto_dict[f"{canine}_{state}"]) #handle the high confidence loci which have a TSD and were not filtered out in any previous step. 
            for SNP_key in SNP_rates_dict.keys():#estimate the rate of divergence 
                if f"{SNP_key}_genome_divergence_in_generations" not in info_dict.keys():
                    info_dict[f"{SNP_key}_genome_divergence_in_generations"] = round((float(info_dict['SNP_total'])*.5/float(info_dict['auto_length']))/float(SNP_rates_dict[SNP_key]),0)
                    #print(info_dict)
                #print(info_dict)
                info_dict[f"{TE_type}_{state}_inverse_rate_{SNP_key}"] = round(float(info_dict[f"{SNP_key}_genome_divergence_in_generations"])/(float(info_dict[f"auto_count_permissive_{state}_{TE_type}"])),3)
                info_dict[f"{TE_type}_{state}_inverse_rate_hc_{SNP_key}"] = round(float(info_dict[f"{SNP_key}_genome_divergence_in_generations"])/(float(info_dict[f"{TE_type}_auto_count_hc_{state}"])),3)
            print('state',state,info_dict)
        dict_of_canines[canine] = info_dict
        print(dict_of_canines)


    #calculate eror for genome divergence. 
    #NOTE: a lower SNP rate means that the number of generations since divergence increases in this model.
    SNP_count = []
    upper_SNP_bound_error_bars_gen_div = [] #these bars go below the point estimate in the chart
    lower_SNP_bound_error_bars_gen_dive = [] #these bars go above the point estimate in the chart
    gen_div = [] #point estimate for generations since divergence
    for canine in misc_info_dict['canines']:        
        if canine == misc_info_dict['ref_name']:
            continue
        SNP_count.append(int(dict_of_canines[canine]["SNP_total"])/1000000)
        upper_SNP_bound_error_bars_gen_div.append(dict_of_canines[canine]["point_SNP_estimate_genome_divergence_in_generations"]-dict_of_canines[canine]["upper_SNP_bound_genome_divergence_in_generations"])
        lower_SNP_bound_error_bars_gen_dive.append(dict_of_canines[canine]["lower_SNP_bound_genome_divergence_in_generations"]-dict_of_canines[canine]["point_SNP_estimate_genome_divergence_in_generations"])
        gen_div.append(dict_of_canines[canine]["point_SNP_estimate_genome_divergence_in_generations"])
    
   # sys.exit() #i'm here i think

    #calculate error between point estimate and upper/lower SNP bounds. NOTE: a lower SNP rate means a higher TE rate in this model because the 2 are anti-correlated.
    #NOTE: The retrotransposition rates are plotted in Births per new TE. This means that the inverse rate being 20 means that one out of 20 births are expectedto have a new TE.
    #Thus, if the snp rate decreases, the TE rate increases, which causes the birhts per TE to go down. 
    
    upper_SNP_bound_error_bars_permissive_rate = [] 
    lower_SNP_bound_error_bars_permissive_rate = []
    inverse_TE_rates = [] #permissive rate list
    
    upper_SNP_bound_error_bars_hc = [] 
    lower_SNP_bound_error_bars_hc = []
    inverse_TE_rates_hc = [] #hc rate list
    
    for canine in misc_info_dict['canines']:
        if canine != misc_info_dict['ref_name']:
            state = "query"
        else:
            canine = misc_info_dict['canines'][0]
            state = "ref"
        upper_SNP_bound_error_bars_permissive_rate.append(dict_of_canines[canine][f"{TE_type}_{state}_inverse_rate_point_SNP_estimate"] - dict_of_canines[canine][f"{TE_type}_{state}_inverse_rate_upper_SNP_bound"])
        lower_SNP_bound_error_bars_permissive_rate.append(dict_of_canines[canine][f"{TE_type}_{state}_inverse_rate_lower_SNP_bound"] - dict_of_canines[canine][f"{TE_type}_{state}_inverse_rate_point_SNP_estimate"])
        inverse_TE_rates.append(dict_of_canines[canine][f"{TE_type}_{state}_inverse_rate_point_SNP_estimate"])
        
        upper_SNP_bound_error_bars_hc.append(dict_of_canines[canine][f"{TE_type}_{state}_inverse_rate_hc_point_SNP_estimate"] - dict_of_canines[canine][f"{TE_type}_{state}_inverse_rate_hc_upper_SNP_bound"])
        lower_SNP_bound_error_bars_hc.append(dict_of_canines[canine][f"{TE_type}_{state}_inverse_rate_hc_lower_SNP_bound"] - dict_of_canines[canine][f"{TE_type}_{state}_inverse_rate_hc_point_SNP_estimate"])
        inverse_TE_rates_hc.append(dict_of_canines[canine][f"{TE_type}_{state}_inverse_rate_hc_point_SNP_estimate"])
    
    #plot_and_save SNP count figure
    fig, ax = mpl.pyplot.subplots()
    mpl.pyplot.bar(misc_info_dict['canines_alt_names'][0:-1],SNP_count,color='tan')
    mpl.pyplot.ylabel("SNV Count (Millions)")
    mpl.pyplot.savefig(f"{misc_info_dict['Subdir_for_plots']}/SNV_count_to_{misc_info_dict['canines_alt_names'][-1]}_{TE_type}.png",dpi=300,bbox_inches="tight")
    mpl.pyplot.close() 
    
    #plot and save gen_div figure
    gen_div_errors = [upper_SNP_bound_error_bars_gen_div,lower_SNP_bound_error_bars_gen_dive]
    fig, ax = mpl.pyplot.subplots()
    mpl.pyplot.bar(misc_info_dict['canines_alt_names'][0:-1],gen_div,color='lightblue',yerr=gen_div_errors,capsize=3)
    mpl.pyplot.ylabel("Generations Since Genome Divergence")
    mpl.pyplot.savefig(f"{misc_info_dict['Subdir_for_plots']}/Generations_to_{misc_info_dict['canines_alt_names'][-1]}_{TE_type}.png",dpi=300,bbox_inches="tight")
    mpl.pyplot.close() 
   
    
   #plot and save permissive 
    
    perm_errors = [upper_SNP_bound_error_bars_permissive_rate,lower_SNP_bound_error_bars_permissive_rate]
    fig, ax = mpl.pyplot.subplots()
    if TE_type == "LINE":
        mpl.pyplot.ylabel(f"{TE_type}-1 rate (1/n live births)")
        block_color = "lightgreen"
    else:
        mpl.pyplot.ylabel(f"{TE_type} rate (1/n live births)")
        block_color = "purple"
    mpl.pyplot.axhline(y=sum(inverse_TE_rates)/len(inverse_TE_rates), color='r', linestyle='--')
    mpl.pyplot.bar(misc_info_dict['canines_alt_names'],inverse_TE_rates,color=block_color,yerr=perm_errors,capsize=3)

    #mpl.pyplot.ylim(0,50)
    mpl.pyplot.savefig(f"{misc_info_dict['Subdir_for_plots']}/{TE_type}_rate_to_{misc_info_dict['canines_alt_names'][-1]}.png",dpi=300,bbox_inches="tight")
    mpl.pyplot.close() 
    
   #plot and save hc 
    hc_errors = [upper_SNP_bound_error_bars_hc,lower_SNP_bound_error_bars_hc]
    fig, ax = mpl.pyplot.subplots()
    mpl.pyplot.axhline(y=sum(inverse_TE_rates_hc)/len(inverse_TE_rates_hc), color='r', linestyle='--')
    mpl.pyplot.bar(misc_info_dict['canines_alt_names'],inverse_TE_rates_hc,color=block_color,yerr=hc_errors,capsize=3)
    if TE_type == "LINE":
        mpl.pyplot.ylabel(f"{TE_type}-1 refined rate (1/n live births)")
    else:
        mpl.pyplot.ylabel(f"{TE_type} refined rate (1/n live births)")
    #mpl.pyplot.ylim(0,50)
    mpl.pyplot.savefig(f"{misc_info_dict['Subdir_for_plots']}/{TE_type}_rate_to_{misc_info_dict['canines_alt_names'][-1]}_refined.png",dpi=300,bbox_inches="tight")
    mpl.pyplot.close() 
    
    print(f"SNP count, generations since divergence, permissive rate, and high confidence rate for {TE_type}s")
    
    
    
    
    #create side by side bar plots for count fo variants (total. This will include variants that have multiple orientations) and side by side point estimates of the rat
    TE_count_ref_auto = []
    TE_count_sample_auto = []
    TE_count_ref_X = []
    TE_count_sample_X = []
    rate_ref = []
    rate_sample = []
    
    for canine in misc_info_dict['canines'][0:-1]:
        TE_count_ref_auto.append(dict_of_canines[canine][f"auto_count_permissive_ref_{TE_type}"])
        TE_count_sample_auto.append(dict_of_canines[canine][f"auto_count_permissive_query_{TE_type}"])
        TE_count_ref_X.append(dict_of_canines[canine][f"x_count_permissive_ref_{TE_type}"])
        TE_count_sample_X.append(dict_of_canines[canine][f"x_count_permissive_query_{TE_type}"])
        rate_ref.append(dict_of_canines[canine][f"{TE_type}_ref_inverse_rate_point_SNP_estimate"])
        rate_sample.append(dict_of_canines[canine][f"{TE_type}_query_inverse_rate_point_SNP_estimate"])

    width = 0.25
    n=len(misc_info_dict['canines'])-1
    r = np.arange(n) 
    if TE_type == "LINE":
        ylabel_mod = "LINE-1"
    else:
        ylabel_mod = "SINEC"
        
    mpl.pyplot.bar(r, TE_count_ref_auto,color = 'b', 
        width = width, edgecolor = 'black', 
        label=misc_info_dict['canines_alt_names'][-1]) 
    mpl.pyplot.bar(r+width,TE_count_sample_auto, color = 'g', 
        width = width, edgecolor = 'black', 
        label='Sample') 

    mpl.pyplot.ylabel(f"Dimorphic Autosomal {ylabel_mod} Count")  
    mpl.pyplot.xticks(r + width/2,misc_info_dict['canines_alt_names'][0:-1]) 
    mpl.pyplot.legend(bbox_to_anchor=(0.0, -.14, 1., .102), ncol=2,loc='lower center', frameon=False, borderaxespad=0,fontsize='8')
    #print((sum(TE_count_ref_auto)+sum(TE_count_sample_auto))/len(new_canine_list))
    mpl.pyplot.savefig(f"{misc_info_dict['Subdir_for_plots']}/{TE_type}_auto_count.png",dpi=300,bbox_inches="tight")
    mpl.pyplot.close() 
    
    mpl.pyplot.bar(r, TE_count_ref_X,color = 'b', 
        width = width, edgecolor = 'black', 
        label=misc_info_dict['canines_alt_names'][-1]) 
    mpl.pyplot.bar(r+width,TE_count_sample_X, color = 'g', 
        width = width, edgecolor = 'black', 
        label='Sample') 

    mpl.pyplot.ylabel(f"Dimorphic X Chromosome {ylabel_mod} Count")  
    mpl.pyplot.xticks(r + width/2,misc_info_dict['canines_alt_names'][0:-1]) 
    mpl.pyplot.legend(bbox_to_anchor=(0.0, -.14, 1., .102), ncol=2,loc='lower center', frameon=False, borderaxespad=0,fontsize='8')
    mpl.pyplot.savefig(f"{misc_info_dict['Subdir_for_plots']}/{TE_type}_X_count.png",dpi=300,bbox_inches="tight")
    mpl.pyplot.close() 
    
    
    mpl.pyplot.bar(r, rate_ref,color = 'b', 
        width = width, edgecolor = 'black', 
        label=misc_info_dict['canines_alt_names'][-1]) 
    mpl.pyplot.bar(r+width,rate_sample, color = 'g', 
        width = width, edgecolor = 'black', 
        label='Sample') 

    mpl.pyplot.ylabel(f"{ylabel_mod} Rate (1/n live births)")  
    mpl.pyplot.xticks(r + width/2,misc_info_dict['canines_alt_names'][0:-1]) 
    mpl.pyplot.legend(bbox_to_anchor=(0.0, -.14, 1., .102), ncol=2,loc='lower center', frameon=False, borderaxespad=0,fontsize='8')
    mpl.pyplot.savefig(f"{misc_info_dict['Subdir_for_plots']}/{TE_type}_Rate_Side_by_Side.png",dpi=300,bbox_inches="tight")
    mpl.pyplot.close() 
    
    
    #create table which can be pruned as output.
    header = ""
    body = []
    for i in range(len(misc_info_dict['canines'])-1):
        canine = misc_info_dict['canines'][i]
        if header == "":
            header = list(dict_of_canines[canine].keys())
            header.insert(0,"Sample_Name")
        body.append(list(dict_of_canines[canine].values()))
        body[-1].insert(0,misc_info_dict['canines_alt_names'][i])
    print(header,body)
    with open(f"{misc_info_dict['Subdir_for_plots']}/Rate_information_post_processing_{TE_type}.txt",'wt') as out:
        out.write("\t".join(header))
        out.write("\n")
        for output in body:
            output = [str(x) for x in output]
            output = "\t".join(output) + "\n"
            out.write(output)
    print(f"Output table for Rate information post processing {TE_type}s")
    #dict_of_canines[canine]
    #print(dict_of_canines)
    
    
def scatterplot_tsd_vs_polya(pd,np,scipy,mpl,misc_info_dict,canine_dict,TE_type):
    for key in canine_dict.keys():
        scatter_TSD_lengths = []
        scatter_Poly_lengths = []
        #print(key,misc_info_dict['canines'][0],misc_info_dict['ref_name'])
        #sys.exit()
        if key.split("_")[0] not in [misc_info_dict['canines'][0]]: #For all v mcan this is Mischka query and ref
            continue
        #print(key)
        df = canine_dict[key][0]
        for TSD_seq in df.Age_TSD_Seq_right:
            if TSD_seq == "N/A":
                scatter_TSD_lengths.append(0)
            else:
                scatter_TSD_lengths.append(len(TSD_seq))

        for poly in df['Adjacent_Poly']:
            if poly == "N/A":
                scatter_Poly_lengths.append(0)
                continue
            length = poly[1]-poly[0]+1
            scatter_Poly_lengths.append(length)
        x = scipy.stats.pearsonr(scatter_TSD_lengths,scatter_Poly_lengths)
        print(x)


        # Generate some test data
        y = scatter_TSD_lengths
        x = scatter_Poly_lengths
        
        print(f"Looking for outliers in {key}")
        for i in range(len(y)):
            if y[i] > 150 or x[i] > 90:
                print(x[i],y[i])
        #[print(z) for z in x if z > 100]
        #print(max(x),max(y))
        fig,ax = mpl.pyplot.subplots(dpi=300)
        heatmap, xedges, yedges = np.histogram2d(x, y, bins=(max(x)+1,max(y)+1))
        extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
        mpl.pyplot.imshow(heatmap.T, norm = mpl.colors.LogNorm(), extent=extent, origin='lower')
        mpl.pyplot.ylabel("TSD lengths")
        mpl.pyplot.xlabel("poly(A) length")
        mpl.pyplot.xlim(0,90)
        mpl.pyplot.ylim(0,150)
        cbar = mpl.pyplot.colorbar(ticks=[1, 10, 100])
        cbar.ax.set_yticklabels(['1', '10', '100'])  # vertically oriented colorbar
        mpl.pyplot.show()
        if TE_type == "LINE":
            TE_name = "LINE1"
        else:
            TE_name = "SINEC"
            
        mpl.pyplot.savefig(f"{misc_info_dict['Subdir_for_plots']}/{key}_{TE_name}_Heatmap_Hallmarks.png",dpi=300,bbox_inches="tight")
        mpl.pyplot.close() 
        #there is an outlier here, but I don't think it is worth keeping in because it distorts the figure
        #sys.exit()

def subfamily_analysis(pd,mpl,canine_dict,misc_info_dict,TE_type):
    dict_cf = {}
    min_TSD_len = 10
    for canine in canine_dict.keys():
        print(canine)
        df = canine_dict[canine][0]
        
        
        dict_keys = {}
        unresolved = 0
        
        dict_keys_HC = {}
        unresolved_HC = 0 
        for row in df.itertuples():
            #print(RM)
            #extract subfamily information
            RM = RM = row.RM_subfamilies
            if TE_type == "SINE":
                replace_name = "SINEC2A1_CF"
                RM = RM.replace(replace_name, "SINEC_Cf").split("\\")
            else:
                RM = RM.split("\\")
            #print(RM)
            
            
            #check if the locus has a high confidence TSD:
            if len(row.Age_TSD_Seq_right) >= min_TSD_len and row.Age_TSD_Seq_right != "N/A":
                long_TSD = True
            else:
                long_TSD = False

            if len(set(RM)) == 1:
                if RM[0] not in dict_keys.keys():
                    dict_keys[RM[0]] = 1
                else:
                    dict_keys[RM[0]] +=1
                if long_TSD == True:
                    if RM[0] not in dict_keys_HC.keys():
                        dict_keys_HC[RM[0]] = 1
                    else:
                        dict_keys_HC[RM[0]] +=1
            else:
                unresolved +=1
                if long_TSD == True:
                    unresolved_HC +=1
                

        dict_of_dicts = [dict_keys,dict_keys_HC]   
        for i in range(len(dict_of_dicts)):
            vals = []
            dictionary = dict_of_dicts[i]
            if TE_type == "SINE":
                dictionary["SINEC_Cf/2A1"] = dictionary.pop("SINEC_Cf")
            sorted_dict = dict(sorted(dictionary.items(), key=lambda item: item[1],reverse=True))
            if i == 0:
                if unresolved != 0:
                    sorted_dict['unresolved'] = unresolved
                    dictionary['unresolved'] = unresolved
            elif i == 1:
                if unresolved_HC != 0:
                    sorted_dict['unresolved'] = unresolved_HC
                    dictionary['unresolved'] = unresolved_HC
            keys_list = sorted_dict.keys()
            for key in keys_list:
                vals.append(dictionary[key])
            print(len(keys_list),len(vals))

            fig = mpl.pyplot.figure(figsize=(12,9),dpi=300)
            plot = mpl.pyplot.bar(keys_list,vals,label=vals)
            mpl.pyplot.bar_label(plot, label_type='edge', fontsize=12)
            mpl.pyplot.xticks(rotation=45, ha='right')
            if TE_type == "SINE":
                mpl.pyplot.xlabel("SINEC subfamily")
                if i == 0:
                    mpl.pyplot.ylabel("Subfamily Dimorphic SINE count")
                    mpl.pyplot.savefig(f"{misc_info_dict['Subdir_for_plots']}/{canine}_{TE_type}C_Subfamily_Fractions.png",dpi=300,bbox_inches="tight")
                elif i ==1:
                    mpl.pyplot.ylabel("Subfamily Dimorphic SINE count (High Confidence TSD Only)")
                    mpl.pyplot.savefig(f"{misc_info_dict['Subdir_for_plots']}/{canine}_{TE_type}C_Subfamily_Fractions_HC.png",dpi=300,bbox_inches="tight")
            else:
                mpl.pyplot.xlabel("LINE-1 subfamily")
                if i == 0:
                    mpl.pyplot.ylabel("Subfamily Dimorphic LINE-1 count")
                    mpl.pyplot.savefig(f"{misc_info_dict['Subdir_for_plots']}/{canine}_{TE_type}1_Subfamily_Fractions.png",dpi=300,bbox_inches="tight")
                elif i ==1:
                    mpl.pyplot.ylabel("Subfamily Dimorphic LINE-1 count (High Confidence TSD Only)")
                    mpl.pyplot.savefig(f"{misc_info_dict['Subdir_for_plots']}/{canine}_{TE_type}1_Subfamily_Fractions_HC.png",dpi=300,bbox_inches="tight")

                
            mpl.pyplot.close() 
            print(keys_list)
            print(vals)
            for i in range(len(keys_list)):  
                if vals[i] > 25:
                    print(canine,list(keys_list)[i],vals[i],round(vals[i]/sum(vals),4)*100)