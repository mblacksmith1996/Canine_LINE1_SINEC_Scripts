import sys
import os
import subprocess

def complement(c):
    if c == 'A':
        return 'T'
    if c == 'T':
        return 'A'
    if c == 'C':
        return 'G'
    if c == 'G':
        return 'C'
    if c == 'a':
        return 't'
    if c == 't':
        return 'a'
    if c == 'c':
        return 'g'
    if c == 'g':
        return 'c'
    # If not ACTGactg simply return same character
    return c    

def revcomp(seq):
    c = ''
    seq = seq[::-1] #reverse
    # Note, this good be greatly sped up using list operations
    seq = [complement(i) for i in seq]
    c = ''.join(seq)
    return c
    
    
def get_info_from_subprocess(cmd):
    proc = subprocess.Popen(cmd,shell=True,stdout=subprocess.PIPE)
    out,err = proc.communicate()
    out = out.decode()
    return out    

def extract_sequence_from_faidx_cmd(cmd):
    out = get_info_from_subprocess(cmd)
    final_seq = ""
    for seq in out.split("\n")[1::]:
        final_seq = final_seq + seq
    return final_seq

def process_age(age_out,locus_information):
    #function originally copied from Checking_Dimorphic_Retrocopies.py
    print(age_out)
    age_out = age_out.split("\n")
    if len(age_out) == 4:
        locus_information["Refined_Filled_start"] = "N/A"
        locus_information["Refined_Filled_end"] = "N/A"
    for i in range(len(age_out)):
        if age_out[i].startswith("Alignment:"):
            #print(age_out[i+1],age_out[i+2])
            #print((age_out[i+1].split()))
            age_out[i+1] = age_out[i+1].replace("]","").replace("[","").replace(","," ")
            age_out[i+2] = age_out[i+2].replace("]","").replace("[","").replace(","," ")
            
            if "EXCISED" not in age_out[i+1]:
                locus_information["Refined_Filled_start"] = "N/A"
                locus_information["Refined_Filled_end"] = "N/A"
                locus_information["Refined_Filled_left_bound"] = "N/A"
                locus_information["Refined_Filled_right_bound"] = "N/A"
                
                locus_information["Refined_Empty_start"] = "N/A"
                locus_information["Refined_Empty_end"] = "N/A"
                locus_information["Refined_Empty_left_bound"] = "N/A"
                locus_information["Refined_Empty_right_bound"] = "N/A"
                locus_information["Age_Insertion_Found"] = False
            else:
                locus_information["Refined_Filled_left_bound"] = int(age_out[i+1].split()[4])+locus_information['Filled_start']-locus_information['age_flank']
                locus_information["Refined_Filled_right_bound"] = int(age_out[i+1].split()[7])+locus_information['Filled_start']-locus_information['age_flank']
                locus_information["Age_Insertion_Found"] = True
                if abs(locus_information["Refined_Filled_left_bound"]-locus_information["Refined_Filled_right_bound"]) == 1:
                    locus_information["Refined_Filled_start"] = "N/A"
                    locus_information["Refined_Filled_end"] = "N/A"

                else:
                    locus_information["Refined_Filled_start"] = int(age_out[i+1].split()[4])+locus_information['Filled_start']-locus_information['age_flank']+1
                    locus_information["Refined_Filled_end"] = int(age_out[i+1].split()[7])+locus_information['Filled_start']-locus_information['age_flank']-1
                
                
                #if locus_information["Reference_orientation"] == "-":
                #    start_modifier = 0
                #    end_modifier = 0
                    #Note I think age gets me most of the way there on its own. It uses original coordinates.
                locus_information["Refined_Empty_left_bound"] = int(age_out[i+2].split()[4])+locus_information['Empty_start']-locus_information['age_flank']
                locus_information["Refined_Empty_right_bound"] = int(age_out[i+2].split()[7])+locus_information['Empty_start']-locus_information['age_flank']

                if abs(locus_information["Refined_Empty_left_bound"]- locus_information["Refined_Empty_right_bound"]) == 1:
                    print("No corresponding deletion in empty site")
                    #print(locus_information)
                    locus_information["Refined_Empty_start"] = "N/A"
                    locus_information["Refined_Empty_end"] = "N/A"
                #sys.exit()
                else:
                    print("Deletion found in empty site")
                    
                    #new version. Always have refined empty left bound <= refined empty right bound.
                    locus_information['Refined_Empty_start'] = min(locus_information["Refined_Empty_left_bound"],locus_information["Refined_Empty_right_bound"]) + 1
                    locus_information['Refined_Empty_end'] = max(locus_information["Refined_Empty_left_bound"],locus_information["Refined_Empty_right_bound"]) -1
            print(locus_information)
            #sys.exit()
        elif age_out[i].startswith("Identity at breakpoints:"):
            #print(age_out[i+1],age_out[i+2])
            age_out[i+1] = age_out[i+1].replace("]","").replace("[","")
            age_out[i+2] = age_out[i+2].replace("]","").replace("[","")
            if age_out[i+1].split()[3] != "0":
                len_tsd = int(age_out[i+1].split()[3])
                #print(len_tsd)
                #handle
                locus_information['left_TSD_filled'] = [str(int(age_out[i+1].split()[5].split(",")[0])+locus_information['Filled_start']-locus_information['age_flank']),str(int(age_out[i+1].split()[5].split(",")[1])+locus_information['Filled_start']-locus_information['age_flank'])]
                locus_information['right_TSD_filled'] = [str(int(age_out[i+1].split()[7].split(",")[0])+locus_information['Filled_start']-locus_information['age_flank']),str(int(age_out[i+1].split()[7].split(",")[1])+locus_information['Filled_start']-locus_information['age_flank'])]
            
            
                #extract left TSD, right tsd, and empty site
                left_cmd = f"samtools faidx {locus_information['Ref']} {locus_information['chr']}:{locus_information['left_TSD_filled'][0]}-{locus_information['left_TSD_filled'][1]}"
                #print(left_cmd)
                left_cmd_out = extract_sequence_from_faidx_cmd(left_cmd).lower()
                #print(left_cmd_out)
                right_cmd = f"samtools faidx {locus_information['Ref']} {locus_information['chr']}:{locus_information['right_TSD_filled'][0]}-{locus_information['right_TSD_filled'][1]}"
                right_cmd_out = extract_sequence_from_faidx_cmd(right_cmd).lower()
                #print(right_cmd)
                
                #Note: Requires same chrom for both query and reference                   
                if locus_information["Reference_orientation"] == "-":
                    if int(locus_information["left_TSD_filled"][0])<= locus_information["Refined_Filled_left_bound"] and int(locus_information["left_TSD_filled"][1]) >= locus_information["Refined_Filled_start"]:
                        #print(locus_information["left_TSD_filled"],locus_information["Refined_Filled_left_bound"])
                        upstream = locus_information["Refined_Filled_left_bound"]-int(locus_information["left_TSD_filled"][0])+1
                        downstream = int(locus_information['left_TSD_filled'][1])-locus_information['Refined_Filled_start']+1
                        print(upstream,downstream)
                        empty_site_cmd1 = f"samtools faidx {locus_information['Query']} -i {locus_information['chr']}:{int(locus_information['Refined_Empty_left_bound'])}-{int(locus_information['Refined_Empty_left_bound'])+upstream-1}"
                        empty_site_cmd2 = f"samtools faidx {locus_information['Query']} -i {locus_information['chr']}:{int(locus_information['Refined_Empty_right_bound'])-downstream+1}-{int(locus_information['Refined_Empty_right_bound'])}"
                        print(empty_site_cmd1)
                        print(empty_site_cmd2)
                        locus_information['empty_coords'] = f"{int(locus_information['Refined_Empty_left_bound'])}-{int(locus_information['Refined_Empty_left_bound'])+upstream-1}\\{int(locus_information['Refined_Empty_right_bound'])-downstream+1}-{int(locus_information['Refined_Empty_right_bound'])}"
                        empty_site_cmd_out1 = extract_sequence_from_faidx_cmd(empty_site_cmd1).lower()
                        empty_site_cmd_out2 = extract_sequence_from_faidx_cmd(empty_site_cmd2).lower()
                        empty_site_cmd_out = empty_site_cmd_out1+empty_site_cmd_out2
                        #print(empty_site_cmd_out)
                    elif int(locus_information["left_TSD_filled"][1]) <= locus_information["Refined_Filled_left_bound"]:#check if upper bound of left TSD is left of the insertion sequence. 
                        empty_site_cmd = f"samtools faidx {locus_information['Query']} -i {locus_information['chr']}:{int(locus_information['Refined_Empty_left_bound'])}-{int(locus_information['Refined_Empty_left_bound'])+len_tsd-1}"
                        empty_site_cmd_out = extract_sequence_from_faidx_cmd(empty_site_cmd).lower()
                        locus_information['empty_coords'] = f"{int(locus_information['Refined_Empty_left_bound'])}-{int(locus_information['Refined_Empty_left_bound'])+len_tsd-1}"
                        #print(empty_site_cmd)
                    else:
                        empty_site_cmd = f"samtools faidx {locus_information['Query']} -i {locus_information['chr']}:{int(locus_information['Refined_Empty_right_bound'])-(len_tsd)+1}-{int(locus_information['Refined_Empty_right_bound'])}"
                        empty_site_cmd_out = extract_sequence_from_faidx_cmd(empty_site_cmd).lower()
                        locus_information['empty_coords'] = f"{int(locus_information['Refined_Empty_right_bound'])-(len_tsd)+1}-{int(locus_information['Refined_Empty_right_bound'])}"
                        #print(empty_site_cmd)
                else:
                    if int(locus_information["left_TSD_filled"][0])<= locus_information["Refined_Filled_left_bound"] and int(locus_information["left_TSD_filled"][1]) >= locus_information["Refined_Filled_start"]:
                        #print(locus_information["left_TSD_filled"],locus_information["Refined_Filled_left_bound"])
                        upstream = locus_information["Refined_Filled_left_bound"]-int(locus_information["left_TSD_filled"][0])+1
                        downstream = int(locus_information['left_TSD_filled'][1])-locus_information['Refined_Filled_start']+1
                        print(upstream,downstream)
                        empty_site_cmd1 = f"samtools faidx {locus_information['Query']} {locus_information['chr']}:{int(locus_information['Refined_Empty_left_bound'])-upstream+1}-{int(locus_information['Refined_Empty_left_bound'])}"
                        empty_site_cmd2 = f"samtools faidx {locus_information['Query']} {locus_information['chr']}:{int(locus_information['Refined_Empty_right_bound'])}-{int(locus_information['Refined_Empty_right_bound'])+downstream-1}"
                        locus_information['empty_coords'] = f"{int(locus_information['Refined_Empty_left_bound'])-upstream+1}-{int(locus_information['Refined_Empty_left_bound'])}\\{int(locus_information['Refined_Empty_right_bound'])}-{int(locus_information['Refined_Empty_right_bound'])+downstream-1}"
                        print(empty_site_cmd1)
                        print(empty_site_cmd2)
                        empty_site_cmd_out1 = extract_sequence_from_faidx_cmd(empty_site_cmd1).lower()
                        empty_site_cmd_out2 = extract_sequence_from_faidx_cmd(empty_site_cmd2).lower()
                        empty_site_cmd_out = empty_site_cmd_out1+empty_site_cmd_out2
                        #print(empty_site_cmd_out)
                        print(empty_site_cmd_out)
                        #sys.exit("i am here")
                    elif int(locus_information["left_TSD_filled"][1]) <= locus_information["Refined_Filled_left_bound"]:#check if upper bound of left TSD is left of the insertion sequence. 
                        empty_site_cmd = f"samtools faidx {locus_information['Query']} {locus_information['chr']}:{int(locus_information['Refined_Empty_left_bound'])-len_tsd+1}-{int(locus_information['Refined_Empty_left_bound'])}"
                        empty_site_cmd_out = extract_sequence_from_faidx_cmd(empty_site_cmd).lower()
                        locus_information['empty_coords'] = f"{int(locus_information['Refined_Empty_left_bound'])-len_tsd+1}-{int(locus_information['Refined_Empty_left_bound'])}"
                        print(empty_site_cmd_out)
                        #print(empty_site_cmd)
                    else:
                        empty_site_cmd = f"samtools faidx {locus_information['Query']} {locus_information['chr']}:{int(locus_information['Refined_Empty_right_bound'])}-{int(locus_information['Refined_Empty_right_bound'])+len_tsd-1}"
                        empty_site_cmd_out = extract_sequence_from_faidx_cmd(empty_site_cmd).lower()
                        locus_information['empty_coords'] = f"{int(locus_information['Refined_Empty_right_bound'])}-{int(locus_information['Refined_Empty_right_bound'])+len_tsd-1}"
                        print(empty_site_cmd_out)
                    #print(empty_site_cmd)
                    #sys.exit("handle")
                
                locus_information["Age_TSD_Seq_left"] = left_cmd_out
                locus_information["Age_TSD_Seq_right"] = right_cmd_out
                locus_information["Age_TSD_Seq_empty"] = empty_site_cmd_out
                #print(left_cmd_out,right_cmd_out,empty_site_cmd_out)
                if left_cmd_out != right_cmd_out or left_cmd_out != empty_site_cmd_out:
                    print(left_cmd_out,right_cmd_out,empty_site_cmd_out)
                    locus_information['mismatch_TSD'] = True
                    #sys.exit("Left tsd either does not match right or does not match empty site")
                else:
                    locus_information['mismatch_TSD'] = False
                #print(locus_information)

            else:
                #no TSD
                locus_information['left_TSD_filled'] = ["N/A","N/A"]
                locus_information['right_TSD_filled'] = ["N/A","N/A"]
                locus_information["Age_TSD_Seq_left"] = "N/A"
                locus_information['Age_TSD_Seq_right'] = "N/A"
                locus_information['Age_TSD_Seq_empty'] = "N/A"
                locus_information['empty_coords'] = "N/A"
                locus_information['mismatch_TSD'] = False
            if age_out[i+2].split()[3] != "0":
                #handle
                locus_information['left_TSD_empty'] = [int(age_out[i+2].split()[5].split(",")[0])+locus_information['Empty_start']-locus_information['age_flank']-1,int(age_out[i+2].split()[5].split(",")[1])+locus_information['Empty_start']-locus_information['age_flank']-1]
                locus_information['right_TSD_empty'] = [int(age_out[i+2].split()[7].split(",")[0])+locus_information['Empty_start']-locus_information['age_flank']-1,int(age_out[i+2].split()[7].split(",")[1])+locus_information['Empty_start']-locus_information['age_flank']-1]
            else:
                #no TSD
                locus_information['left_TSD_empty'] = ["N/A","N/A"]
                locus_information['right_TSD_empty'] = ["N/A","N/A"]
            #sys.exit("bottom of function")
            #print(comp)
            #print(locus_information["Refined_Filled_start"],locus_information["Refined_Filled_end"])
            #sys.exit()
        else:
            continue
        
def process_blat_file(input_blat,locus_information,left,right): #ToDo, what if multiple with different orientation, what if does not reach the end of the element?
    with open(input_blat,'rt') as infile:
        file_contents = infile.readlines()

    if len(file_contents) < 6:
        locus_information["percent_identical"] = "N/A"
        locus_information["match_count"] = "N/A"
        locus_information['element_blat_orientation'] = "N/A"
        locus_information['blat_aln_end'] = "N/A"
        locus_information['blat_aln_end_ref_coord'] = "N/A"
        locus_information['multiple_blat_segments'] = False
        #sys.exit("No blat alignment found") #Todo figure out what happens here
        return
    else:
        #obtain general information
        first_output = file_contents[5].split()
        print(first_output)
        locus_information['multiple_blat_segments'] = False
        #if first_output[8] != "+":
        #    sys.exit("Reverse blat file orientation")
        if len(file_contents) == 6: #exactly one output. This is the easiest processing 
            
            #get the blocks. This is because if there are gaps in the alignment in which ct repeats are present in the reference element but absent in the query, I do not want the gapped sequence counted as length in the alignment.
            len_blocks = [int(x) for x in first_output[18].split(",")[:-1]]
            blocks_query = [int(x) for x in first_output[19].split(",")[:-1]]
            blocks_target = [int(x) for x in first_output[20].split(",")[:-1]]
            if first_output[8] == "+":
                locus_information['element_blat_orientation'] = "Forward"
            else:
                locus_information['element_blat_orientation'] = "Reverse"
            if locus_information['element_blat_orientation'] == "Forward":
                locus_information['blat_aln_end'] = int(first_output[12])-1
            else:
                locus_information['blat_aln_end'] = int(first_output[11])
            locus_information['blat_aln_end_ref_coord'] = locus_information['blat_aln_end'] + left
            #query_ranges = []
            target_ranges = []
            print(blocks_query)
            for i in range(len(blocks_query)):
                #query_ranges.append([blocks_query[i]+1,blocks_query[i]+len_blocks[i]])
                target_ranges.append([blocks_target[i]+1,blocks_target[i]+len_blocks[i]])
            #print(query_ranges)
            target_gaps = []
            print(target_ranges)
            for i in range(len(target_ranges)-1):
                target_gaps.append([target_ranges[i][1]+1,target_ranges[i+1][0]-1])
            #todo edge case for 0 gaps
            #Reduce the length of the ct by the number of bp overlapping with gaps.
            ct_region = [132,149]
            ct_length = ct_region[1]-ct_region[0]+1
            overlap = []
            for gap in target_gaps:
                if gap[0] >ct_region[1] or gap[1] < ct_region[0]: #does not overlap ct
                    continue 
                elif gap[0] <= ct_region[0] and gap[1] >= ct_region[1]: #overlaps the entirety of the ct
                    gap[0] = ct_region[0]
                    gap[1] = ct_region[1]
                elif gap[0] <= ct_region[0]: #starts at left side of ct and does not reach end
                    gap[0] = ct_region[0]
                elif gap[1] >= ct_region[1]: #starts in the ct repeat and may reaches the end
                    gap[1] = ct_region[1]
                elif gap[0] > ct_region[0] and gap[1] < ct_region[1]: #starts and ends in the middle. I doubt this would happen, but if it does, no changes are required.
                    pass
                else:
                    sys.exit("Unexpected blat parse")
                for i in range(gap[0],gap[1]+1):
                    overlap.append(i)
            overlap = set(overlap)
            if overlap == []:
                ct_overlap_final = 0
            else:
                ct_overlap_final = len(overlap)
            print(first_output)
            
            

        else: #Todo. Get this online.
            locus_information["percent_identical"] = "N/A"
            locus_information["match_count"] = "N/A"
            locus_information['element_blat_orientation'] = "N/A"
            locus_information['blat_aln_end'] = "N/A"
            locus_information['blat_aln_end_ref_coord'] = "N/A"
            locus_information['multiple_blat_segments'] = True
            return 
            #sys.exit("Multiple blat alignments. Requires seperate processing") #Todo impliment this
            

        
        length_TE = int(first_output[14])-ct_overlap_final
        length_insertion = int(first_output[10])
        print(length_insertion,length_TE)
        
        match_count = int(first_output[0])
        fraction_identical = match_count/length_TE
        print(fraction_identical)
        print(match_count)
        locus_information["percent_identical"] = fraction_identical
        locus_information["match_count"] = match_count
        #print(locus_information)
        
       
        
def homopolymer_extraction(locus_information):   
    #extract relevant locus.
    if locus_information["RM_orientation"] == "Forward":
        nucleotide = "a"
    elif locus_information["RM_orientation"] == "Reverse":
        nucleotide = "t"
    else:
        sys.exit("Neither forward nor reverse orientation")
        
    if locus_information["right_TSD_filled"][0] != "N/A":
        start = int(locus_information['left_TSD_filled'][1])+1
        samtools_extract_coords = f"{locus_information['chr']}:{start}-{int(locus_information['right_TSD_filled'][0])-1}"
    else:
        start = int(locus_information['Refined_Filled_start'])
        samtools_extract_coords = f"{locus_information['chr']}:{start}-{int(locus_information['Refined_Filled_end'])}"
        
    ##extract relevant locus.
    #if locus_information["element_blat_orientation"] == "Forward":
    #    start = locus_information['blat_aln_end_ref_coord']+1
    #    nucleotide = "a"
    #    if locus_information["right_TSD_filled"][0] != "N/A":
    #        samtools_extract_coords = f"{locus_information['chr']}:{start}-{int(locus_information['right_TSD_filled'][0])-1}"
    #    else:
    #        samtools_extract_coords = f"{locus_information['chr']}:{start}-{int(locus_information['Refined_Filled_end'])}"
    #elif locus_information["element_blat_orientation"] == "Reverse":
    #    nucleotide = "t"
    #    if locus_information["left_TSD_filled"][1] != "N/A":
    #        start = int(locus_information['left_TSD_filled'][1])+1
    #        samtools_extract_coords = f"{locus_information['chr']}:{start}-{locus_information['blat_aln_end_ref_coord']-1}"
    #    else:
    #        start = locus_information['Refined_Filled_start']
    #        samtools_extract_coords = f"{locus_information['chr']}:{start}-{locus_information['blat_aln_end_ref_coord']-1}"
    #
    #else:
    #    sys.exit("Neither forward nor reverse orientation")
        

    print(samtools_extract_coords)
    cmd = f"samtools faidx {locus_information['Ref']} {samtools_extract_coords}"
    poly_extract = extract_sequence_from_faidx_cmd(cmd).lower()
    print(poly_extract,start)
    
    poly_list = []
    current_bp = start
    print(start)
    for i in range(len(poly_extract)):
        if poly_list == []:
            if poly_extract[i] != nucleotide:
                pass
            else:
                poly_list.append([current_bp])
        elif poly_extract[i] != nucleotide:
            if len(poly_list[-1]) == 1:
                poly_list[-1].append(current_bp-1)
            else:
                pass
        elif poly_extract[i] == nucleotide:
            if len(poly_list[-1]) == 1:
                pass
            else:
                poly_list.append([current_bp])
        current_bp+=1
    if len(poly_list) == 0:
        print(poly_list)
        poly_list = [['N/A']]
    elif len(poly_list[-1]) == 1:
        poly_list[-1].append(current_bp-1)
        #print(poly_extract[i])
    print(poly_list)
    if poly_list == []:
        poly_list = "N/A"
    locus_information["poly_a_list"] = poly_list
    
def process_file(input_file):
    """
    Extract the query, reference, and between sequences from a smith-waterman alignment
    
    Args:
        input_file (str): Path to file from which the alignment is being extracted
    
    Returns:
        seq1: (list) containing the first 13 characters of the a sequence, start position the aligned portion of the A sequence, the aligned nucleotides of the A sequence, the end position of the aligned A sequence
        seq2: (list) containing the first 13 characters of the b sequence, start position the aligned portion of the B sequence, the aligned nucleotides of the B sequence, the end position of the aligned B sequence
    """
    with open(input_file) as infile:
                file_contents = infile.readlines()[32:-3:]
                i = 0
                seq1 = ["","","",""]
                seq2 = ["","","",""]
                while i < len(file_contents):
                    content = file_contents[i].split()
                    if i%4 == 0:
                        if seq1[0] == "":
                            seq1[0] = (content[0])
                            seq1[1] = (content[1])
                        seq1[2] = seq1[2] + content[2]
                        if i+4 == len(file_contents):
                            seq1[3] = (content[3])
                    elif i%4 == 2:
                        if seq2[0] == "":
                            seq2[0] = (content[0])
                            seq2[1] = (content[1])
                        seq2[2] = seq2[2] + content[2]
                        if i+2 == len(file_contents):
                            seq2[3] = (content[3])
                    i+=1
    return(seq1,seq2)

def investigate_TSD_validity(flank_dist,dist,seq1,seq2,length,search_internal,interior_dist,internal_TSD_upstream=False):
    """
    Detects if target site duplications are near within a predetermined distance of the locus of interest.
    Additionally, determines if flanks need to be extended (this occurs when the TSD is within the predetermined distance, but can be extended beyond the currently extracted sequence
    
    Args:
        flank_dist: (int) the number of base pairs that are extracted beyond the locus of interest 
        dist: (int) detected TSDs must be within dist bp of the start of the locus to be considered valid_TSD
        seq1: (list) containing the first 13 characters of the a sequence, start position the aligned portion of the A sequence, the aligned nucleotides of the A sequence, the end position of the aligned A sequence
        seq2: (list) containing the first 13 characters of the a sequence, start position the aligned portion of the B sequence, the aligned nucleotides of the B sequence, the end position of the aligned B sequence
        length: (int) length of the insertion. This will be used in future version in which non-reference hallmark detection is possible
        search_internal: (bool) Determines if TSDs can be detected within the locus of interest. Should be set to False for retrogene detection. Later, True will be added for non-reference hallmark detection
        interior_dist: (int) The number of base pairs within the locus of interest that will be used to search for TSDs. Can be used when search_internal is False
        
    Returns:
        bool: True if the flank has been extended and the alignment is to be re-performed.
        bool: True if one or both of the sequences are invalid (as determined by proximity to the locus of interest). If True, the above bool is also set to False     
    """
    alignment_fail = []

    if search_internal == True:
        #sys.exit("search_internal will be enabled in a future update. Please set search_internal == False")
         
        #to be added at a later date for non-reference retro-intersetions
        for seq in [seq1,seq2]:
            print(flank_dist)
            print(seq)
            if internal_TSD_upstream == False:
                if abs(flank_dist - int(seq[1])) > dist and abs(flank_dist - int(seq[3])) > dist:
                    if seq == seq1:
                        alignment_fail.append("seq1")
                    else:
                        alignment_fail.append("seq2")
            else:
                #print('here')
                if abs(int(seq[1])-interior_dist) > dist and abs(int(seq[3])-interior_dist) > dist:
                        if seq == seq1:
                            alignment_fail.append("seq1")
                        else:
                            alignment_fail.append("seq2")
                #sys.exit()
        if len(alignment_fail) >= 1:
            print(f"At least 1 alignment failed to be withing {dist}bp of the start/end of the inserted sequence")
            print(seq1,seq2)
            return False, False
        else:
            print('here',internal_TSD_upstream)
            if internal_TSD_upstream == False:
                if (int(seq1[1]) <= 5) or (int(seq2[1]) <= 5):
                    #print(length)
                    if flank_dist > length/2:
                        False,False
                        #sys.exit()
                    #print(seq1,seq2)
                    print("TSD may be longer than extracted flank. Doubling flank and trying again")
                    subprocess.run("rm water.txt",shell=True)
                    return True, False
            else:
                if (flank_dist+interior_dist - int(seq1[3]) < 5) or (flank_dist+interior_dist - int(seq2[3]) < 5):
                    #print(length)
                    if flank_dist > length/2:
                        False,False
                        #sys.exit()
                    #print(seq1,seq2)
                    print("TSD may be longer than extracted flank. Doubling flank and trying again")
                    subprocess.run("rm water.txt",shell=True)
                    return True, False
    #else:
    #    print(seq1,seq2)
    #    print(dist,flank_dist)
    else:
        #sys.exit()
        if (int(seq1[1]) <= 5) or (flank_dist+interior_dist - int(seq2[3]) < 5):
            print("TSD may be longer than extracted flank. Extending flank and trying again")
            subprocess.run("rm water.txt",shell=True)
            return True, False
    return False, True
    
        
def run_water(flank_dist, dist, coords, extraction_path,ref,orientation, search_internal, chrom_length,internal_TSD_upstream=False,file_prefix=""):
    """
    Function to extract two sequences from a fasta file of interest and then perform a smith-waterman alignment using the EMBOSS module
    
    Args:
        flank_dist: (int) the number of base pairs that are extracted beyond the locus of interest 
        dist: (int) detected TSDs must be within dist bp of the start of the locus to be considered valid_TSD
        coords: (list) containing the sequence name, start base, and end base of the locus of interest. 1 based 
        extraction_path: (str) path to the fasta containing the sequence being extracted
        ref: (str) path to reference. As of now not used, but likely will be for non-reference detection
        orientation: (str) Forward or Reverse, the orientation of the retroelement in reference to the chromosome or contig
        search_internal: (bool) Determines if TSDs can be detected within the locus of interest. Should be set to False for retrogene detection. Later, True will be added for retrogene hallmark detection
        chrom_length: the length of the chromosome or contig in which the locus of interest resides
        
    Returns:
        seq1: (list) containing the first 13 characters of the a sequence, start position the aligned portion of the A sequence, the aligned nucleotides of the A sequence, the end position of the aligned A sequence, returns "" if no TSD detected
        seq2: (list) containing the first 13 characters of the b sequence, start position the aligned portion of the B sequence, the aligned nucleotides of the B sequence, the end position of the aligned B sequence, returns "" if no TSD detected
    """
    if file_prefix != "" and file_prefix[-1] != "_" and file_prefix[-1] != "/":
        file_prefix = file_prefix + "_"
    interior_dist = 0
    re_run = True
    while re_run:
        if search_internal == True:
            interior_dist = 5
            if internal_TSD_upstream == True:
                extract_up = f"{coords[0]}:{max(coords[1]-interior_dist,1)}-{coords[1]+(flank_dist-1)}"
                extract_down = f"{coords[0]}:{coords[2]-(interior_dist-1)}-{coords[2]+flank_dist}"
            elif internal_TSD_upstream == False:
                extract_up = f"{coords[0]}:{max(coords[1]-flank_dist,1)}-{coords[1]+(interior_dist -1)}"
                extract_down = f"{coords[0]}:{coords[2]-(flank_dist-1)}-{coords[2]+interior_dist}"
            #sys.exit("search_internal will be enabled in a future update. Please set search_internal == False")
            #to be added at a later date for non-reference retro-intersetions
            
            #extract_up = f"{coords[0]}:{max(coords[1]-flank_dist,1)}-{coords[1]+(flank_dist -1)}"
            #extract_down = f"{coords[0]}:{coords[2]-(flank_dist-1)}-{coords[2]+flank_dist}"
        elif search_internal != True:
            dist = flank_dist
            extract_up = f"{coords[0]}:{max(coords[1]-dist,1)}-{coords[1]+(interior_dist-1)}"
            extract_down = f"{coords[0]}:{coords[2]-(interior_dist-1)}-{coords[2]+dist}"

        command = f"samtools faidx {extraction_path} {extract_up} > {file_prefix}upstream_and_5_start.fa ; samtools faidx {extraction_path} {extract_down} > {file_prefix}downstream_and_3_end.fa"
        print(command)
        subprocess.run(command, shell=True)    
        gap_open = 10
        gap_extend = 10
        custom_matrix_path = os.path.dirname(__file__)
        cmd = f"water {file_prefix}upstream_and_5_start.fa {file_prefix}downstream_and_3_end.fa -gapopen {gap_open} --gapextend {gap_extend} -datafile {custom_matrix_path}/Custom_Matrix -outfile {file_prefix}water.txt"
        #print(cmd)
        subprocess.run(cmd, shell=True)
        seq1,seq2 = process_file(f"{file_prefix}water.txt")   
        re_run, valid_TSD = investigate_TSD_validity(flank_dist,dist,seq1,seq2,coords[2]-coords[1]+1,search_internal,interior_dist,internal_TSD_upstream)
        
        #Check for edge cases
        if coords[1]-flank_dist <= 1:
            re_run = False
        elif coords[2]+dist >= chrom_length:
            re_run = False
        elif len(seq1[2]) == 1:
            re_run = False
        elif "N" in seq1[2] or "N" in seq2[2]:
            sys.exit("Masked based in aligned sequence.")
        if re_run == True:
            flank_dist = flank_dist+5
            continue
            
        if valid_TSD != True:
            print("No Valid TSD detected")
            return "",""
        else:
            #print(seq1,seq2)
            if seq1[2].lower() == seq2[2].lower():
                print("Seq1 and Seq2 are identical")
            return seq1,seq2
            
def Refine_Coords(seq1,seq2,extract):
    if seq1 != "":
        print(seq1,seq2,extract)
        #print(int(seq1[0].split("-")[0])+int(seq1[3])+1,int(seq2[0].split("-")[0])-2+int(seq1[1]))
        start_after_TSD = int(seq1[0].split("-")[0])+int(seq1[3])
        end_before_TSD = int(seq2[0].split("-")[0])-1+int(seq1[1])-1
        RM_coords = [extract[0],start_after_TSD,end_before_TSD]
        TSDs = [extract[0],int(seq1[0].split("-")[0])+int(seq1[1])-1,seq1[2],int(seq1[0].split("-")[0])+int(seq1[3])-1,int(seq2[0].split("-")[0])+int(seq2[1])-1,seq2[2],int(seq2[0].split("-")[0])+int(seq2[3])-1]
        #print(TSDs)
    else:
        RM_coords = extract
        TSDs = ["N/A","N/A","N/A","N/A","N/A","N/A","N/A"]    
    return(RM_coords,TSDs)
    
def Repeat_Positions(repeat_masker_output, RM_position, RM_Class):
    #load in the relevant portions of the out file.
    if repeat_masker_output.endswith(".bed"):
        command = f"cat {repeat_masker_output} | grep -w '^{RM_position[0]}'"
        print(command)
        proc = subprocess.Popen(command,shell=True,stdout=subprocess.PIPE)
        rm_out,err = proc.communicate()
        rm_out = rm_out.decode()
        rm_out = rm_out.split("\n")
        #print(rm_out)

        start_index = 1
        end_index = 2
        offset = 1
        ME_Type = 4
        ME_Subtype = 3
        O = 8

    else:
        command = f"cat {repeat_masker_output}"
        proc = subprocess.Popen(command,shell=True,stdout=subprocess.PIPE)
        rm_out,err = proc.communicate()
        rm_out = rm_out.decode()
        rm_out = rm_out.split("\n")
        #print(rm_out)

        start_index = 5
        end_index = 6
        offset = 0
        ME_Type = 10
        O = 8

    print(f"RM for Repeat_Positions is {repeat_masker_output}")
    #new method of determining % L1 content.
    lower_bound = RM_position[1]
    upper_bound = RM_position[2]

    L1s = []
    #with open(repeat_masker_output, 'rt') as file:
    for line in rm_out:
        if line.startswith("There were no repetitive sequences detected"):
            break
        if line.startswith("   SW") or line.startswith("score") or line == "\n" or line == "":
            continue
        line = line.split()
        #Define relevant parameters
        #print(line)
        RM_start_in_seq = int(line[start_index])+ offset
        RM_end_in_seq = int(line[end_index])
        Element_Class = line[ME_Type]
        Element_Orientation = line[O]
        Element_Subclass = line[ME_Subtype]

        #create list containing the location of the sequence relative to the LINE-1
        #ToDo Is the LINE-1 internal coordinates 0 based or 1 based within the bedfiles files.
        if repeat_masker_output.endswith(".bed") and Element_Orientation == "+":
            orientation = "Forward"
            Position_within_Repeat = [int(line[9].split(',')[0]),int(line[9].split(',')[1])]
        elif repeat_masker_output.endswith(".bed") and Element_Orientation == "C":
            Position_within_Repeat = [int(line[9].split(',')[2]),int(line[9].split(',')[1])]
            orientation = "Reverse"
        elif Element_Orientation == "+":
            orientation = "Forward"
            Position_within_Repeat = [int(line[11]),int(line[12])]
        elif Element_Orientation == "C":
            Position_within_Repeat = [int(line[13]),int(line[12])]
            orientation = "Reverse"
        line_start = 0
        line_end = 0

        if Element_Class == RM_Class:
            #print(f"Lower and Upper Bounds of the detected insertion are {lower_bound} and {upper_bound}")
            #print(f"Start and end of detected repeat are {RM_start_in_seq} {RM_end_in_seq}")
            #print(line)
            #if RM is outside of the insertion
            if RM_end_in_seq < int(lower_bound)+20:
                #print("L1 outide insertion")
                continue
            elif RM_start_in_seq > int(upper_bound)-20:
                #print("L1 outide insertion")
                break
            #if RM extends beyond both ends of the insertion    
            elif RM_start_in_seq < int(lower_bound) and RM_end_in_seq > int(upper_bound):
                #print(line)
                line_start = int(lower_bound)
                line_end = int(upper_bound)
                #print("L1 extends beyond insertion boundaries")

            #if the RM covers one of the boundaries of the insertion
            elif RM_start_in_seq <= int(lower_bound) and RM_end_in_seq <= int(upper_bound):
                #print(line)
                line_start = int(lower_bound)
                line_end = RM_end_in_seq
                #print("L1 overlaps one boundary")
            elif RM_start_in_seq >= int(lower_bound) and RM_start_in_seq <= int(upper_bound) and RM_end_in_seq >= int(upper_bound):
                #print(line)
                line_start = RM_start_in_seq
                line_end = int(upper_bound)
                #print("L1 overlaps one boundary")

            #if the RM is entirely within the insertion
            elif RM_start_in_seq > int(lower_bound) and RM_start_in_seq < int(upper_bound) and RM_end_in_seq < int(upper_bound):
                #print(line)
                line_start = RM_start_in_seq
                line_end = RM_end_in_seq
                #print("L1 is within insertion boundary")
            else:
                print("Erroneous Read Detected")
                print(line)
                sys.exit("Erroneous RM Detected")
            #print("Appending L1")
            #print(line)
            L1s.append([line_start, line_end, orientation, Position_within_Repeat,[Element_Subclass]])

    #print(initial_range)
    #print(L1s)


    #combine nearby elements.
    diff = 50 #TODO Refine the segment merging to take either the longest segment or the segment closest to the end of the LINE-1. It currently reports the segment closest to the 3' end of the insertion itself which is suboptimal.
    exit_loop = False
    if len(L1s) == 0:
        print("Length of L1s list is 0. Terminating")
        return "N/A", "N/A", "N/A", "N/A", "N/A"
    print(L1s)
    while not exit_loop:
        #print('hi',len(L1s))
        merged_L1 = []
        for i in range(len(L1s)):
            #print('hi')
            if i+1 == len(L1s):
                exit_loop = True
                break            
            #sys.exit()
            if L1s[i][2] == L1s[i+1][2]:
                #print(f"Orientations are the same",L1s[i][2],L1s[i+1][2])
                #sys.exit()
                if L1s[i][2] == "Forward" and abs(L1s[i+1][3][0] - L1s[i][3][1]) < diff:
                    #print(f"L1s are within range")
                    merged_L1 = [L1s[i][0],L1s[i+1][1],L1s[i][2],[L1s[i][3][0],L1s[i+1][3][1]]]
                    #print(merged_L1)
                elif L1s[i][2] == "Reverse" and abs(L1s[i][3][1] - L1s[i+1][3][0]) < diff:
                    #print(f"L1s are within range",L1s[i][3][1],L1s[i+1][3][0],diff)
                    merged_L1 = [L1s[i][0],L1s[i+1][1],L1s[i][2],[L1s[i][3][0],L1s[i+1][3][1]]]
                    #print(merged_L1)
                #else:
                    #print(f"L1s are close in contig distance, but do not overlap based on LINE-1 coordinates. The distance between the elements is larger than {diff}")

            else:
                #print(f"Orientations are different",L1s[i][2],L1s[i+1][2])
                #sys.exit()
                if max(int(L1s[i][3][0]),int(L1s[i][3][1])) > max(int(L1s[i+1][3][0]),int(L1s[i+1][3][1])):
                    #sys.exit()
                    if L1s[i][2] == "Forward":
                            #sys.exit()
                            continue
                    if abs(L1s[i][3][1] - L1s[i+1][3][1]) > diff:
                        #print("L1s are too far apart. Moving on")
                        continue
                    else:
                        #print("L1s are within range")
                        merged_L1 = [L1s[i][0],L1s[i+1][1],L1s[i][2],[L1s[i][3][0],L1s[i+1][3][0]]]
                        #print(merged_L1)
                        #sys.exit()
                    #sys.exit()
                else:
                    #diff2 = abs(max(L1s[i][0],L1s[i][1]) -  min(L1s[i+1][0],L1s[i+1][1]))
                    if L1s[i][2] == "Forward":
                        #sys.exit()
                        continue
                    if abs(L1s[i][3][0] - L1s[i+1][3][0]) > diff:
                        #print("L1s are too far apart. Moving on")
                        continue
                    else:
                        #print("L1s are within range")
                        merged_L1 = [L1s[i][0],L1s[i+1][1],L1s[i][2],[L1s[i][3][1],L1s[i+1][3][1]]]
                        #print(merged_L1)
                        #sys.exit()
            if merged_L1 != []:
                merged_L1.append(L1s[i][4])
                merged_L1[4].append(L1s[i+1][4][0])
                L1s[i] = merged_L1
                L1s.pop(i+1)
                #sys.exit()
                break
                            #print(merged_L1)
    #print(L1s)
    #return L1s
    #sys.exit()

    #sys.exit()
    #determine L1 content between TSDs
    #print(f"L1s is {L1s}")
    if len(L1s) > 1:
        L1_pos = "["
        L1_content = 0
        orientations = []
        for i in range(len(L1s)):
            item = L1s[i]
            #if i > 0:
            #    if item[0] <= L1s[i-1][1]:
            #        item[0] = L1s[i-1][1] + 1
            orientations.append(item[2])
            L1_pos = L1_pos + f"[{item[3][0]},{item[3][1]}]"
            #print(f"{L1_content} + ({item[1]}-{item[0]})+1")
            L1_content = L1_content + (item[1]-item[0])+1
        L1_pos = L1_pos + "]"
        print(L1_pos)
        transduction = ["N/A","N/A"]
        if len(set(orientations)) == 1:
            orientation = orientations[0]
        else: #TODO report the segments even if multiple orientations are detected
            orientation = "Multiple_L1_Orientations_detected"
        if orientation == "Multiple_L1_Orientations_detected":
            transduction = ["N/A","N/A"]
        elif orientation == "Forward":
            if int(L1s[-1][1]) - int(upper_bound) == 0:
                transduction = ["N/A","N/A"]
            else:
                transduction = [L1s[-1][1]+1,upper_bound]
            #print(transduction)
        else:
            if int(lower_bound) - int(L1s[0][0]) == 0:
                transduction = ["N/A","N/A"]
            else:
                transduction = [lower_bound,L1s[0][0]-1]
            #print(transduction)

    else:
        if L1s == []:
            L1_content = "No_L1_detected"
            transduction = ["N/A","N/A"]
            orientation = "No_L1_detected"
            L1_pos = "No_L1_detected"
            #print('goop')
        else:
            orientation = L1s[0][2]
            L1_content = (L1s[0][1]-L1s[0][0])+1
            L1_pos = f"[{L1s[0][3][0]},{L1s[0][3][1]}]"
            print(L1_content)
            if L1_content == 1:
                L1_content = 0
                transduction = ["N/A","N/A"]
            else:
                #putative 3' transduction detection
                if orientation == "Forward":
                    if int(L1s[0][1]) - int(upper_bound) == 0:
                        transduction = ["N/A","N/A"]
                    else:
                        transduction = [L1s[0][1]+1,upper_bound]
                    #print(transduction)
                else:
                    if int(lower_bound) - int(L1s[0][0]) == 0:
                        transduction = ["N/A","N/A"]
                    else:
                        transduction = [lower_bound,L1s[0][0]-1]
                    #print(transduction)
    
    print(f"Unfixed L1s is {L1s}")
    final_L1s = []
    for L1 in L1s:
        revised_L1 = L1[0:3]
        revised_L1.append(L1[3][0])
        revised_L1.append(L1[3][1])
        print("\\".join(L1[4]))
        revised_L1.append("\\".join(L1[4]))
        final_L1s.append(revised_L1)
    if len(final_L1s) == 0:
        final_L1s.append("N/A")
        
    print(final_L1s)
    #if len(L1s) > 1:
    #sys.exit()
    return L1_content, transduction, orientation, L1_pos, final_L1s
def extract_for_Poly_A(transduction,extraction_path,extracted_seq,extra_seq,orientation,internal,interior_dist,DATA_DICT_TEMP,exterior_cutoff=0):
    """
    Creates the samtools faidx command which will extract sequence for poly(A) detection.
        
    Args:
        transduction: Currently not used, will be needed in a future version for detecting non-reference retroelements
        extraction_path: (str) path to the fasta containing the sequence being extracted
        extracted_seq: (list) containing the sequence name, start base, and end base of the locus of interest. 1 based 
        extra_seq: (int) the number of base pairs that are extracted beyond the locus of interest 
        orientation: (str) Forward or Reverse, the orientation of the retroelement in reference to the chromosome or contig
        internal: (bool) Determines if Poly(A)s can be detected within the locus of interest. Should be set to False for retrogene detection. Later, True will be added for retrogene hallmark detection
        interior_dist: (int) The number of base pairs within the locus of interest that will be used to search for Poly(A)s. Can be used when search_internal is False
        exterior_cutoff: (int) The minimum/maximum (depending on orientation) distance that the poly(A) can be expanded into. This prevents expansion into TSDs. Default = 0
        
    Returns:
        faidx_cmd: (str) the samtools faidx command that will be used to extract sequence. If the TSD is directly adjacent to the locus of interest, no poly(A) will be found and faidx_cmd is equal to "No Poly A"
        start_in_contig: (int) The distance into the contig/chromosome of the first base extracted. 1 based.
    """
    
    if internal == True: #Todo consider if I should only do the interior_TE_dist if the 3' end of a TE is present.
        #print(DATA_DICT_TEMP)
        #sys.exit()
        if DATA_DICT_TEMP["orientation"] != "Forward" and DATA_DICT_TEMP["orientation"] != "Reverse":
            faidx_cmd = "No Poly A"
            start_in_contig = "No Poly A"
            return faidx_cmd, start_in_contig
        elif DATA_DICT_TEMP["orientation"] == "Forward":
            if transduction[0] != "N/A":
                start = max(DATA_DICT_TEMP["transduction"][0]-DATA_DICT_TEMP["Interior_TE_dist"],DATA_DICT_TEMP["L1_list"][-1][0]+1)
                end = DATA_DICT_TEMP["transduction"][1]
            else:
                start = max(DATA_DICT_TEMP["L1_list"][-1][1]-(DATA_DICT_TEMP["Interior_TE_dist"]-1),DATA_DICT_TEMP["L1_list"][-1][0]+1)
                if DATA_DICT_TEMP["TSDs"][4] != "N/A":
                    end = DATA_DICT_TEMP["TSDs"][4]-1
                else:
                    end = DATA_DICT_TEMP["Filled_end"]
        elif DATA_DICT_TEMP["orientation"] == "Reverse":
            if transduction[0] != "N/A":
                end = min(DATA_DICT_TEMP["transduction"][1]+DATA_DICT_TEMP["Interior_TE_dist"],DATA_DICT_TEMP["L1_list"][0][1]-1)
                start = DATA_DICT_TEMP["transduction"][0]
            else:
                end = min(DATA_DICT_TEMP["L1_list"][0][0]+(DATA_DICT_TEMP["Interior_TE_dist"]-1),DATA_DICT_TEMP["L1_list"][0][1]-1)
                if DATA_DICT_TEMP["TSDs"][3] != "N/A":
                    start = DATA_DICT_TEMP["TSDs"][3]+1
                else:
                    start = DATA_DICT_TEMP["Filled_start"]+1
            print("start",start,"end",end)
        else:
            sys.exit()
        faidx_cmd = f"samtools faidx {extraction_path} {extracted_seq[0]}:{start}-{end}"
        start_in_contig = start

        #sys.exit("internal == True will be enabled in a future update. Please set internal == False")        
        #if transduction[0] != "N/A":
        #    faidx_cmd = f"samtools faidx {extraction_path} {extracted_seq[0]}:{transduction[0]}-{transduction[1]}"
        #    start_in_contig = transduction[0]
        #else:
        #    faidx_cmd = "No Poly A"
        #    start_in_contig = "No Poly A"
            #don't forget to filter to not go earlier than the start of the element.
            #if orientation == "Forward":
            #    faidx_cmd = f"samtools faidx {extraction_path} {extracted_seq[0]}:{max(extracted_seq[2]-extra_seq,extracted_seq[1])}-{extracted_seq[2]}"
            #    start_in_contig = max(extracted_seq[2]-extra_seq,extracted_seq[1])
            #elif orientation == "Reverse":
            #    faidx_cmd = f"samtools faidx {extraction_path} {extracted_seq[0]}:{extracted_seq[1]}-{min(extracted_seq[1]+extra_seq,extracted_seq[2])}"
            #    start_in_contig = extracted_seq[1]#to be added at a later date for non-reference retro-intersetions
    elif internal == False:
        if exterior_cutoff == 0 and orientation == "Forward":
            exterior_cutoff = 100000000000
        if orientation == "Forward":
            faidx_cmd = f"samtools faidx {extraction_path} {extracted_seq[0]}:{extracted_seq[2]-(interior_dist-1)}-{min(extracted_seq[2]+extra_seq,exterior_cutoff)}"
            extract_length = min(extracted_seq[2]+extra_seq,exterior_cutoff)-(extracted_seq[2]-(interior_dist-1))
            start_in_contig = extracted_seq[2]-(interior_dist-1)
        elif orientation == "Reverse":
            faidx_cmd = f"samtools faidx {extraction_path} {extracted_seq[0]}:{max(extracted_seq[1]-extra_seq,exterior_cutoff)}-{extracted_seq[1]+(interior_dist-1)}"
            extract_length = (extracted_seq[1]+(interior_dist-1)) - max(extracted_seq[1]-extra_seq,exterior_cutoff)
            start_in_contig = max(extracted_seq[1]-extra_seq,exterior_cutoff)
        
        #occurs if the TSD starts at the first base of sequence. Thus no Poly(A) can be detected.
        if extract_length == -1:
            faidx_cmd = "No Poly A"
        elif extract_length < -1:
            sys.exit("Invalid length for poly(A) discovery")
    else:
        sys.exit()
    #print(start_in_contig)
    return faidx_cmd, start_in_contig
    
def Detect_Poly_As(transduction,extraction_path,extracted_seq,extra_seq,max_dist,orientation,internal,chrom_length,DATA_DICT_TEMP,exterior_cutoff=0):
    """
    Detects Poly(A) tails located within extra_seq base pairs of a locus of interest. 
        
    Args:
        transduction: Currently not used, will be needed in a future version for detecting non-reference retroelements
        extraction_path: (str) path to the fasta containing the sequence being extracted
        extracted_seq: (list) containing the sequence name, start base, and end base of the locus of interest. 1 based 
        extra_seq: (int) the number of base pairs that are extracted beyond the locus of interest 
        max_dist: (int) detected Poly(A)s must be within dist bp of the start of the locus to be considered a valid poly(A)
        orientation: (str) Forward or Reverse, the orientation of the retroelement in reference to the chromosome or contig
        internal: (bool) Determines if Poly(A)s can be detected within the locus of interest. Should be set to False for retrogene detection. Later, True will be added for retrogene hallmark detection
        chrom_length: the length of the chromosome or contig in which the locus of interest resides
        exterior_cutoff: (int) The minimum/maximum (depending on orientation) distance that the poly(A) can be expanded into. This prevents expansion into TSDs. Default = 0
        
    Returns:
        final_polys: (list) a list of lists where each sublist contains in order the start coordinate of each detected poly A, its sequence, the end coordinate of each poly A. All coordinates are 1 based. If no poly A is detected, an empty list is returned.
        terminated_early: (bool) Returns True if the final poly(A) would be extended if not for a limiting factor such as extending into the beginning/end of a chromosome or into a previously detected TSD. Primarily used for debugging. 
    """
    re_run = True
    if orientation != "Forward" and orientation != "Reverse":
        print("Multiple different orientations. Cannot resolve poly(A) tail")
        poly_a = []
        
    else:   
        while re_run == True:
            interior_dist = 0
            re_run = False
            

            faidx_cmd, start_in_contig = extract_for_Poly_A(transduction,extraction_path,extracted_seq,extra_seq,orientation,internal,interior_dist,DATA_DICT_TEMP,exterior_cutoff)
            print(faidx_cmd)
            trailing_sequence = ""
            
            if faidx_cmd == "No Poly A":
                poly_a = []
                break
                
            proc = subprocess.Popen(faidx_cmd,shell=True,stdout=subprocess.PIPE)
            faidx,err = proc.communicate()
            faidx = faidx.decode()
            
            if orientation == "Reverse":
                A_or_T = "T"
            elif orientation == "Forward":
                A_or_T = "A"
            poly_a = []
            poly_a_in_progress = ""
            start = 0
            end = None
            min_poly = 5
            for line in faidx.split("\n")[1::]:
                trailing_sequence = trailing_sequence + (line.rstrip())
            trailing_sequence = trailing_sequence.upper()
            for i in range(len(trailing_sequence)):
                if trailing_sequence[i] == A_or_T and poly_a_in_progress == "":
                    poly_a_in_progress = poly_a_in_progress + (trailing_sequence[i])
                    start = i + start_in_contig
                elif trailing_sequence[i] == A_or_T and not i == len(trailing_sequence) -1:
                    poly_a_in_progress  = poly_a_in_progress + (trailing_sequence[i])
                elif trailing_sequence[i] != A_or_T and poly_a_in_progress != "" or i == len(trailing_sequence)-1:
                    upcoming_poly = 0
                    end_of_extract = False
                    for j in range(i+1,i+6):
                        if j >= len(trailing_sequence):
                            end_of_extract = True
                            break
                        if trailing_sequence[j] == A_or_T:
                            upcoming_poly +=1
                    if upcoming_poly >= 4 or end_of_extract == True and i+1 != len(trailing_sequence):
                        poly_a_in_progress = poly_a_in_progress + trailing_sequence[i]
                    else:
                        #trim the ends of the poly if they have mutations in them.
                        if i+1 == len(trailing_sequence) and trailing_sequence[i] == A_or_T:
                            poly_a_in_progress = poly_a_in_progress + trailing_sequence[i]
                        
                        if not len(poly_a_in_progress) < min_poly:
                            end_poly = 3
                            trim_length_start = 0
                            while poly_a_in_progress[0:end_poly] != end_poly*A_or_T:
                                if len(poly_a_in_progress) < min_poly:
                                    poly_a_in_progress = ""
                                    break
                                else:
                                    poly_a_in_progress = poly_a_in_progress[1::]
                                    start+=1
                            while poly_a_in_progress[-1:-1*end_poly-1:-1] != end_poly*A_or_T:
                                if len(poly_a_in_progress) < min_poly:
                                    poly_a_in_progress = ""
                                    break
                                else:
                                    poly_a_in_progress = poly_a_in_progress[0:-1]
                            if len(poly_a_in_progress) >= min_poly:
                                if orientation == "Forward":
                                    #start = start - extra_seq
                                    end = start + len(poly_a_in_progress)-1
                                    if (extracted_seq[2] + extra_seq) - (end) < 5 and (exterior_cutoff == 0 or exterior_cutoff == "0"):
                                        if extracted_seq[2] + extra_seq >= chrom_length or internal == True:
                                            re_run = False
                                            print("Re-run not happening. At end of chromosome")
                                        else:
                                            re_run = True
                                            extra_seq +=5
                                            print(poly_a_in_progress)
                                            print("Re-run happening, Forward")
                                elif orientation == "Reverse":
                                    end = start + len(poly_a_in_progress)-1
                                    if abs((extracted_seq[1] - extra_seq) - (start)) <= 5 and (exterior_cutoff == 0 or exterior_cutoff == "0"):
                                        if extracted_seq[1] - extra_seq > 1 or internal != True:
                                            re_run = True
                                            extra_seq +=5
                                            print(poly_a_in_progress)
                                            print("Re-run happening, Reverse")
                                        else:
                                            print("Re-run not happening. Too close to start of chrom")
                                poly_a.append([start,poly_a_in_progress,end])                               
                        poly_a_in_progress = ""
    print(poly_a)
    final_polys = []
    
    if internal == False:
        terminated_early = False
        if len(poly_a) > 0:
            for poly in poly_a:
                if orientation == "Reverse":
                    if poly[2] < extracted_seq[1] - max_dist:
                        print("Out of range")
                    else:
                        final_polys.append(poly)
                        if abs(extracted_seq[1] - extra_seq) - poly[0] <= 5:
                            terminated_early = True
                elif orientation == "Forward":
                    if poly[0] > extracted_seq[2] + max_dist:
                        print("Out of range")
                    else:
                        final_polys.append(poly)
                        if (extracted_seq[2] + extra_seq) - poly[2] < 5:
                            terminated_early = True
                else:
                    continue
                    sys.exit()
            print(final_polys)
        return final_polys, terminated_early
    else:
        return poly_a
    
def transduction_detection(locus_info):
    
    if locus_info["strict_processing"] == True:
        if len(locus_info["poly_a"]) < 2 or locus_info["poly_a"] == [["N/A","N/A","N/A"]]:
            return ["N/A","N/A"]
        else:
            Transduction = [locus_info["poly_a"][0][2]+1,locus_info["poly_a"][-1][0]-1]
            return Transduction
    else:
        #print(locus_info["poly_a"],locus_info["TSDs"])
        if locus_info["orientation"] == "Forward":
            if "N/A" in locus_info["transduction"]:
                Transduction = locus_info["transduction"]
            elif "N/A" in locus_info["poly_a"][0] and "N/A" in locus_info["TSDs"]:
                Transduction = locus_info["transduction"]
            elif "N/A" in locus_info["poly_a"][0]:
                Transduction = [locus_info["transduction"][0],locus_info["TSDs"][4]-1]
            elif "N/A" in locus_info["TSDs"]:
                if locus_info["Filled_end"] - locus_info["poly_a"][-1][2] <= 6:
                    if locus_info["poly_a"][-1][0] == locus_info["transduction"][0]:
                        Transduction = ["N/A","N/A"]
                    else:
                        Transduction = [locus_info["transduction"][0],locus_info["poly_a"][-1][0]-1]
                else:
                    Transduction = locus_info["transduction"]
            #print(locus_info["TSDs"][4],locus_info["poly_a"][-1][2])
            elif locus_info["TSDs"][4] - locus_info["poly_a"][-1][2] <= 6: #cannot be more than 5bp between the TSD and poly(A)
                #print(locus_info["transduction"])
                if locus_info["poly_a"][-1][0] == locus_info["transduction"][0]:
                    Transduction = ["N/A","N/A"]
                else:
                    Transduction = [locus_info["transduction"][0],locus_info["poly_a"][-1][0]-1]
            elif locus_info["TSDs"][4] - locus_info["poly_a"][-1][2] > 6:
                Transduction = [locus_info["transduction"][0],locus_info["TSDs"][4]-1]
            else:
                print("Impossible Transduction")
                sys.exit()
        elif locus_info["orientation"] == "Reverse":
            print(locus_info["poly_a"],locus_info["TSDs"])
            if "N/A" in locus_info["transduction"]:
                Transduction = locus_info["transduction"]
            elif "N/A" in locus_info["poly_a"][0] and "N/A" in locus_info["TSDs"]:
                Transduction = locus_info["transduction"]
            elif "N/A" in locus_info["poly_a"][0]:
                Transduction = [locus_info["TSDs"][3]+1,locus_info["transduction"][1]]           
            elif "N/A" in locus_info["TSDs"]:
                if locus_info["poly_a"][0][0] - locus_info["Filled_start"] <= 6:
                    if locus_info["poly_a"][0][2] == locus_info["transduction"][1]:
                        Transduction = ["N/A","N/A"]
                    else:
                        Transduction = [locus_info["poly_a"][0][2]+1,locus_info["transduction"][1]]
                else:
                    Transduction = locus_info["transduction"]
            elif locus_info["poly_a"][0][0] - locus_info["TSDs"][3] <= 6: #cannot be more than 5bp between the TSD and poly(A)
                #print(locus_info["transduction"])
                if locus_info["poly_a"][0][2] == locus_info["transduction"][1]:
                    Transduction = ["N/A","N/A"]
                else:
                    Transduction = [locus_info["poly_a"][0][2]+1,locus_info["transduction"][1]]
            elif locus_info["poly_a"][0][0] - locus_info["TSDs"][3] > 6:
                Transduction = [locus_info["TSDs"][3]+1,locus_info["transduction"][1]]
            else:
                print("Impossible Transduction")
                sys.exit()
        else:
            return["N/A","N/A"] #need to know which orientation the LINE-1 is to detect a 3' transduction
        return Transduction
def identify_endo_site(orientation,TSDs,ref):
    print(orientation,TSDs)
    if orientation == "Forward" or orientation == "+":
        #extract coordinate for the start.
        print('start')
        start = TSDs[0].replace(":","-").split("-")[1]
        #print(start)
        
        locus_extract = f"{TSDs[0].replace(':','-').split('-')[0]}:{int(start)-2}-{int(start)+4}"
        print(locus_extract)
        
    elif orientation == "Reverse" or orientation == "-":
        end = start = TSDs[2].replace(":","-").split("-")[2]
        locus_extract = f"{TSDs[0].replace(':','-').split('-')[0]}:{int(end)-4}-{int(end)+2}"
        print(locus_extract)
        #sys.exit()
    
    faidx_cmd = f"samtools faidx {ref} {locus_extract}"
    proc = subprocess.Popen(faidx_cmd,shell=True,stdout=subprocess.PIPE)
    faidx,err = proc.communicate()
    faidx = faidx.decode()
    print(faidx)
    
    if orientation == "Forward" or orientation == "+":
        final_seq = revcomp(faidx.split("\n")[1])
    else:
        final_seq = faidx.split("\n")[1]
        #print(final_seq)
        #sys.exit()
    print(final_seq)
    #sys.exit()
    return final_seq
