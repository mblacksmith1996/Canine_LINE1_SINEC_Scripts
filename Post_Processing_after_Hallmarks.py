import Python_Functions_Post_Processing.Blacksmith_Post_Processing_Utils  as PFPP
import sys
import pkg_resources
import os

import numpy as np
import pandas as pd
import matplotlib as mpl
import logomaker
import scipy.stats
import argparse

#parse input args
parser = argparse.ArgumentParser(prog="Post_Processing_after_Hallmarks.py",description="Processes the fully annotated LINE-1 and SINEC files")
parser.add_argument("--input_directory",dest="in_dir",required=True,help="The path to the directory containing files to be parsed")
args = parser.parse_args()




#extract version information for modules (if available)
loaded_modules = PFPP.automatically_extract_modules(globals().copy(),sys.modules.keys())
for module in loaded_modules:
    if module == "builtins":
        continue
    try:   
        print(f"Version of {module} is {pkg_resources.get_distribution(module).version}")
    except:
        print(f"No version information for {module}")

#manually create dictionary of all alternate sample name. Not all samples are in all comparisons, but this will handle all cases for this dataset.
alt_name_dict = {"Mischka":"GSD1","Nala":"GSD2","Tasha":"BOX","Zoey":"GDN","Sandy":"DNG","Clu":"A_WOLF","mCanlor":"G_WOLF","Yella2":"LAB1","Ros":"LAB2","BD":"BMD1","OD":"BMD2"}


#determine the canines in the sample and which one is the refernce.
samples = []
ref = []

list_of_dir_contents = os.listdir(args.in_dir)
for file in list_of_dir_contents:
    if not file.endswith("_processed.bed"):
        continue
    file = file.split("_")
    if file[0] not in samples:
        samples.append(file[0])
    if file[1] not in ref:
        ref.append(file[1])
if len(ref) > 1:
    sys.exit("More than one reference file detected")
    
#sort the names according to dictionary order. This keeps the output in the order that is best for visualizing. 
final_sample_names = []
canine_alt_names = []
for key in alt_name_dict.keys():
    if key in samples:
        final_sample_names.append(key)
        canine_alt_names.append(alt_name_dict[key])
final_sample_names.append(ref[0])
canine_alt_names.append(alt_name_dict[ref[0]])


#create miscellanious information dictionary 
misc_info_dict = {"Directory_of_Data": args.in_dir,\
"canines":final_sample_names,\
"ref_name":ref[0],\
"canines_alt_names":canine_alt_names}#,\
misc_info_dict["Subdir_for_plots"] = f"{misc_info_dict['Directory_of_Data']}/Fig_Save"

#print(misc_info_dict)
#sys.exit()
if not os.path.exists(misc_info_dict["Subdir_for_plots"]):
    os.mkdir(misc_info_dict["Subdir_for_plots"])

canine_dict_LINE,canine_dict_SINE = "",""
canine_dict_LINE, canine_dict_SINE = PFPP.create_dictionary_of_LINE1_and_SINEC_dataframes(misc_info_dict,pd,sys)
#sys.exit()

Categories_dict_LINE, Auto_only_dict_LINE = PFPP.make_histograms(pd,mpl,logomaker,"LINE",canine_dict_LINE, misc_info_dict)
Categories_dict_SINE, Auto_only_dict_SINE = PFPP.make_histograms(pd,mpl,logomaker,"SINE",canine_dict_SINE, misc_info_dict)

PFPP.create_hallmark_stack_barplot(np,mpl,canine_dict_LINE,Categories_dict_LINE,"LINE",misc_info_dict)
PFPP.create_hallmark_stack_barplot(np,mpl,canine_dict_SINE,Categories_dict_SINE,"SINE",misc_info_dict)

PFPP.produce_rate_estimates(pd,np,mpl,misc_info_dict,canine_dict_LINE,Auto_only_dict_LINE,"LINE")
PFPP.produce_rate_estimates(pd,np,mpl,misc_info_dict,canine_dict_SINE,Auto_only_dict_SINE,"SINE")

PFPP.scatterplot_tsd_vs_polya(pd,np,scipy.stats,mpl,misc_info_dict,canine_dict_LINE,"LINE")
PFPP.scatterplot_tsd_vs_polya(pd,np,scipy.stats,mpl,misc_info_dict,canine_dict_SINE,"SINE")

PFPP.subfamily_analysis(pd,mpl,canine_dict_LINE,misc_info_dict,"LINE")
PFPP.subfamily_analysis(pd,mpl,canine_dict_SINE,misc_info_dict,"SINE")