Github containing scripts associated with the publication "A phylogenetic estimate of canine retrotransposition rates based on genome assembly comparisons"

This github houses the most essential scripts associated with producing the data.

Note that Detecting_Hallmarks_Functions.py is a heavily modified revamp of the script by the same name used by Nguyen, Blacksmith, and Kidd PMID: 38946312

2025_08_06
This README is designed to be the most straightforward way to view how the data for the LINE-1 and SINE manuscript was created.
Any changes to scripts will be in the regular README.txt file, this will simply summarize everystep of generating the data.


Scripts:
convert-RM-to-bed.py converts a RepeatMasker outfile into a bedlike file compatable with other scripts.
Prepare_RM_Files.py Prepares the bedlike RepeatMasker files for use in the variant detection and characterization steps.
Performing_Alignments.py is used to align one reference genome to another.
Process_Paftools_Output.py processes the output from the alignment, isolating SNPs and SNVs as well as identifying the autosomal genome length of the alignment.
Create_Detect_Hallmarks_CMDs.py creates the commands that Detect_Hallmarks_Revamp.py uses to characterize the variants. Detect_Hallmarks_Revamp.py relies upon Detecting_Hallmarks_Functions.py.
Post_Processing_after_Hallmarks.py handles post processing as well as some miscellanious analyses such as TSD and poly(A) distributions, estimating the rate of retrotransposition, and subfamily fraction calculations. This script utilizes Blacksmith_Post_Processing_Utils.py.
python Aggregate_Calls_Revised_List.py  aggregates loci together while tracking which genomes possess each variant.
Loci_By_Chrom_Revised_List.ipynb processes the output from Aggregate_Calls_Revised_List.py and generates UpSet plots.
Transduction_Detection.py identifies putative LINE-1 3' transductions.
Transduction_Blat_Driver.sh is a command line utility to run the blat alignment of the transductions against the G_WOLF genome
Process_Blat_PSLX_Output.ipynb processes the blat pslx outfiles, calculates transduction fractions, and identifies parentless transductions.
Compare_Alignments.ipynb generates counts of LINE-1s and SINECs per chromosome per sample and plots the result.
Poly_A_Check.ipynb is used to calculate the fractions of variants with TSDs and Poly(A)s using less restrictive parameters.
Manuscript_Statistics.ipynb handles some general processing and calculates the mean and standard deviation for estimates such as the rate of retrotransposition.
Identify_Target_Site_Deletions.py as the name implies can be used to identify which variants possess target site deletions, and reports some summary statistics.
python Identify_heterozygous_dimorphic_SINEs.py identifies heterozygous variants in GSD1.
Process_Locus_for_Diagram.py can be used to easy comparisons of filled sites and empty sites from a list of variants. The script relies upon Produce_Miropeats_Image.py and requires manual curation to ensure that variants are visualized properly.
