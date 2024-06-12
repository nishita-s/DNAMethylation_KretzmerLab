# After Running RepeatMasker on my input fasta file, I got an output file with all the read information and the Transposable Elements
# I converted the output file to a bed file using RMout_to_bed.pl (https://github.com/4ureliek/Parsing-RepeatMasker-Outputs/blob/master/RMout_to_bed.pl)
# Using the following script, I split the obtained bedfile into smaller bedfiles with reads of the same TE types


import pandas as pd
import os as os
import re

os.chdir(path="/path/to/directory")

# Reading the output bedfile from Repeat Masker
df_NFZout_bed = pd.read_csv("NFZv2.0.fasta.out.bed", sep = "\t", header = None)

# df_NFZout_bed.head()

# Take the third column and get the item between the 10th and 11th semicolons
# Split it using separator "/" into TE family and TE subfamily

df_TE_Families = df_NFZout_bed[3].str.split(";", expand=True)
df_TE_Families = df_TE_Families[10]
df_TE_Families = df_TE_Families.str.split("/", expand=True)

#df_TE_Families.head()

# Array of TE_Families (not subfamilies) 
TE_Families = df_TE_Families[0].unique()
# array(['LINE', 'DNA', 'LTR', 'unclear', 'SINE', 'Unknown'], dtype=object)

# Attaching columns with TE family information separately
df_merged = pd.concat([df_NFZout_bed, df_TE_Families], axis=1, ignore_index=True) 

df_DNA = df_merged[df_merged[6]== "DNA"]
df_LINE = df_merged[df_merged[6]== "LINE"]
df_SINE = df_merged[df_merged[6] == "SINE"]
df_LTR = df_merged[df_merged[6]== "LTR"]
df_Unknown = df_merged[df_merged[6]== "Unknown"]
df_unclear = df_merged[df_merged[6]== "unclear"]

# Now we have split the main dataframe into separate bedfiles based on the TE Family
# We drop the last two columns about TE information before saving the dataframe into a bedfile

df_DNA_bed = df_DNA.drop(df_DNA.columns[-2:], axis=1)
df_LINE_bed = df_LINE.drop(df_LINE.columns[-2:], axis=1)
df_SINE_bed = df_SINE.drop(df_SINE.columns[-2:], axis=1)
df_LTR_bed = df_LTR.drop(df_LTR.columns[-2:], axis=1)
df_Unknown_bed = df_Unknown.drop(df_Unknown.columns[-2:], axis=1)
df_unclear_bed = df_unclear.drop(df_unclear.columns[-2:], axis=1)

# Combining "Unknown" and "Unclear" TEs because what even is the difference (ToT)
df_uncertain_bed = pd.concat([df_Unknown_bed, df_unclear_bed], ignore_index=True)
# df_uncertain_bed.head()

# Saving the dataframes into bedfiles

df_DNA_bed.to_csv("./path/DNA_TE_NFZ.bed", sep='\t', index=False, header=False)
df_LINE_bed.to_csv("./path/LINE_TE_NFZ.bed", sep='\t', index=False, header=False)
df_SINE_bed.to_csv("./path/SINE_TE_NFZ.bed", sep='\t', index=False, header=False)
df_LTR_bed.to_csv("./path/LTR_TE_NFZ.bed", sep='\t', index=False, header=False)
df_uncertain_bed.to_csv("./path/uncertain_TE_NFZ.bed", sep='\t', index=False, header=False)
