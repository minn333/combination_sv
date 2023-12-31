import pandas as pd
import os

# combi = "intersect"  # choose this if intersect
combi = "union" # choose this if union

# model = "mdg"  # chooses this if combination of 3 callers
model = "mdgls"  # choose this if combination of 5 callers

infile = combi + ".HG002." + model + ".intv"
outfile = combi + ".HG002." + model + ".venn."


output_dir = "...\\combination_sv\\step4_performance\\models\\for_venn\\"
if not os.path.exists(output_dir):
    os.mkdir(output_dir)

input_dir = "...\\combination_sv\\step4_performance\\models\\for_intv\\"



df = pd.read_csv(input_dir + infile + ".tsv", sep="\t")
manta_df = df[df['CALLERS'].str.contains("manta")]
delly_df = df[df['CALLERS'].str.contains("delly")]
gridss_df = df[df['CALLERS'].str.contains("gridss")]
lumpy_df = df[df['CALLERS'].str.contains("lumpy")]
svaba_df = df[df['CALLERS'].str.contains("svaba")]


pd.DataFrame.to_csv(manta_df, output_dir + outfile + 'manta.tsv', sep='\t', header = False, na_rep='.', index=False)
pd.DataFrame.to_csv(delly_df, output_dir + outfile + 'delly.tsv', sep='\t', header = False, na_rep='.', index=False)
pd.DataFrame.to_csv(gridss_df, output_dir + outfile + 'gridss.tsv', sep='\t', header = False, na_rep='.', index=False)
pd.DataFrame.to_csv(lumpy_df, output_dir + outfile + 'lumpy.tsv', sep='\t', header = False, na_rep='.', index=False)
pd.DataFrame.to_csv(svaba_df, output_dir + outfile + 'svaba.tsv', sep='\t', header = False, na_rep='.', index=False)
