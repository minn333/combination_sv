import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

input_dir = "...\\combination_sv\\step4_performance\\models\\for_intv\\"

output_dir = "...\\combination_sv\\step4_performance\\models\\for_intv\\"
if not os.path.exists(output_dir):
    os.mkdir(output_dir)

df_sigsv_mdg = pd.read_table(input_dir + "intersect.HG002.mdg.intv.tsv", sep='\t', low_memory=False)
# print(df)
data_1 = df_sigsv_mdg["INTERVAL"]
new_data_1 = [item for item in data_1 if not(pd.isnull(item)) == True]
# print(new_data_1)

df_sigsv_nocnv = pd.read_table(input_dir + "intersect.HG002.mdgls.intv.tsv", sep='\t', low_memory=False)
data_2 = df_sigsv_nocnv["INTERVAL"]
new_data_2 = [item for item in data_2 if not(pd.isnull(item)) == True]

df_sigsv0_mdg = pd.read_table(input_dir + "union.HG002.mdg.intv.tsv", sep='\t', low_memory=False)
data_3 = df_sigsv0_mdg["INTERVAL"]
new_data_3 = [item for item in data_3 if not(pd.isnull(item)) == True]

df_sigsv0_nocnv = pd.read_table(input_dir + "union.HG002.mdgls.intv.tsv", sep='\t', low_memory=False)
data_4 = df_sigsv0_nocnv["INTERVAL"]
new_data_4 = [item for item in data_4 if not(pd.isnull(item)) == True]

data = [new_data_1, new_data_2, new_data_3, new_data_4]


fig = plt.figure(figsize =(10, 7))
ax = fig.add_subplot(111)


# Creating axes instance
bp = ax.boxplot(data, patch_artist = True,
                notch ='True', vert = 90)

colors = ['#0000FF', '#00FF00',
          '#FFFF00', '#FF00FF']

for patch, color in zip(bp['boxes'], colors):
    patch.set_facecolor(color)


# changing color and linewidth of
# whiskers
for whisker in bp['whiskers']:
    whisker.set(color ='#8B008B',
                linewidth = 1.5,
                linestyle =":")

# changing color and linewidth of
# caps
for cap in bp['caps']:
    cap.set(color ='#8B008B',
            linewidth = 2)
    

# changing color and linewidth of
# medians
for median in bp['medians']:
    median.set(color ='red',
               linewidth = 3)

# changing style of fliers
for flier in bp['fliers']:
    flier.set(marker ='D',
              color ='#e7298a',
              alpha = 0.5)

# # x-axis labels
ax.set_xticklabels(['III-interesection', 'V-intersection',
                    'III-union', 'V-union'])

# # Adding title
plt.title("Distances of the Neighbor SVs")
ax.set_ylabel("Nucleotides")
ax.set_xlabel("Calling Strategy")
plt.suptitle("Distances of the Neighbor SVs")
# ax2.set_ylabel("BASES")
ax.set_ylim(0, 600)
# ax.set_ylim(0, 50)

plt.show()