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


# fig = plt.figure(figsize =(10, 7))
# ax = fig.add_subplot(111)

# If we were to simply plot pts, we'd lose most of the interesting
# details due to the outliers. So let's 'break' or 'cut-out' the y-axis
# into two portions - use the top (ax1) for the outliers, and the bottom
# (ax2) for the details of the majority of our data
fig, (ax1, ax2)= plt.subplots(nrows=2, ncols=1, sharex=True) #(2, 1, sharex=True)
fig.subplots_adjust(hspace=0.05) # adjust space between axes


# zoom-in / limit the view to different portions of the data
ax2.set_ylim(0, 50) # most of the data
ax1.set_ylim(100, 600) # outliers only


# Creating axes instance
bp1 = ax1.boxplot(data, patch_artist = True,
                notch ='True', vert = 90)
bp2 = ax2.boxplot(data, patch_artist = True,
                notch ='True', vert = 90)
 
colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
 
for patch, color in zip(bp1['boxes'], colors):
    patch.set_facecolor(color)

for patch, color in zip(bp2['boxes'], colors):
    patch.set_facecolor(color)

# changing color and linewidth of
# whiskers
for whisker in bp1['whiskers']:
    whisker.set(color ='#8B008B',
                linewidth = 1.5,
                linestyle =":")
    
for whisker in bp2['whiskers']:
    whisker.set(color ='#8B008B',
                linewidth = 1.5,
                linestyle =":")
 
# changing color and linewidth of
# caps
for cap in bp1['caps']:
    cap.set(color ='#8B008B',
            linewidth = 2)
    
for cap in bp2['caps']:
    cap.set(color ='#8B008B',
            linewidth = 2)
 
# changing color and linewidth of
# medians
for median in bp1['medians']:
    median.set(color ='yellow',
               linewidth = 3)

for median in bp2['medians']:
    median.set(color ='yellow',
               linewidth = 3)
 
# changing style of fliers
for flier in bp1['fliers']:
    flier.set(marker ='D',
              color ='#e7298a',
              alpha = 0.5)

for flier in bp2['fliers']:
    flier.set(marker ='D',
              color ='#e7298a',
              alpha = 0.5)


# hide the spines between ax and ax2
ax1.spines.bottom.set_visible(False)
ax2.spines.top.set_visible(False)
ax1.xaxis.set_ticks_position("none")
# ax2.xaxis.set_ticks_position("none")
ax2.tick_params(labeltop=False)  # don't put tick labels at the top
ax2.xaxis.tick_bottom()
     
# # x-axis labels
ax1.set_xticklabels(['1','2' , '3', '4',  'III-interesection', 'V-intersection',
                    'III-union', 'V-union'])

# # Adding title
plt.suptitle("Distances of the Neighbor SVs")
# ax2.set_ylabel("BASES")
ax1.set_ylabel("Nucleotides", loc='bottom')
ax2.set_xlabel("Calling Strategy")
# # ax.set_xlim(0, 600)
# # ax.set_xlim(0, 40)
 
# # Removing top axes and right axes
# # ticks
# ax1.get_xaxis().tick_bottom()
# ax1.get_yaxis().tick_left()
     
# show plot
plt.savefig(output_dir + "models.intv_man.png")
plt.show()


def get_summary_statistics(dataset):
    
    mean = np.round(np.mean(dataset), 2)
    median = np.round(np.median(dataset), 2)
    min_value = np.round(dataset.min(), 2)
    max_value = np.round(dataset.max(), 2)
    quartile_1 = np.round(dataset.quantile(0.25), 2)
    quartile_3 = np.round(dataset.quantile(0.75), 2)
    # Interquartile range
    iqr = np.round(quartile_3 - quartile_1, 2)
    print('Min: %s' % min_value)
    print('Mean: %s' % mean)
    print('Max: %s' % max_value)
    print('25th percentile: %s' % quartile_1)
    print('Median: %s' % median)
    print('75th percentile: %s' % quartile_3)
    print('Interquartile range (IQR): %s' % iqr)
    # print('Setosa summary statistics')
    
print('\n\nSummary statistics (intersect.HG002.mdg)')
get_summary_statistics(data_1)
print('\n\nSummary statistics (intersect.HG002.mdgls)')
get_summary_statistics(data_2)
print('\n\nSummary statistics (union.HG002.mdg)')
get_summary_statistics(data_3)
print('\n\nSummary statistics (union.HG002.mdgls)')
get_summary_statistics(data_4)