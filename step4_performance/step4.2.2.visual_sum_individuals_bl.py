import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


# folder = "HG002.annotated.refine"
# folder = "HG002.annotated.refine.replaced"
# folder = "HG002.annotated.refine.replaced.filtered"  ## 篩選svlen>49
folder = "HG002.annotated.refine.replaced.filteredrange"  ## 篩選49 < svlen < 10001


input_dir = "...\\combination_sv\\step4_performance\\" + folder + "\\stat\\"

output_dir = "...\\combination_sv\\step4_performance\\" + folder + "\\stat\\"
if not os.path.exists(output_dir):
    os.mkdir(output_dir)


inputs = []
outputs = []
fig_titles = []
for input in os.listdir(input_dir):
    if input.lower()[:].find("summary_man.csv") != -1:
        input_path = input_dir + input
        word = input.split('.')
        output = word[0] + "." + word[1]
        fig_title = word[0] + "_" + word[1]
        output_path = output_dir + output
        inputs.append(input_path)
        outputs.append(output_path)
        fig_titles.append(fig_title)
print(inputs)
print(outputs)

for i in range(len(inputs)):
    open_input = inputs[i]
    write_output = outputs[i]
    print(open_input)
    print(write_output)
    
data_set = pd.read_csv(open_input)
df = pd.DataFrame(data_set)
svtype_cols = [6,7,8,10,9]
sv_cols = [1,2,3,4]
df_svtype = df[df.columns[svtype_cols]]
df_sv = df[df.columns[sv_cols]]
# print(df)

x = np.arange(len(df))  # the label locations
# print(x)


########################################

width = 0.125  # the width of the bars
multiplier = 0

# If we were to simply plot pts, we'd lose most of the interesting
# details due to the outliers. So let's 'break' or 'cut-out' the y-axis
# into two portions - use the top (ax1) for the outliers, and the bottom
# (ax2) for the details of the majority of our data
fig, (ax1, ax2)= plt.subplots(nrows=2, ncols=1, sharex=True) #(2, 1, sharex=True)
fig.subplots_adjust(hspace=0.05) # adjust space between axes


# zoom-in / limit the view to different portions of the data
if folder.endswith("replaced"):
    ax2.set_ylim(0, 11000) # most of the data
    ax1.set_ylim(20000, 30000) # outliers only
if folder.endswith("filteredrange"):
    ax2.set_ylim(0, 600) # most of the data
    ax1.set_ylim(700, 6000) # outliers only


# fig, ax = plt.subplots(layout='constrained') #(2, 1, sharex=True)

for svtype, count in df_svtype.items():
    offset = width * multiplier
    rects = ax2.bar(x + offset, count, width, label=svtype) 
    rects = ax1.bar(x + offset, count, width, label=svtype) 
    ax2.bar_label(rects, padding=3, rotation=90)
    ax1.bar_label(rects, padding=3, rotation=90)
    multiplier += 1

# Add some text for labels, title and custom x-axis tick labels, etc.
# hide the spines between ax and ax2
ax1.spines.bottom.set_visible(False)
ax2.spines.top.set_visible(False)
ax1.xaxis.set_ticks_position("none")
# ax2.xaxis.set_ticks_position("none")
ax2.tick_params(labeltop=False)  # don't put tick labels at the top
ax2.xaxis.tick_bottom()


# ax2.set_ylabel("BASES")


if folder.endswith("replaced"):
    plt.suptitle("Structural Variants of All Sizes")
    ax1.legend(loc='upper right')
    ax1.set_ylabel("Count", loc='bottom')
    ax2.set_xlabel("CALLERS")
if folder.endswith("filteredrange"):
    plt.suptitle('Structural Variants in the Range of Sizes, from 50 to 10K Nucleotides')
    ax1.legend(loc='upper right')
    ax1.set_ylabel("Count", loc='bottom')
    ax2.set_xlabel("CALLERS")


ax2.set_xticks(x + width, df["CALLSET"])

plt.savefig(write_output + ".svtype.png")
plt.show()


