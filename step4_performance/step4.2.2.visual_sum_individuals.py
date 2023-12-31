import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


# folder = "HG002.annotated.refine"
# folder = "HG002.annotated.refine.replaced"
# folder = "HG002.annotated.refine.replaced.filtered"  ## 篩選svlen>49
folder = "HG002.annotated.refine.replaced.filteredrange"  ## 篩選49 < svlen < 10001


input_dir = "...\\combination_sv\\sv_man_package500\\svlen\\" + folder + "\\stat\\"

output_dir = "...\\combination_sv\\sv_man_package500\\svlen\\" + folder + "\\stat\\"
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

fig, ax = plt.subplots(layout='constrained') #(2, 1, sharex=True)

for svtype, count in df_svtype.items():
    offset = width * multiplier
    rects = ax.bar(x + offset, count, width, label=svtype) 
    ax.bar_label(rects, padding=3, rotation=90)
    multiplier += 1

# Add some text for labels, title and custom x-axis tick labels, etc.
if folder.endswith("replaced"):
    ax.set_title('Structural Variants of All Sizes')
    ax.legend(loc='upper right')
    ax.set_xlabel("CALLERS")
    ax.set_ylabel("Count")
    ax.set_ylim(0, 30000)
if folder.endswith("filteredrange"):
    ax.set_title('Structural Variants in the Range of Sizes, from 50 to 10K Nucleotides')
    ax.legend(loc='upper right')
    ax.set_xlabel("CALLERS")
    ax.set_ylabel("COUNT")
    ax.set_ylim(0, 6000)


ax.set_xticks(x + width, df["CALLSET"])

plt.savefig(write_output + ".svtype.png")
plt.show()


