import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

folder = "HG002.annotated.refine.replaced.filteredrange"  ## 49 < svlen < 10001


input_dir = "...\\combination_sv\\step4_performance\\" + folder + "\\performance\\"

output_dir = "...\\combination_sv\\step4_performance\\" + folder + "\\performance\\"
if not os.path.exists(output_dir):
    os.mkdir(output_dir)



inputs = []
outputs = []
fig_titles = []
for input in os.listdir(input_dir):
    if input.lower()[:].find("performance_man.csv") != -1:
        input_path = input_dir + input
        word = input.split('.')
        output = word[0] + "." + word[1] + "." + word[-2]
        output_path = output_dir + output
        fig_title = word[0] + "_" + word[1]
        inputs.append(input_path)
        outputs.append(output_path)
        fig_titles.append(fig_title)
print(inputs)
print(output_path)
print(fig_titles)

for i in range(len(inputs)):
    open_input = inputs[i]
    # write_output = outputs[i]
    # print(open_input)
    # print(write_output)

data_set = pd.read_csv(open_input)
df = pd.DataFrame(data_set)

prec_recall_col = [10,11] # DEL
# prec_recall_col = [16,17] # INS
# prec_recall_col = [10,11,16,17]
f1_col = [12] # DEL
# f1_col = [18] # INS
# f1_col = [12,18] 
svtype_cols = [2,5] # DEL
# svtype_cols = [3,6] # INS
# svtype_cols = [4,5,6] 
performance_cols = [7,8,9] # DEL
# performance_cols = [13,14,15] # INS
# performance_cols = [7,8,9,13,14,15]

df_prec_recall = df[df.columns[prec_recall_col]]
df_f1 = df[df.columns[f1_col]]
df_svtype = df[df.columns[svtype_cols]]
df_performance = df[df.columns[performance_cols]]
# print(df_prec_recall)
# print(df_f1)
# print(df_svtype)
# print(df_performance)

x = np.arange(len(df))  # the label locations
# print(x)

print(df["CALLERS"])

##############################

width = 0.125  # the width of the bars
multiplier = 0

fig, ax1 = plt.subplots(layout='constrained') #(2, 1, sharex=True)
ax2 = ax1.twinx()

for callset, prec_recall in df_prec_recall.items():
    offset = width * multiplier
    rects = ax1.bar(x + offset, prec_recall, width, label=callset) 
    ax1.bar_label(rects, padding=3, rotation=90)
    ax2.bar_label(rects, padding=3, rotation=90)
    multiplier += 1
    
color = ["green"]
for callset, f1 in df_f1.items():
    offset = width * multiplier
    rects2 = ax2.bar(x + offset, f1, width, label=callset, color = color) 
    ax1.bar_label(rects, padding=3, rotation=90)
    ax2.bar_label(rects2, padding=3, rotation=90)
    multiplier += 1

# Add some text for labels, title and custom x-axis tick labels, etc.
ax1.set_title('Precision, Recall and F1 score of Single Callers')
ax1.legend(loc='upper left')
ax2.legend(loc='upper right')
ax1.set_ylabel("Precision and Recall (%)")
ax2.set_ylabel("F1 Score")
ax1.set_xlabel("Single Caller")
ax1.set_xticks(x + width, df["CALLERS"])
ax1.set_ylim(0, 100)
ax2.set_ylim(0.0, 1.0)
plt.savefig(output_path + ".PRF_man.png")
plt.show()


##############################

width = 0.25  # the width of the bars
multiplier = 0

fig, ax = plt.subplots(layout='constrained') #(2, 1, sharex=True)

for callset, f1 in df_f1.items():
    offset = width * multiplier
    rects = ax.bar(x + offset, f1, width, label=callset) 
    ax.bar_label(rects, padding=3, rotation=90)
    multiplier += 1

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_title('F1 Score of Single Callers')
ax.legend(loc='upper right')
ax.set_ylabel("SCORE")
ax.set_xlabel("Single Caller")
ax.set_xticks(x, df["CALLERS"])
ax.set_ylim(0, 1)

# plt.savefig(output_path + ".f1_man.png")
plt.show()

##############################

width = 0.125  # the width of the bars
multiplier = 0

fig, ax = plt.subplots(layout='constrained') #(2, 1, sharex=True)


for callset, count in df_svtype.items():
    offset = width * multiplier
    rects = ax.bar(x + offset, count, width, label=callset) 
    ax.bar_label(rects, padding=3, rotation=90)
    multiplier += 1


# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_title('COUNT of Single Callers')
ax.legend(loc='upper right')
ax.set_ylabel("COUNT")
ax.set_xlabel("Single Caller")
ax.set_xticks(x + width/2, df["CALLERS"])
ax.set_ylim(0, 7000)

# plt.savefig(output_path + ".testcount_man.png")
plt.show()

######################################

width = 0.125  # the width of the bars
multiplier = 0

fig, ax = plt.subplots(layout='constrained') #(2, 1, sharex=True)


for callset, count in df_performance.items():
    offset = width * multiplier
    rects = ax.bar(x + offset, count, width, label=callset) 
    ax.bar_label(rects, padding=3, rotation=90)
    multiplier += 1

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_title('True Positive, False Positive and False Negative of Single Callers')
ax.legend(loc='upper right', ncol=6)
ax.set_ylabel("COUNT")
ax.set_xlabel("Single Caller")
ax.set_xticks(x + width, df["CALLERS"])
ax.set_ylim(0, 6000)

# plt.savefig(output_path + ".fpfn_man.png")
plt.show()

##############################
