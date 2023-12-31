import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

input_dir = "...\\combination_sv\\step4_performance\\giab\\stat\\"

output_dir = "...\\combination_sv\\step4_performance\\giab\\stat\\"
if not os.path.exists(output_dir):
    os.mkdir(output_dir)
inputs = []
for input in os.listdir(input_dir):
    if input.lower()[:].find("summary.csv") != -1:
        input_path = input_dir + input
        word = input.split('.')
        output = word[0]
        output_path = output_dir + output
        inputs.append(input_path)
print(inputs)
# print(output_path)

for i in range(len(inputs)):
    open_input = inputs[i]
    # write_output = outputs[i]
    # print(open_input)
    # print(write_output)
    
data_set = pd.read_csv(open_input)
df = pd.DataFrame(data_set)
svtype_cols = [6,7,8,9,10,11]
sv_cols = [1,2,3,4,5]
df_svtype = df[df.columns[svtype_cols]]
df_sv = df[df.columns[sv_cols]]
# print(df)

x = np.arange(len(df))  # the label locations
# print(x)

###############################

width = 0.125  # the width of the bars
multiplier = 0

fig, ax = plt.subplots(layout='constrained') #(2, 1, sharex=True)

for svtype, count in df_svtype.items():
    offset = width * multiplier
    rects = ax.bar(x + offset, count, width, label=svtype) 
    ax.bar_label(rects, padding=3, rotation=90)
    multiplier += 1

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_title('COUNT by different callsets')
ax.legend(loc='upper right')
ax.set_xlabel("CALLSETS")
ax.set_ylabel("COUNT")
ax.set_ylim(0, 10000)
ax.set_xticks(x + width, df["CALLSET"])


plt.savefig(output_path + "_svtypes.png")
plt.show()


########################

width = 0.125  # the width of the bars
multiplier = 0


fig, ax = plt.subplots(layout='constrained') #(2, 1, sharex=True)

for sv_properties, count in df_sv.items():
    offset = width * multiplier
    rects = ax.bar(x + offset, count, width, label=sv_properties) #
    ax.bar_label(rects, padding=3, rotation=90)
    multiplier += 1

# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_title('COUNT by different callsets')
ax.legend(loc='upper right')
ax.set_xlabel("CALLSETS")
ax.set_ylabel("COUNT")
ax.set_ylim(0, 15000)
ax.set_xticks(x + width, df["CALLSET"])


plt.savefig(output_path + "_svprop.png")
plt.show()
