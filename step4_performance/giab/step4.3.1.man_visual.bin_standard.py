import os
import pandas as pd
import matplotlib.pyplot as plt  
import seaborn as sns 


input_dir = "...\\combination_sv\\step4_performance\\giab\\stat\\"

output_dir = "...\\combination_sv\\step4_performance\\giab\\stat\\"
if not os.path.exists(output_dir):
    os.mkdir(output_dir)


inputs = []
for input in os.listdir(input_dir):
    if input.lower()[:].find(".data.csv") != -1:
        input_path = input_dir + input
        word = input.split('.')
        output = word[0]
        output_path = output_dir + output
        inputs.append(input_path)
print(inputs)
print(output_path)

for i in range(len(inputs)):
    open_input = inputs[i]

    df = pd.read_csv(open_input)
    df['SVLEN'] = pd.to_numeric(df['SVLEN'], errors='coerce')
    print(df.dtypes)
    df['SVLEN'] = df['SVLEN'].abs()


#set up bins 
bin = [50, 1000, 10000] # 49 < SVLEN < 10001
# bin = [0, 50, 1000, 10000, 100000, 1000000]  ## not filtered
#use pd.cut function can attribute the values into its specific bins 
category = pd.cut(df.SVLEN, bin, right=False) 
category = category.to_frame() 
category.columns = ['SV Length'] 
#concatenate age and its bin 
df_new = pd.concat([df,category],axis = 1) 


#draw histogram plot 
ax = sns.countplot(x = 'SV Length', data = df_new, hue = 'SVTYPE', hue_order=['DEL','INS','DUP','INV'])
# ax = sns.countplot(x = 'SV Length', data = df_new, palette = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728'], hue = 'SVTYPE', 
#             hue_order=['DEL','INS','DUP','INV'])
ax.set(title="NIST_Tier1_v0.6, Structural Variants in the Range of Sizes, from 50 to 10K Nucleotides")

for label in ax.containers:
    ax.bar_label(label)
        
plt.show()
plt.savefig(output_path + ".standardrange.png")
# plt.savefig(output_path + ".benchmark.png")



