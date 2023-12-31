import os
import pandas as pd
import matplotlib.pyplot as plt  
import seaborn as sns 


# folder = "HG002.annotated.refine"
# folder = "HG002.annotated.refine.replaced"
# folder = "HG002.annotated.refine.replaced.filtered"  ## 篩選svlen>49
folder = "HG002.annotated.refine.replaced.filteredrange"  ## 篩選49 < svlen < 10000


input_dir = "...\\combination_sv\\step4_performance\\" + folder + "\\stat\\"

output_dir = "...\\combination_sv\\step4_performance\\" + folder + "\\stat\\"
if not os.path.exists(output_dir):
    os.mkdir(output_dir)


# fig_titles = list(["DELLY", "DRAGEN", "GRIDSS", "LUMPY", "Manta", "SvABA"])
inputs = []
outputs = []
fig_titles = []
for input in os.listdir(input_dir):
    if input.lower()[:].find("data.csv") != -1:
        input_path = input_dir + input
        word = input.split('.')
        output = word[0] + "." + word[1] + "." + word[-3]
        fig_title = word[-3] + "_" + word[1]
        output_path = output_dir + output
        inputs.append(input_path)
        outputs.append(output_path)
        fig_titles.append(fig_title)
# print(inputs)
# print(outputs)
# print(fig_titles)

new_fig_titles = list(["DELLY", "DRAGEN", "GRIDSS", "LUMPY", "Manta", "SvABA"])

dfs = []
for file in inputs:
    df = pd.read_csv(file)
    dfs.append(df)

# print(dfs)
    # print(df.shape)

for i, df in enumerate(dfs):
    df['SVLEN'] = pd.to_numeric(df['SVLEN'], errors='coerce')
    # print(df.dtypes)
    df['SVLEN'] = df['SVLEN'].abs()

    #set up bins 
    bin = [50, 1000, 10000] # filtered ranged
    # bin = [0, 50, 1000, 10000, 100000, 1000000]  ## not filtered
    #use pd.cut function can attribute the values into its specific bins 
    category = pd.cut(df.SVLEN, bin, right=False) 
    category = category.to_frame() 
    category.columns = ['SV Length'] 

    #concatenate age and its bin 
    df_new = pd.concat([df,category],axis = 1) 
    # print(df_new)


    #draw histogram plot 
    ax = sns.countplot(x = 'SV Length', data = df_new, hue = 'SVTYPE', hue_order=['DEL','INS','DUP','INV'])
    # ax = sns.countplot(x = 'SV Length', data = df_new, palette = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728'], hue = 'SVTYPE', 
    #             hue_order=['DEL','INS','DUP','INV'])
    ax.set(title=new_fig_titles[i])

    
    for label in ax.containers:
        ax.bar_label(label)

    plt.show()
    plt.savefig(outputs[i] + ".png")
    