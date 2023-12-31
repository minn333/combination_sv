import os
import pandas as pd
import matplotlib.pyplot as plt  
import seaborn as sns 

color = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2',
 '#7f7f7f', '#bcbd22', '#17becf']

folder = "models\\"


input_dir = "...\\combination_sv\\step4_performance\\" + folder + "stat\\"

output_dir = "...\\combination_sv\\step4_performance\\" + folder + "stat\\"
if not os.path.exists(output_dir):
    os.mkdir(output_dir)


inputs = []
outputs = []
fig_titles = []
for input in os.listdir(input_dir):
    if input.lower()[:].find("data.csv") != -1:
        input_path = input_dir + input
        word = input.split('.')
        output = word[0] + "." + word[1] + "." + word[-3]
        fig_title = word[-3] + "_range_" + word[1]
        output_path = output_dir + output
        inputs.append(input_path)
        outputs.append(output_path)
        fig_titles.append(fig_title)
# print(inputs)
print(outputs)
print(fig_titles)


new_fig_titles = list(["DRAGEN", "III-intersect", "V-intersect" , "III-union", "V-union"])

dfs = []
for file in inputs:
    df = pd.read_csv(file)
    dfs.append(df)

# print(dfs)
    # print(df.shape)

for i, df in enumerate(dfs):
    df['SVLEN'] = pd.to_numeric(df['SVLEN'], errors='coerce')
    print(df.dtypes)
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


    #draw histogram plot 
    ax = sns.countplot(x = 'SV Length', data = df_new, hue = 'SVTYPE', hue_order=['DEL','INS','DUP','INV'])
    # ax = sns.countplot(x = 'SV Length', data = df_new, palette = 'hls', hue = 'SVTYPE', hue_order=['DEL','INS','DUP','INV'])
    ax.set(title=new_fig_titles[i])
    
    for label in ax.containers:
        ax.bar_label(label)
    plt.show()
    plt.savefig(outputs[i] + ".png")
