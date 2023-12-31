import os
from sv_function import *


chr_list = ["chr1","chr2","chr3", "chr4", "chr5", "chr6","chr7","chr8", "chr9", "chr10",
    "chr11","chr12","chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20",
    "chr21", "chr22", "chrX", "chrY"]

input_dir = "...\\combinaton_sv\\step1_refine_sv\\HG002.annotated\\"

output_dir = "...\\combinaton_sv\\step1_refine_sv\\HG002.annotated.refine\\"
if not os.path.exists(output_dir):
    os.mkdir(output_dir)

keyword = "gridss"
    
inputs = []
outputs = []
for input in os.listdir(input_dir):
    if input.lower()[:].find(keyword + ".vcf") != -1:
        input_path = input_dir + input
        word = input.split('.')
        output = word[0] + "." + word[1] + ".refine." + word[-2] + "." + word[-1]
        output_path = output_dir + output
        inputs.append(input_path)
        outputs.append(output_path)
print(inputs)
print(outputs)

for i in range(len(inputs)):
    open_input = inputs[i]
    write_output = outputs[i]
    # print(open_input)
    # print(write_output)
    
    with open(open_input, "r") as input_vcf:
        headers = []
        svs=[]
        for line in input_vcf:
            if(line[0]=="#" or line[0:2]=="##"):
                headers.append(line)
            if(line[0]!="#" and line[0:2]!="##"):
                fields = line.split("\t")
                if(fields[0] in chr_list and fields[6].find("PASS") !=-1) and not fields[2].endswith("h"): # gridss
                # if(fields[0] in chr_list and fields[6].find("PASS") !=-1) and not fields[2].endswith(":2"): # svaba
                # if(fields[0] in chr_list and not fields[2].endswith("_2")): # lumpy
                # if(fields[0] in chr_list and fields[6].find("PASS") !=-1): # delly, cnvpytor
                # if(fields[0] in chr_list and fields[6].find("PASS") !=-1) and not fields[2].endswith(":1"): # manta, dragen
                    sv = write_vcf(fields[0],fields[1],fields[2],fields[3],fields[4],\
                    fields[5],fields[6],fields[7],fields[8],fields[9])
                    svs.append(sv)
        # print(svs)

        # svs.sort(key=lambda x: (len(x.chr), x.chr, x.pos))
        
    with open(write_output, "w") as filtered_f:
        for header in headers:
            filtered_f.write(header.strip() + "\n")
        for m in svs:
            filtered_f.write(str(m).strip() + "\n") 