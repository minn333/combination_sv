import os
from sv import *


input_dir = "...\\combination_sv\\step1_refine_sv\\HG002.annotated.refine\\"

output_dir = "...\\combination_sv\\step4_performance\\HG002.annotated.refine.replaced.filteredrange\\"
if not os.path.exists(output_dir):
    os.mkdir(output_dir)
    
keyword = "dragen"
# keyword = "manta"
# keyword = "delly"
# keyword = "gridss"
# keyword = "lumpy"
# keyword = "svaba"
    
inputs = []
outputs = []
for input in os.listdir(input_dir):
    if input.lower()[:].find(keyword + ".vcf") != -1:
        input_path = input_dir + input
        word = input.split('.')
        output = word[0] + "." + word[1] + "." + word[2] + ".replaced.filteredrange." + word[-2] + "." + word[-1]
        output_path = output_dir + output
        inputs.append(input_path)
        outputs.append(output_path)
# print(inputs)
# print(outputs)

for i in range(len(inputs)):
    open_input = inputs[i]
    write_output = outputs[i]
    # print(open_input)
    print(write_output)
    
    with open(open_input,"r") as input_vcf:
        headers = []
        svs=[]
        for line in input_vcf:
            if(line[0]=="#" or line[0:2]=="##"):
                headers.append(line)
            if(line[0]!="#" and line[0:2]!="##"):
                fields = line.split("\t")
                info_field = line.split("\t")[7]
                id_field = line.split("\t")[2]
                info_value = info_field.split(";")

### Except Dragen
                if(info_field.find("SVTYPE") !=-1 and info_field.find("SIMPLE_TYPE") != -1):
                    oldstring = ""
                    newstring = ""
                    for x in info_value:
                        if "SVTYPE" in x:
                            oldstring = x.split("=")[1].strip()
                        if "SIMPLE_TYPE" in x:
                            newstring = x.split("=")[1].strip()
                            # print(newstring)
                    fields[7] = fields[7].replace(oldstring,newstring)    
                    
                # print(fields[7]) 
                
                if info_field.find("SVLEN") !=-1:
                    for x in info_value:
                        if "SVLEN" in x:
                            svlen_value = abs(int(x.split("=")[1].strip()))
                            # print(svlen_value)
                            if svlen_value > 49 and svlen_value < 10001:                             
                                sv = write_vcf(fields[0],fields[1],fields[2],fields[3],fields[4],\
                                    fields[5],fields[6],fields[7],fields[8],fields[9])
                                svs.append(sv)
            # print(svs)

# # For Dragen
#                 if info_field.find("SVLEN") !=-1:
#                     for x in info_value:
#                         if "SVLEN" in x:
#                             svlen_value = abs(int(x.split("=")[1].strip()))
#                             # print(svlen_value)
#                             if svlen_value > 49 and svlen_value < 10001:                                
#                                 sv = write_vcf(fields[0],fields[1],fields[2],fields[3],fields[4],\
#                                     fields[5],fields[6],fields[7],fields[8],fields[9])
#                                 svs.append(sv)

        svs.sort(key=lambda x: (len(x.chr), x.chr, x.pos))
        
    with open(write_output, "w") as replaced_f:
        for header in headers:
            replaced_f.write(header.strip() + "\n")
        for m in svs:
            replaced_f.write(str(m).strip() + "\n") 