import os
from sv_function import *

class SV():
    def __init__(self,chr,pos,id,ref,alt,qual,filter,info,format,others):

        self.chr = chr   #field[0]
        self.pos = int(pos) #field[1]
        self.id = id #field[2]
        self.ref = ref #field[3]
        self.alt = alt #field[4]
        self.qual = qual #field[5]
        self.filter = filter #field[6]
        self.info = info #field[7]
        self.format = format #field[8]
        self.others = others #field[9]
    
    def __str__(self):
        return self.chr + "\t" + str(self.pos) + "\t" +self.id + "\t" +self.ref + "\t" +self.alt + "\t" +self.qual + "\t" + \
               self.filter + "\t" + self.info + "\t" + self.format + "\t" + self.others

chr_list = ["1","2", "3", "4", "5", "6","7","8", "9", "10", "11","12","13", "14", "15", "16",
            "17", "18", "19", "20", "21", "22", "X", "Y"]

input_dir = "...\\combination_sv\\step4_performance\\giab\\"

output_dir = "...\\combination_sv\\step4_performance\\giab\\"
if not os.path.exists(output_dir):
    os.mkdir(output_dir)


inputs = []
outputs = []
for input in os.listdir(input_dir):
    if input.lower()[:].find("benchmarkcalls.vcf") != -1:
        input_path = input_dir + input
        word = input.split('.')
        output = word[0] + word[1] + ".refine.standardrange." + word[-1]
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
                info_field = line.split("\t")[7]
                id_field = line.split("\t")[2]
                info_value = info_field.split(";")
                if fields[0] in chr_list: # NISTv0.6
                    if info_field.find("SVLEN") !=-1:
                        for x in info_value:
                            if "SVLEN" in x:
                                svlen_value = abs(int(x.split("=")[1].strip()))
                                # print(svlen_value)
                                if (49 < svlen_value < 10001): 
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