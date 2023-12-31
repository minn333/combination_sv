import os
import pandas as pd
from sv_function import *


folder = "models\\"


input_dir = "...\\combination_sv\\step4_performance\\" + folder

output_dir = "...\\combination_sv\\step4_performance\\" + folder + "stat\\"
if not os.path.exists(output_dir):
    os.mkdir(output_dir)



callsets = []
inputs = []
outputs = []
for input in os.listdir(input_dir):
    if input.lower()[:].find(".vcf") != -1:
        input_path = input_dir + input
        word = input.split('.')
        output = word[0] + "." + word[1] + "." + word[-2]
        output_path = output_dir + output
        count_file_path = output_dir + word[0]
        callset = word[1] + "_" + word[-2]
        callsets.append(callset)
        inputs.append(input_path)
        outputs.append(output_path)

# print(inputs)
# print(outputs)
# print(callsets)
# print(count_file_path)

svcount_list = []; svlen_count_list = []; 
GT_count_list = []; svtype_count_list = []; simple_type_count_list = []
dels_list = []; inss_list = []; dups_list = []; ctxs_list = []; invs_list = []; 
bnds_list = []; unspecified_list = []


for i in range(len(inputs)):
    open_input = inputs[i]
    write_output = outputs[i]
    # print(open_input)
    # print(write_output)

    svcallers = []; svlens = []; svtypes = []; simple_types = []; GT_values = []
    with open(open_input ,"r") as input_vcf:
        sv_count = 0; caller_count= 0; svlen_count = 0; GT_count = 0; svtype_count = 0; simple_type_count = 0
        del_count = 0; ins_count = 0; dup_count = 0; ctx_count = 0; inv_count = 0; bnd_count = 0; 
        for line in input_vcf:
            if(line[0]!="#" and line[0:2]!="##"):
                sv_count += 1 
                fields = line.split("\t")
                info_field = fields[7]
                format_field = fields[8]
                info_value = info_field.split(";")
                format_value = format_field.split(":")
                
                if(info_field.find("SVLEN")!=-1):
                    svlen_count += 1
                if(info_field.find("SVTYPE")!=-1):
                    svtype_count += 1
                    for x in info_value:
                        if "SVTYPE" in x:
                            svtype = x.split("=")[1].strip()
                            if svtype == "DEL":
                                del_count += 1
                            if svtype == "INS":
                                ins_count += 1
                            if svtype == "DUP":
                                dup_count += 1
                            if svtype == "CTX":
                                ctx_count += 1
                            if svtype == "INV":
                                inv_count += 1
                            if svtype == "BND":
                                bnd_count += 1

                if(info_field.find("SIMPLE_TYPE")!=-1):
                    simple_type_count += 1
                    
                if(format_field.find("GT")!=-1):
                    GT_count += 1
                
                svcaller = findInfo(info_value,"CALLER")
                svcallers.append(svcaller)
                svlen = findInfo(info_value,"SVLEN")
                svlens.append(svlen)
                svtype = findInfo(info_value,"SVTYPE")
                svtypes.append(svtype)
                simple_type = findInfo(info_value,"SIMPLE_TYPE")
                simple_types.append(simple_type)
                GT_value = findFormat(format_field, fields[9])
                GT_values.append(GT_value)

    svcount_list.append(sv_count), svlen_count_list.append(svlen_count), GT_count_list.append(GT_count),
    svtype_count_list.append(svtype_count); simple_type_count_list.append(simple_type_count)
    dels_list.append(del_count), inss_list.append(ins_count), dups_list.append(dup_count), 
    ctxs_list.append(ctx_count), invs_list.append(inv_count), bnds_list.append(bnd_count)
    
    data = {"CALLER":svcallers, "SVLEN":svlens, "SVTYPE":svtypes, "SIMPLE_TYPE":simple_types, "GT": GT_values}
    df = pd.DataFrame(data)
    df.to_csv(write_output + ".data.csv", index=False)

df = pd.DataFrame({"CALLSET":callsets, "SV_COUNT":svcount_list, "SVLEN_COUNT":svlen_count_list,"GT_COUNT":GT_count_list, 
                   "SVTYPE_COUNT":svtype_count_list, "SIMPLE_TYPE_COUNT":simple_type_count_list,
                   "DEL":dels_list, "INS":inss_list, "DUP":dups_list, 
                   "CTX":ctxs_list, "INV":invs_list, "BND":bnds_list})
df.to_csv(count_file_path + ".summary.csv", index=False)

