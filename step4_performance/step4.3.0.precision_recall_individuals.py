import os    
import pandas as pd   
import matplotlib.pyplot as plt
from sv_function import *


folder = "HG002.annotated.refine.replaced.filteredrange"  ## 篩選49 < svlen < 10001


test_dir = "...\\combination_sv\\step4_performance\\" + folder + "\\"

output_dir = "...\\combination_sv\\step4_performance\\" + folder + "\\performance\\"
if not os.path.exists(output_dir):
    os.mkdir(output_dir)
    

#############################################################################

truth_dir = "...\\combination_sv\\sv_man_package500\\svlen\\standard\\"

output_truth_dir = "...\\combination_sv\\sv_man_package500\\svlen\\standard\\truth\\"
if not os.path.exists(output_truth_dir):
    os.mkdir(output_truth_dir)


truth_count_list = []; DEL_truth_count_list = []; INS_truth_count_list = []
test_count_list = []; DEL_test_count_list = []; INS_test_count_list = []

DEL_tp_list = []; DEL_fp_list = []; DEL_fn_list = []; 
DEL_recall_list=[]; DEL_precision_list = [];  DEL_f1_list = []
INS_tp_list = []; INS_fp_list = []; INS_fn_list = []; 
INS_recall_list=[]; INS_precision_list = [];  INS_f1_list = []


with open(truth_dir + "HG002_SVs_Tier1_v06.refine.standardrange.vcf", "r") as truth_f:  #篩選svlen>49
    del_truth_vcf=[]
    ins_truth_vcf=[]
    truth_count = 0; DEL_truth_count = 0; INS_truth_count = 0
    for line in truth_f.readlines():
        if(line[0]!="#" and line[0:2]!="##"):
            truth_count += 1
            item = line.split('\t')
            info_value = item[7].split(";")
            svtype = findInfo(info_value,"SVTYPE")
            for x in info_value:
                if "SVTYPE" in x:
                    svtype = x.split("=")[1].strip()
                    if svtype == "DEL":
                        del_truth_vcf.append(Variant("chr"+item[0], item[1], svtype))
                        DEL_truth_count += 1
            for x in info_value:
                if "SVTYPE" in x:
                    svtype = x.split("=")[1].strip()
                    if svtype == "INS":
                        ins_truth_vcf.append(Variant("chr"+item[0], item[1], svtype))
                        INS_truth_count += 1

with open(output_truth_dir + "HG002.truthdelrange.vcf", "w") as del_truth_f:
    for m in del_truth_vcf:
        del_truth_f.write(str(m).strip() + "\n") 

with open(output_truth_dir + "HG002.truthinsrange.vcf", "w") as ins_truth_f:
    for m in ins_truth_vcf:
        ins_truth_f.write(str(m).strip() + "\n")



callers = []
inputs = []
outputs = []
for input in os.listdir(test_dir):
    if input.lower()[:].find(".vcf") != -1:
        input_path = test_dir + input
        word = input.split('.')
        output = word[0] + ".testrange"
        caller = word[-2]
        output_path = output_dir + output
        inputs.append(input_path)
        callers.append(caller)
        outputs.append(output_path)
# print(inputs)
# print(callers)
# print(outputs)
# print(output_path)

# ######################


for i in range(len(inputs)):
    open_input = inputs[i]
    write_output = outputs[i]
    # print(open_input)
    # print(write_output)
    
    del_test_vcf = []
    ins_test_vcf = []
    with open(open_input, "r") as test_f:
        test_count= 0; DEL_test_count = 0; INS_test_count = 0
        for line in test_f.readlines():
            if(line[0]!="#" or line[0:2]!="##"):
                test_count += 1
                item = line.split('\t')
                info_value = item[7].split(";")
                svtype = findInfo(info_value,"SVTYPE")
                for x in info_value:
                    if "SVTYPE" in x:
                        svtype = x.split("=")[1].strip()
                        if svtype == "DEL":
                            del_test_vcf.append(Variant(item[0], item[1], svtype))
                            DEL_test_count += 1
                for x in info_value:
                    if "SVTYPE" in x:
                        svtype = x.split("=")[1].strip()
                        if svtype == "INS":
                            ins_test_vcf.append(Variant(item[0], item[1], svtype))
                            INS_test_count += 1

    with open(write_output + "delrange." + callers[i] + ".vcf", "w") as del_test_f:
        for m in del_test_vcf:
            del_test_f.write(str(m).strip() + "\n") 
    
    with open(write_output + "insrange." + callers[i] + ".vcf", "w") as ins_test_f:
        for m in ins_test_vcf:
            ins_test_f.write(str(m).strip() + "\n") 
                            

##########################

    if len(del_test_vcf)!=0:
        del_tp=0
        for v1 in del_test_vcf:
            for v2 in del_truth_vcf:
                if(v1.chr==v2.chr): 
                    if v2.pos > v1.pos-500 and v2.pos < v1.pos+500: 
                        del_tp+=1
                        break
                    

        del_fp = len(del_test_vcf)-del_tp
        del_fn = len(del_truth_vcf)-del_tp
        del_recall = round((del_tp / (del_tp + del_fn))*100,2)
        del_precision = round((del_tp / (del_tp + del_fp))*100,2)
        del_f1 = round(((2*del_precision*del_recall) / (del_precision+del_recall))*0.01, 2)

    else:
        del_tp = 0
        del_fp = 0
        del_fn = 0
        del_recall = 0 
        del_precision = 0
        del_f1 = 0

#######################


    if len(ins_test_vcf)!=0:
        ins_tp=0
        for v1 in ins_test_vcf:
            for v2 in ins_truth_vcf:
                if(v1.chr==v2.chr): 
                    if v2.pos >= v1.pos-500 and v2.pos < v1.pos+500: 
                        ins_tp+=1
                        break

        ins_fp = len(ins_test_vcf)-ins_tp
        ins_fn = len(ins_truth_vcf)-ins_tp
        ins_recall = round((ins_tp / (ins_tp + ins_fn))*100,2)   
        ins_precision = round((ins_tp / (ins_tp + ins_fp))*100,2)
        ins_f1 = round(((2*ins_precision*ins_recall) / (ins_precision+ins_recall))*0.01,2)


    else:
        ins_tp = 0
        ins_fp = 0
        ins_fn = 0
        ins_recall = 0
        ins_precision = 0
        ins_f1 = 0


    truth_count_list.append(truth_count)
    DEL_truth_count_list.append(DEL_truth_count)
    INS_truth_count_list.append(INS_truth_count)
    test_count_list.append(test_count)
    DEL_test_count_list.append(DEL_test_count)
    INS_test_count_list.append(INS_test_count)

    DEL_tp_list.append(del_tp); DEL_fp_list.append(del_fp); DEL_fn_list.append(del_fn)
    DEL_recall_list.append(del_recall); DEL_precision_list.append(del_precision)
    DEL_f1_list.append(del_f1)
    INS_tp_list.append(ins_tp); INS_fp_list.append(ins_fp); INS_fn_list.append(ins_fn)
    INS_recall_list.append(ins_recall); INS_precision_list.append(ins_precision)
    INS_f1_list.append(ins_f1)

    
df = pd.DataFrame({"CALLERS":callers, "TRUTH_COUNT":truth_count_list, 
                   "DEL_truth_COUNT":DEL_truth_count_list,"INS_truth_count":INS_truth_count_list, 
                "TEST_COUNT":test_count_list, "DEL_test_COUNT":DEL_test_count_list, 
                "INS_test_COUNT":INS_test_count_list,
                "DEL_tp":DEL_tp_list,"DEL_fp":DEL_fp_list,"DEL_fn":DEL_fn_list,
                "DEL_recall":DEL_recall_list,"DEL_precision":DEL_precision_list,
                "DEL_F1":DEL_f1_list,
                "INS_tp":INS_tp_list,"INS_fp":INS_fp_list,"INS_fn":INS_fn_list,
                "INS_recall":INS_recall_list,"INS_precision":INS_precision_list,
                "INS_F1":INS_f1_list})

df.to_csv(output_path + ".performance.csv", index=False)

