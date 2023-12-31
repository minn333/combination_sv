import os
from sv_function import *
import sys

inputdir= sys.argv[1]
outputdir= sys.argv[2]
mergefile= sys.argv[3]
firstfile= sys.argv[4]
allfile= sys.argv[5]
statfile= sys.argv[6]

#print(samplename)

vcfdir = "...\\combination_sv\\step3_combi_union\\HG002.annotated.refine.replaced.filteredrange\\" + inputdir + "\\"

output_dir = "...\\combination_sv\\step3_combi_union\\HG002.annotated.refine.replaced.filteredrange\\" + inputdir + "\\" + outputdir + "\\"

if not os.path.exists(output_dir):
    os.mkdir(output_dir)



vcfs = []
for vcf in os.listdir(vcfdir):
    if vcf.lower()[-4:].find(".vcf") != -1:
        sample = vcf.split(".")[0]
        if(vcf==sample+ "." + mergefile + ".vcf"):
            continue
        vcfs.append(vcf)


for vcf in range(len(vcfs)):
    merge_vcf(vcfdir + vcfs[vcf], output_dir +sample + "." + mergefile + ".vcf")

with open(output_dir +sample + "." + mergefile + ".vcf","r") as merge_vcf:
    svs=[]
    for line in merge_vcf:
        if(line[0]!="#" and line[0:2]!="##"):
            fields = line.split("\t")
            if(len(fields)<9):  # if the count of field is less than 9 fields
                fields.append("")  ## add fields[8]
                fields.append("")  ## add fields[9]

                
            caller = fields[7].split(";")[-1].split("=")[1].strip()
            
            svtype = "unknown"
            for info in fields[7].split(";"):
                if(info.find("SVTYPE")!=-1):
                    svtype = info.split('=')[1]
                    # print(svtype)
            
            for i in range(len(fields)):
                fields[i]=fields[i].strip()
                if len(fields[i])==0:
                    fields[i]="."                  
            
            sv = SV(fields[0],fields[1],fields[2],fields[3],fields[4],\
                fields[5],fields[6],fields[7],fields[8],fields[9],caller,svtype)
            svs.append(sv)
            
    svs.sort(key=lambda x: (len(x.chr), x.chr, x.pos))

    with open(output_dir +sample + "." + firstfile + ".vcf","w") as significantsv_first_f:
        with open(output_dir +sample + "." + allfile + ".vcf","w") as significantsv_all_f:
            write_header(significantsv_first_f, sample)
            write_header(significantsv_all_f, sample)
        
            for m in svs:
                if(not m.used):
                    m.used=True
                    m.candidate_sv(svs)

            for m in svs:
                if(len(m.neighbor_sv)>0) and (m.callerNumber()>0):  # caller > 1 
                    for mm in m.neighbor_sv:
                        significantsv_first_f.write(str(mm).strip() + "\n") 
                        break
            
            for m in svs:
                if(len(m.neighbor_sv)>0) and (m.callerNumber()>0):  # caller > 1 
                    for mm in m.neighbor_sv:
                        significantsv_all_f.write(str(mm).strip() + "\n") 


    ## for statistics                    
    with open(output_dir +sample + "." + statfile + ".tsv","w") as sv_stats:  


        #neighbor_sv (pos +\\-500)       
        neighbor_sv= neighborsv_stats(-1,svs)
        neighbor_sv_0= neighborsv_stats(0,svs)
        neighbor_sv_1= neighborsv_stats(1,svs)
        neighbor_sv_2= neighborsv_stats(2,svs)
        neighbor_sv_3= neighborsv_stats(3,svs)
        neighbor_sv_4= neighborsv_stats(4,svs)
        neighbor_sv_5= neighborsv_stats(5,svs)
        neighbor_sv_ge6= neighborsv_stats(6,svs)

        ########
        neighborsv0call0_count= neighborsv_call_stats(0,0,svs)
        neighborsv1call1_count= neighborsv_call_stats(1,1,svs)
        neighborsv2call1_count= neighborsv_call_stats(2,1,svs)
        neighborsv2call2_count= neighborsv_call_stats(2,2,svs)
        neighborsv3call1_count= neighborsv_call_stats(3,1,svs)
        neighborsv3call2_count= neighborsv_call_stats(3,2,svs)
        neighborsv3call3_count= neighborsv_call_stats(3,3,svs)
        neighborsv4call1_count= neighborsv_call_stats(4,1,svs)
        neighborsv4call2_count= neighborsv_call_stats(4,2,svs)
        neighborsv4call3_count= neighborsv_call_stats(4,3,svs)
        neighborsv4call4_count= neighborsv_call_stats(4,4,svs)
        neighborsv5call1_count= neighborsv_call_stats(5,1,svs)
        neighborsv5call2_count= neighborsv_call_stats(5,2,svs)
        neighborsv5call3_count= neighborsv_call_stats(5,3,svs)
        neighborsv5call4_count= neighborsv_call_stats(5,4,svs)
        neighborsv5call5_count= neighborsv_call_stats(5,5,svs)
        neighborsvge6call1_count= neighborsv_call_stats(6,1,svs)
        neighborsvge6call2_count= neighborsv_call_stats(6,2,svs)
        neighborsvge6call3_count= neighborsv_call_stats(6,3,svs)
        neighborsvge6call4_count= neighborsv_call_stats(6,4,svs)
        neighborsvge6call5_count= neighborsv_call_stats(6,5,svs)
        neighborsvge6callge6_count= neighborsv_call_stats(6,6,svs)
        significant_sv= neighborsv_call_stats(-1,-1,svs)

        sv_stats.write("total sv\t"+ str(len(svs)) +"\n")
        sv_stats.write("total neighbor_sv\t"+ str(neighbor_sv) + "\n")
        sv_stats.write("neighbor_sv_eq0\t"+ str(neighbor_sv_0) + "\n")
        sv_stats.write("neighbor_sv_eq1\t"+ str(neighbor_sv_1) + "\n")
        sv_stats.write("neighbor_sv_eq2\t"+ str(neighbor_sv_2) + "\n")
        sv_stats.write("neighbor_sv_eq3\t"+ str(neighbor_sv_3) + "\n")
        sv_stats.write("neighbor_sv_eq4\t"+ str(neighbor_sv_4) + "\n")
        sv_stats.write("neighbor_sv_eq5\t"+ str(neighbor_sv_5) + "\n")
        sv_stats.write("neighbor_sv_ge6\t"+ str(neighbor_sv_ge6) + "\n")
        sv_stats.write("neighborsv0_0call\t"+ str(neighborsv0call0_count) +"\n")
        sv_stats.write("neighborsv1_1call\t"+ str(neighborsv1call1_count) +"\n")
        sv_stats.write("neighborsv2_1call\t"+ str(neighborsv2call1_count) +"\n")
        sv_stats.write("neighborsv2_2call\t"+ str(neighborsv2call2_count) +"\n")
        sv_stats.write("neighborsv3_1call\t"+ str(neighborsv3call1_count) +"\n")
        sv_stats.write("neighborsv3_2call\t"+ str(neighborsv3call2_count) +"\n")
        sv_stats.write("neighborsv3_3call\t"+ str(neighborsv3call3_count) +"\n")
        sv_stats.write("neighborsv4_1call\t"+ str(neighborsv4call1_count) +"\n")
        sv_stats.write("neighborsv4_2call\t"+ str(neighborsv4call2_count) +"\n")
        sv_stats.write("neighborsv4_3call\t"+ str(neighborsv4call3_count) +"\n")
        sv_stats.write("neighborsv4_4call\t"+ str(neighborsv4call4_count) +"\n")
        sv_stats.write("neighborsv5_1call\t"+ str(neighborsv5call1_count) +"\n")
        sv_stats.write("neighborsv5_2call\t"+ str(neighborsv5call2_count) +"\n")
        sv_stats.write("neighborsv5_3call\t"+ str(neighborsv5call3_count) +"\n")
        sv_stats.write("neighborsv5_4call\t"+ str(neighborsv5call4_count) +"\n")
        sv_stats.write("neighborsv5_5call\t"+ str(neighborsv5call5_count) +"\n")
        sv_stats.write("neighborsvge6_1call\t"+ str(neighborsvge6call1_count) +"\n")
        sv_stats.write("neighborsvge6_2call\t"+ str(neighborsvge6call2_count) +"\n")
        sv_stats.write("neighborsvge6_3call\t"+ str(neighborsvge6call3_count) +"\n")
        sv_stats.write("neighborsvge6_4call\t"+ str(neighborsvge6call4_count) +"\n")
        sv_stats.write("neighborsvge6_5call\t"+ str(neighborsvge6call5_count) +"\n")
        sv_stats.write("neighborsvge6_ge6call\t"+ str(neighborsvge6callge6_count) +"\n")
        sv_stats.write("total significant_sv\t"+ str(significant_sv))
