import os
from sv_function import *
import sys
from operator import attrgetter


# combi = "intersect"  # choose this if intersect
combi = "union" # choose this if union

# model = "mdg"  # chooses this if combination of 3 callers
model = "mdgls"  # choose this if combination of 5 callers


folder = combi + ".HG002." + model
outfile = combi + ".HG002." + model + ".intv"
infile = "HG002.mergedsv." + model

# keywords = ["manta", "delly", "gridss"]
# keywords = ["manta", "delly", "gridss", "lumpy", "svaba"]


output_dir = "...\\combination_sv\\step4_performance\\models\\for_intv\\"
if not os.path.exists(output_dir):
    os.mkdir(output_dir)

vcfdir = "...\\combination_sv\\step4_performance\\models\\" + folder + "\\"

# print(output_dir)
# print(vcfdir)
# print(outfile)
# print(vcfdir + infile + ".vcf")
# print(output_dir + outfile + ".tsv")


with open(vcfdir + infile + ".vcf","r") as merge_vcf:
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

    with open(output_dir + outfile + ".tsv","w") as intv_f:    
        intv_f.write( "ID" + "\t" + "INTERVAL" + "\t" + "COUNT" + "\t" + "DETAILS" + "\t" + "CALLERS" + "\n")    
        for m in svs:
            if(not m.used):
                m.used=True
                m.candidate_sv(svs)


        for m in svs:
            # if(len(m.neighbor_sv)>1) and (m.callerNumber()>1):  # choose this if intersect
            if(len(m.neighbor_sv)>0) and (m.callerNumber()>0):  # choose this if union
                max_attr = max(m.neighbor_sv, key=attrgetter('pos'))
                min_attr = min(m.neighbor_sv, key=attrgetter('pos'))
                intv = max_attr.pos - min_attr.pos
                s = ""
                c = ""
                for mm in m.neighbor_sv:
                    s += str(mm.pos) + ","
                    c += str(mm.caller) + ","
                s = s[0:len(s)-1]
                c = c[0:len(c)-1]
                intv_f.write(str(m.chr) + "-" + str(m.pos) + "\t" + str(intv) + "\t" + str(len(m.neighbor_sv)) \
                    + "\t" + s + "\t" + c + "\n")