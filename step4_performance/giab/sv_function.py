class SV():
    def __init__(self,chr,pos,id,ref,alt,qual,filter,info,format,others,caller,svtype):

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
        self.caller = caller
        self.svtype = svtype
        self.used = False
        self.neighbor_sv = []
        self.neighbor_sv_caller = []
    
    def __str__(self):
        return self.chr + "\t" + str(self.pos) + "\t" +self.id + "\t" +self.ref + "\t" +self.alt + "\t" +self.qual + "\t" + \
               self.filter + "\t" + self.info + "\t" + self.format + "\t" + self.others
               
    def candidate_sv(self, svs):
        self.neighbor_sv.append(self)
        self.neighbor_sv_caller.append(self.caller)
        for i in svs:
            if(not i.used):
                if i.chr==self.chr and (i.pos<=self.pos+500 and i.pos>=self.pos-500 and self.svtype!="unknown" and self.svtype==i.svtype):
                    i.used=True
                    self.neighbor_sv.append(i)
                    self.neighbor_sv_caller.append(i.caller)
    
    def callerNumber(self):
        return len(set(self.neighbor_sv_caller))

class write_vcf():
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

class Variant():
    def __init__(self,chr,pos,svtype):
        self.chr = chr
        self.pos = int(pos)
        self.svtype = svtype
        
    def __str__(self):
        return self.chr + "\t" + str(self.pos) + "\t" + self.svtype


def findInfo(infovalue, something):
    isFind=False
    for s in infovalue:
        if s.find(something)!=-1:
            result = s.split("=")[1]
            isFind=True
            break
    if(isFind):
        return result
    else:
        return "."


def findFormat(formatfield, field9):
    if formatfield.split(":")[0] == "GT":
        return field9.split(":")[0]
    else:
        return "."



def write_header(significantsv_f, sample):
    significantsv_f.write("##fileformat=VCFv4.2\n")
    significantsv_f.write("##reference=ucsc.hg19.fasta\n")
    significantsv_f.write("##contig=<ID=chrM,length=16571>\n")
    significantsv_f.write("##contig=<ID=chr1,length=249250621>\n")
    significantsv_f.write("##contig=<ID=chr2,length=243199373>\n")
    significantsv_f.write("##contig=<ID=chr3,length=198022430>\n")
    significantsv_f.write("##contig=<ID=chr4,length=191154276>\n")
    significantsv_f.write("##contig=<ID=chr5,length=180915260>\n")
    significantsv_f.write("##contig=<ID=chr6,length=171115067>\n")
    significantsv_f.write("##contig=<ID=chr7,length=159138663>\n")
    significantsv_f.write("##contig=<ID=chr8,length=146364022>\n")
    significantsv_f.write("##contig=<ID=chr9,length=141213431>\n")
    significantsv_f.write("##contig=<ID=chr10,length=135534747>\n")
    significantsv_f.write("##contig=<ID=chr11,length=135006516>\n")
    significantsv_f.write("##contig=<ID=chr12,length=133851895>\n")
    significantsv_f.write("##contig=<ID=chr13,length=115169878>\n")
    significantsv_f.write("##contig=<ID=chr14,length=107349540>\n")
    significantsv_f.write("##contig=<ID=chr15,length=102531392>\n")
    significantsv_f.write("##contig=<ID=chr16,length=90354753>\n")
    significantsv_f.write("##contig=<ID=chr17,length=81195210>\n")
    significantsv_f.write("##contig=<ID=chr18,length=78077248>\n")
    significantsv_f.write("##contig=<ID=chr19,length=59128983>\n")
    significantsv_f.write("##contig=<ID=chr20,length=63025520>\n")
    significantsv_f.write("##contig=<ID=chr21,length=48129895>\n")
    significantsv_f.write("##contig=<ID=chr22,length=51304566>\n")
    significantsv_f.write("##contig=<ID=chrX,length=155270560>\n")
    significantsv_f.write("##contig=<ID=chrY,length=59373566>\n")
    significantsv_f.write("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">\n")
    significantsv_f.write("##INFO=<ID=SVLEN,Number=.,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">\n")
    significantsv_f.write("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">\n")
    significantsv_f.write("##INFO=<ID=STRANDS,Number=.,Type=String,Description=\"Strand orientation of the adjacency in BEDPE format (DEL:+-, DUP:-+, INV:++/--)\">\n")
    significantsv_f.write("##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise structural variation\">\n")
    significantsv_f.write("##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"Confidence interval around POS for imprecise variants\">\n")
    significantsv_f.write("##INFO=<ID=CIEND,Number=2,Type=Integer,Description=\"Confidence interval around END for imprecise variants\">\n")
    significantsv_f.write("##INFO=<ID=CIPOS95,Number=2,Type=Integer,Description=\"Confidence interval (95%) around POS for imprecise variants\">\n")
    significantsv_f.write("##INFO=<ID=CIEND95,Number=2,Type=Integer,Description=\"Confidence interval (95%) around END for imprecise variants\">\n")
    significantsv_f.write("##INFO=<ID=MATEID,Number=.,Type=String,Description=\"ID of mate breakends\">\n")
    significantsv_f.write("##INFO=<ID=EVENT,Number=1,Type=String,Description=\"ID of event associated to breakend\">\n")
    significantsv_f.write("##INFO=<ID=SECONDARY,Number=0,Type=Flag,Description=\"Secondary breakend in a multi-line variants\">\n")
    significantsv_f.write("##INFO=<ID=SU,Number=.,Type=Integer,Description=\"Number of pieces of evidence supporting the variant across all samples\">\n")
    significantsv_f.write("##INFO=<ID=PE,Number=.,Type=Integer,Description=\"Number of paired-end reads supporting the variant across all samples\">\n")
    significantsv_f.write("##INFO=<ID=SR,Number=.,Type=Integer,Description=\"Number of split reads supporting the variant across all samples\">\n")
    significantsv_f.write("##INFO=<ID=BD,Number=.,Type=Integer,Description=\"Amount of BED evidence supporting the variant across all samples\">\n")
    significantsv_f.write("##INFO=<ID=EV,Number=.,Type=String,Description=\"Type of LUMPY evidence contributing to the variant call\">\n")
    significantsv_f.write("##INFO=<ID=PRPOS,Number=.,Type=String,Description=\"LUMPY probability curve of the POS breakend\">\n")
    significantsv_f.write("##INFO=<ID=PREND,Number=.,Type=String,Description=\"LUMPY probability curve of the END breakend\">\n")
    significantsv_f.write("##INFO=<ID=CIGAR,Number=A,Type=String,Description=\"CIGAR alignment for each alternate indel allele\">\n")
    significantsv_f.write("##INFO=<ID=HOMLEN,Number=.,Type=Integer,Description=\"Length of base pair identical homology at event breakpoints\">\n")
    significantsv_f.write("##INFO=<ID=HOMSEQ,Number=.,Type=String,Description=\"Sequence of base pair identical homology at event breakpoints\">\n")
    significantsv_f.write("##INFO=<ID=SVINSLEN,Number=.,Type=Integer,Description=\"Length of insertion\">\n")
    significantsv_f.write("##INFO=<ID=SVINSSEQ,Number=.,Type=String,Description=\"Sequence of insertion\">\n")
    significantsv_f.write("##INFO=<ID=LEFT_SVINSSEQ,Number=.,Type=String,Description=\"Known left side of insertion for an insertion of unknown length\">\n")
    significantsv_f.write("##INFO=<ID=RIGHT_SVINSSEQ,Number=.,Type=String,Description=\"Known right side of insertion for an insertion of unknown length\">\n")
    significantsv_f.write("##INFO=<ID=BND_DEPTH,Number=1,Type=Integer,Description=\"Read depth at local translocation breakend\">\n")
    significantsv_f.write("##INFO=<ID=MATE_BND_DEPTH,Number=1,Type=Integer,Description=\"Read depth at remote translocation mate breakend\">\n")
    significantsv_f.write("##INFO=<ID=JUNCTION_QUAL,Number=1,Type=Integer,Description=\"If the SV junction is part of an EVENT (ie. a multi-adjacency variant), this field provides the QUAL value for the adjacency in question only\">\n")
    significantsv_f.write("##INFO=<ID=SPAN,Number=1,Type=Integer,Description=\"Distance between the breakpoints. -1 for interchromosomal\">\n")
    significantsv_f.write("##INFO=<ID=INSERTION,Number=1,Type=String,Description=\"Sequence insertion at the breakpoint.\">\n")
    significantsv_f.write("##INFO=<ID=SCTG,Number=1,Type=String,Description=\"Identifier for the contig assembled by svaba to make the SV call\">\n")
    significantsv_f.write("##INFO=<ID=EVDNC,Number=1,Type=String,Description=\"Evidence for variant. ASSMB assembly only, ASDIS assembly+discordant. DSCRD discordant only, TSI_L templated-sequence insertion (local, e.g. AB or BC of an ABC), TSI_G global (e.g. AC of ABC)\">\n")
    significantsv_f.write("##INFO=<ID=BX,Number=.,Type=String,Description=\"Table of BX tag counts for supporting reads\">\n")
    significantsv_f.write("##INFO=<ID=NM,Number=1,Type=Integer,Description=\"Number of mismatches of this alignment fragment to reference\">\n")
    significantsv_f.write("##INFO=<ID=SUBN,Number=1,Type=Integer,Description=\"Number of secondary alignments associated with this contig fragment\">\n")
    significantsv_f.write("##INFO=<ID=DISC_MAPQ,Number=1,Type=Integer,Description=\"Mean mapping quality of discordant reads mapped here\">\n")
    significantsv_f.write("##INFO=<ID=REPSEQ,Number=1,Type=String,Description=\"Repeat sequence near the event\">\n")
    significantsv_f.write("##INFO=<ID=READNAMES,Number=.,Type=String,Description=\"IDs of ALT reads\">\n")
    significantsv_f.write("##INFO=<ID=MAPQ,Number=1,Type=Integer,Description=\"Mapping quality (BWA-MEM) of this fragement of the contig (-1 if discordant only)\">\n")
    significantsv_f.write("##INFO=<ID=MATEMAPQ,Number=1,Type=Integer,Description=\"Mapping quality of the partner fragment of the contig\">\n")
    significantsv_f.write("##INFO=<ID=MATENM,Number=1,Type=Integer,Description=\"Number of mismatches of partner alignment fragment to reference\">\n")
    significantsv_f.write("##INFO=<ID=NUMPARTS,Number=1,Type=Integer,Description=\"If detected with assembly, number of parts the contig maps to. Otherwise 0\">\n")
    significantsv_f.write("##INFO=<ID=CHR2,Number=1,Type=String,Description=\"Chromosome for POS2 coordinate in case of an inter-chromosomal translocation\">\n")
    significantsv_f.write("##INFO=<ID=POS2,Number=1,Type=Integer,Description=\"Genomic position for CHR2 in case of an inter-chromosomal translocation\">\n")
    significantsv_f.write("##INFO=<ID=SRMAPQ,Number=1,Type=Integer,Description=\"Median mapping quality of split-reads\">\n")
    significantsv_f.write("##INFO=<ID=SRQ,Number=1,Type=Float,Description=\"Split-read consensus alignment quality\">\n")
    significantsv_f.write("##INFO=<ID=CONSENSUS,Number=1,Type=String,Description=\"Split-read consensus sequence\">\n")
    significantsv_f.write("##INFO=<ID=CE,Number=1,Type=Float,Description=\"Consensus sequence entropy\">\n")
    significantsv_f.write("##INFO=<ID=CT,Number=1,Type=String,Description=\"Paired-end signature induced connection type\">\n")
    significantsv_f.write("##INFO=<ID=PRECISE,Number=0,Type=Flag,Description=\"Precise structural variation\">\n")
    significantsv_f.write("##INFO=<ID=SVMETHOD,Number=1,Type=String,Description=\"Type of approach used to detect SV\">\n")
    significantsv_f.write("##ALT=<ID=DEL,Description=\"Deletion\">\n")
    significantsv_f.write("##ALT=<ID=DUP,Description=\"Duplication\">\n")
    significantsv_f.write("##ALT=<ID=INV,Description=\"Inversion\">\n")
    significantsv_f.write("##ALT=<ID=DUP:TANDEM,Description=\"Tandem duplication\">\n")
    significantsv_f.write("##ALT=<ID=INS,Description=\"Insertion of novel sequence\">\n")
    significantsv_f.write("##ALT=<ID=CNV,Description=\"Copy number variable region\">\n")
    significantsv_f.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
    significantsv_f.write("##FORMAT=<ID=SU,Number=1,Type=Integer,Description=\"Number of pieces of evidence supporting the variant\">\n")
    significantsv_f.write("##FORMAT=<ID=PE,Number=1,Type=Integer,Description=\"Number of paired-end reads supporting the variant\">\n")
    significantsv_f.write("##FORMAT=<ID=SR,Number=1,Type=Integer,Description=\"Number of split reads supporting the variant\">\n")
    significantsv_f.write("##FORMAT=<ID=BD,Number=1,Type=Integer,Description=\"Amount of BED evidence supporting the variant\">\n")
    significantsv_f.write("##FORMAT=<ID=FT,Number=1,Type=String,Description=\"Sample filter, 'PASS' indicates that all filters have passed for this sample\">\n")
    significantsv_f.write("##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">\n")
    significantsv_f.write("##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification\">\n")
    significantsv_f.write("##FORMAT=<ID=PR,Number=.,Type=Integer,Description=\"Spanning paired-read support for the ref and alt alleles in the order listed\">\n")  
    significantsv_f.write("##FORMAT=<ID=SR,Number=.,Type=Integer,Description=\"Split reads for the ref and alt alleles in the order listed, for reads where P(allele|read)>0.999\">\n")
    significantsv_f.write("##FORMAT=<ID=SL,Number=1,Type=Float,Description=\"Alignment-quality Scaled log-odds, where LO is LO * (MAPQ - 2*NM)/60\">\n")
    significantsv_f.write("##FORMAT=<ID=LO,Number=1,Type=Float,Description=\"Log-odds that this variant is real vs artifact\">\n")
    significantsv_f.write("##FORMAT=<ID=DR,Number=1,Type=Integer,Description=\"Number of discordant-supported reads for this variant\">\n")
    significantsv_f.write("##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth of coverage: Number of reads covering site.\">\n")
    significantsv_f.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Most likely genotype\">\n")
    significantsv_f.write("##FORMAT=<ID=AD,Number=1,Type=Integer,Description=\"Allele depth: Number of reads supporting the variant\">\n")
    significantsv_f.write("##FORMAT=<ID=GQ,Number=1,Type=String,Description=\"Genotype quality (currently not supported. Always 0)\">\n")
    significantsv_f.write("##FORMAT=<ID=LR,Number=1,Type=Float,Description=\"Log-odds that this variant is REF vs AF=0.5\">\n")
    significantsv_f.write("##FORMAT=<ID=PL,Number=.,Type=Float,Description=\"Normalized likelihood of the current genotype\">\n")
    significantsv_f.write("##FORMAT=<ID=SR,Number=1,Type=Integer,Description=\"Number of spanning reads for this variant\">\n")
    significantsv_f.write("##FORMAT=<ID=GL,Number=G,Type=Float,Description=\"Log10-scaled genotype likelihoods for RR,RA,AA genotypes\">\n")
    significantsv_f.write("##FORMAT=<ID=RC,Number=1,Type=Integer,Description=\"Raw high-quality read counts or base counts for the SV\">\n")
    significantsv_f.write("##FORMAT=<ID=RCL,Number=1,Type=Integer,Description=\"Raw high-quality read counts or base counts for the left control region\">\n")
    significantsv_f.write("##FORMAT=<ID=RCR,Number=1,Type=Integer,Description=\"Raw high-quality read counts or base counts for the right control region\">\n")
    significantsv_f.write("##FORMAT=<ID=RDCN,Number=1,Type=Integer,Description=\"Read-depth based copy-number estimate for autosomal sites\">\n")
    significantsv_f.write("##FORMAT=<ID=DR,Number=1,Type=Integer,Description=\"# high-quality reference pairs\">\n")
    significantsv_f.write("##FORMAT=<ID=DV,Number=1,Type=Integer,Description=\"# high-quality variant pairs\">\n")
    significantsv_f.write("##FORMAT=<ID=RR,Number=1,Type=Integer,Description=\"# high-quality reference junction reads\">\n")
    significantsv_f.write("##FORMAT=<ID=RV,Number=1,Type=Integer,Description=\"# high-quality variant junction reads\">\n")
    significantsv_f.write("##FILTER=<ID=Ploidy,Description=\"For DEL & DUP variants, the genotypes of overlapping variants (with similar size) are inconsistent with diploid expectation\">\n")
    significantsv_f.write("##FILTER=<ID=MaxDepth,Description=\"Depth is greater than 3x the median chromosome depth near one or both variant breakends\">\n")
    significantsv_f.write("##FILTER=<ID=MaxMQ0Frac,Description=\"For a small variant (<1000 bases), the fraction of reads in all samples with MAPQ0 around either breakend exceeds 0.4\">\n")   
    significantsv_f.write("##FILTER=<ID=NoPairSupport,Description=\"For variants significantly larger than the paired read fragment size, no paired reads support the alternate allele in any sample.\">\n")
    significantsv_f.write("##FILTER=<ID=MinQUAL,Description=\"QUAL score is less than 20\">\n")
    significantsv_f.write("##FILTER=<ID=SampleFT,Description=\"No sample passes all the sample-level filters (at the field FORMAT/FT)\">\n")
    significantsv_f.write("##FILTER=<ID=MinGQ,Description=\"GQ score is less than 15 (filter applied at sample level)\">\n")
    significantsv_f.write("##FILTER=<ID=HomRef,Description=\"homozygous reference call (filter applied at sample level)\">\n")
    significantsv_f.write("##FILTER=<ID=LOWSPANDSCRD,Description=\"Discordant-only cluster is too small given isize distribution to call confidently\">\n")
    significantsv_f.write("##FILTER=<ID=MULTIMATCH,Description=\"Low MAPQ and this contig fragment maps well to multiple locations\">\n")
    significantsv_f.write("##FILTER=<ID=PASS,Description=\"Strong assembly support, strong discordant support, or combined support. Strong MAPQ\">\n")
    significantsv_f.write("##FILTER=<ID=TOOSHORT,Description=\"Contig alignment for part of this rearrangement has <= 25bp match to reference\">\n")
    significantsv_f.write("##FILTER=<ID=WEAKSUPPORTHIREP,Description=\"Fewer then 7 split reads for variant with >= 10 bases of repeat sequence (need to be more strict)\">\n")
    significantsv_f.write("##FILTER=<ID=LOCALMATCH,Description=\"Contig realigned to assembly region without clipping\">\n")
    significantsv_f.write("##FILTER=<ID=DUPREADS,Description=\"Contig built from what appear to be duplicate reads (split reads all same contig cov))\">\n")
    significantsv_f.write("##FILTER=<ID=LOWSPLITSMALL,Description=\"Fewer than 4 split reads for small events ( < 1500 bp)\">\n")
    significantsv_f.write("##FILTER=<ID=SIMPLESEQUENCE,Description=\"Major portion of one contig mapping falls in a simple sequence, as given by -R flag. Assembly-only filter\">\n")
    significantsv_f.write("##FILTER=<ID=HIGHHOMOLOGY,Description=\"Contig realigns with > 25% of readlength of homology. High probaility of assembly/mapping artifact\">\n")
    significantsv_f.write("##FILTER=<ID=WEAKDISC,Description=\"Fewer than 7 supporting discordant reads and no assembly support\">\n")
    significantsv_f.write("##FILTER=<ID=COMPETEDISC,Description=\"Discordant cluster found with nearly same breakpoints, but different strands for DSCRD event\">\n")
    significantsv_f.write("##FILTER=<ID=LOWSUPPORT,Description=\"Fewer than 2 split reads or < 4 total alt reads for ASDISC\">\n")
    significantsv_f.write("##FILTER=<ID=NOLOCAL,Description=\"Contig realigned to region outside of local assembly region, and no disc support.\">\n")
    significantsv_f.write("##FILTER=<ID=LOWSPAN,Description=\"Discordant read cluster (no split read support), and less than 10kb span and < 12 reads\">\n")
    significantsv_f.write("##FILTER=<ID=LOWAS,Description=\"Alignment score of one end is less than 80% of contig length, or number of mismatch bases (NM) on one end is >= 10\">\n")
    significantsv_f.write("##FILTER=<ID=NODISC,Description=\"Rearrangement was not detected independently by assembly\">\n")
    significantsv_f.write("##FILTER=<ID=LOWQINVERSION,Description=\"Assembly-only inversion of span < 300 and < 6 split reads. Common artifact in Illumina data\">\n")
    significantsv_f.write("##FILTER=<ID=LOWMATCHLEN,Description=\"Assembly contig has fewer than 40 bases mapping uniquely to a reference locus (<100 if complex mapping or \">\n")
    significantsv_f.write("##FILTER=<ID=SINGLEBX,Description=\"Variant is supported by only a single BX tag (if run with 10X Genomics data)\">\n")
    significantsv_f.write("##FILTER=<ID=LOWMAPQDISC,Description=\"Both clusters of reads failed to achieve mean mapq of > 30 for DSCRD\">\n")
    significantsv_f.write("##FILTER=<ID=LOWMAPQ,Description=\"Assembly contig has non 60/60 mapq and no discordant support\">\n")
    significantsv_f.write("##FILTER=<ID=LOWICSUPPORT,Description=\"Less than 60bp of contig match on one end of an inter-chromosomal break\">\n")
    significantsv_f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT" + "\t" + sample + "\n")

def neighborsv_stats(neighbor_svNumber, svs):
    neighborsv_count=0
    for m in svs:
        if neighbor_svNumber==-1:
            if (len(m.neighbor_sv)>1):
                neighborsv_count +=1
        if neighbor_svNumber>=6:
            if (len(m.neighbor_sv)>=neighbor_svNumber):
                neighborsv_count +=1
        if neighbor_svNumber<6:
            if (len(m.neighbor_sv)==neighbor_svNumber):
                neighborsv_count +=1
    return neighborsv_count

def neighborsv_call_stats(neighbor_svNumber, callerNumber, svs):
    neighborsv_call_count=0
    for m in svs:
        if neighbor_svNumber>=6 and callerNumber<6:
            if (len(m.neighbor_sv)>=neighbor_svNumber and m.callerNumber()==callerNumber):
                neighborsv_call_count +=1
        if neighbor_svNumber>=6 and callerNumber>=6:
            if (len(m.neighbor_sv)>=neighbor_svNumber and m.callerNumber()>=callerNumber):
                neighborsv_call_count +=1
        if neighbor_svNumber<6 and callerNumber>=6:
            if (len(m.neighbor_sv)==neighbor_svNumber and m.callerNumber()>=callerNumber):
                neighborsv_call_count +=1
        if neighbor_svNumber<5 and callerNumber<6:
            if (len(m.neighbor_sv)==neighbor_svNumber and m.callerNumber()==callerNumber):
                neighborsv_call_count +=1
        if neighbor_svNumber==-1 and callerNumber==-1:
            if (len(m.neighbor_sv)>1 and m.callerNumber()>1):
                neighborsv_call_count +=1
    return neighborsv_call_count
