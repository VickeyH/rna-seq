#!/usr/bin/env python
#coding:utf-8
import argparse,os,sys,datetime,commands,glob,re

def mkdir(dir,sampls):
    for n in sampls:
        d = os.path.join(dir,n)
        if not os.path.exists(d):
            os.makedirs(d)
        else:
            continue
    
def listdir(path):
    rel_dir = []
    for a,b,c in os.walk(path):
        for i in c:
            rel_dir.append(os.path.join(a,i))
    return rel_dir
            
def check_makeflow(makeflow):
    input_file,out_file,result_file,need_file,no_file = set(),set(),set(),set(),set()
    with open(makeflow,"r") as f:
        for line in f:
            line = line.strip()
            if ":" in line and line.startswith("/"):
                out = line.split(":")
                for i in out[0].split():
                    if "results" in i:
                        result_file.add(i.replace("\\",""))
                    else:
                        out_file.add(i.replace("\\",""))
                for i in out[1].split():
                    if "results" in i:
                        result_file.add(i.replace("\\",""))
                    else:
                        input_file.add(i.replace("\\",""))
            else:
                continue
    need_file = input_file.difference(out_file,result_file)
    for i in need_file:
        if os.path.exists(i):
            continue
	else:
		no_file.add(i)
    res_str = ""    
    if len(no_file):
        for i in sorted(list(no_file)):
            res_str += i + "\t Not Exists!\n"
        return False,res_str
    else:
        res_str += "All need files exists in %s, OK!\n"%makeflow
        return True,res_str  
        
def map_sample_and_raw(sample_list,path):
    raw_d = listdir(path)
    raw_n = map(lambda x:re.sub("\W+","_",x),map(os.path.basename,raw_d))
    tmp_d  = dict(zip(raw_n,raw_d))
    sample_raw={}
    for sn in sample_list:
        for p in raw_n:
            if "_%s_"%sn in p:
                sample_raw.setdefault(sn,[]).append(tmp_d[p])
        if not sample_raw[sn]:
            print "no %s in rawdata directory,please check" %sn
            sys.exit(0)
    for k in sample_raw:
        kv =  []
        for n in sample_raw[k]:
            bn = re.sub("\W",lambda m:"\\"+m.group(0),os.path.basename(n))
            bn = bn.replace("\\.",".")
            kv.append(os.path.join(os.path.dirname(n),bn))
        kv.sort()
        sample_raw[k] = kv
    return sample_raw

def write_cat_makeflow(sample_dict,out_dir,makeflow):
    if os.path.exists(makeflow):
        print "%s file exists,please check" %makeflow
        sys.exit(0)
    fo =open(makeflow,"a+")
    for s in sample_dict:
        for r in "12":
            fo.write("CATEGORY=%s\n"%s)
            fo.write(os.path.abspath(os.path.join(out_dir,s)) + '/R%s.fastq.gz : %s\n' %(r," ".join(sample_dict[s][int(r)-1::2])))
            fo.write("\tcat %s >> %s\n\n" %(" ".join(sample_dict[s][int(r)-1::2]),os.path.abspath(os.path.join(out_dir,s)) + '/R%s.fastq.gz'%r))           
    fo.close()       


def write_trimmomatic_qc_makeflow(samples,makeflow,input_d,output_d):
    fo =open(makeflow,"a+")
    for sn in samples:
        r1 = os.path.join(input_d,sn,"R1.fastq.gz")
        r2 = os.path.join(input_d,sn,"R2.fastq.gz")
        od = os.path.join(output_d,sn)
        logfile = os.path.join(od,"%s.trim.log"%sn)
        fo.write("CATEGORY=QC_%s\n"%sn)
        fo.write(os.path.join(od,"R1.clean.fq.gz ") + os.path.join(od,"R1.clean.unpaired.fq.gz ") + os.path.join(od,"R2.clean.fq.gz ") + os.path.join(od,"R2.clean.unpaired.fq.gz : ") + r1 + " " + r2 + "\n")
        fo.write("\tjava -jar /lustre/work/yongdeng/software/protokaryon/soft/QC/Trimmomatic/trimmomatic-0.36.jar PE -threads 8 -phred33 -trimlog %s %s %s %s %s %s %s ILLUMINACLIP:/lustre/work/yongdeng/software/protokaryon/soft/QC/Trimmomatic/adapters/adapter:2:30:10 LEADING:15 TRAILING:15 SLIDINGWINDOW:4:15 MINLEN:50 " %(logfile,r1,r2,os.path.join(od,"R1.clean.fq.gz"),os.path.join(od,"R1.clean.unpaired.fq.gz"),os.path.join(od,"R2.clean.fq.gz"),os.path.join(od,"R2.clean.unpaired.fq.gz")) + "> " + os.path.join(od,'trim.out') + " 2> " + os.path.join(od,'trim.err') + "\n\n")       
    fo.close()
    
def write_IlluQC_qc_makeflow(samples,makeflow,input_d,output_d):
    fo =open(makeflow,"a+")
    for sn in samples:
        r1 = os.path.join(input_d,sn,"R1.fastq.gz")
        r2 = os.path.join(input_d,sn,"R2.fastq.gz")
        od = os.path.join(output_d,sn)
        fo.write("CATEGORY=QC_%s\n"%sn)
        fo.write(" ".join([os.path.join(od,i) for i in ["R1.clean.fq.gz","R2.clean.fq.gz","summary.png","Filter_stat","R1_avgQual.png","R2_avgQual.png"]]) + " : " + r1 + " " + r2 + "\n")
        fo.write("\tperl /lustre/work/yongdeng/software/protokaryon/soft/QC/IlluQC_PRLL_CapitalBio_Adapter.pl -pe %s %s /lustre/work/yongdeng/software/protokaryon/soft/QC/adapter 5 -c 8 -z g -o %s > %s 2> %s\n\n"%(r1,r2,od,os.path.join(od,"IlluQC_N.out"),os.path.join(od,"IlluQC_N.err")))
    fo.close()


def write_align_makeflow(samples,genome_dir,makeflow,input_d,output_d):
    if not os.path.exists(genome_dir):
        print "%s dir not exists, please check!" %genome_dir
        sys.exit(0)
    gtf = glob.glob(os.path.join(genome_dir,"*.gtf"))[0]
    fa = glob.glob(os.path.join(genome_dir,"*.fa"))[0]
    if not fa or not gtf:
        print "please check the fa and gtf file in %s" %genome_dir
        sys.exists(0)
    genome_prefix = os.path.splitext(fa)[0]
    fo =open(makeflow,"a+")
    fo.write("CATEGORY=Bowtie2-build\n")
    fo.write("%s : %s\n"%(genome_dir+"/build-genome.log",fa))
    fo.write("\t/lustre/work/yongdeng/software/protokaryon/soft/Align/bowtie2-2.3.3.1/bowtie2-build --threads 8 %s %s &> %s\n\n" %(fa,genome_prefix,genome_dir+"/build-genome.log"))
    fo.write("CATEGORY=create_refFlat\n")
    fo.write(gtf + ".refFlat : " + gtf + "\n")
    fo.write("\tperl /lustre/work/yongdeng/software/protokaryon/soft/Align/create_refFlat.pl %s %s &> %s\n\n"%(gtf,gtf + ".refFlat",genome_dir+"/refflat.err"))
    for sn in samples:
        r1 = os.path.join(input_d,sn,"R1.clean.fq.gz")
        r2 = os.path.join(input_d,sn,"R2.clean.fq.gz")
        od = os.path.join(output_d,sn)
        fo.write("CATEGORY=Bowtie2-align_%s\n"%sn)
        fo.write(os.path.join(od,"bowtie.err ") + os.path.join(od,"accepted_hits.sam : ") + r1 + " " + r2 + " " + genome_dir + "/build-genome.log\n")
        fo.write("@BATCH_OPTIONS = -l h_vmem=50G\n")
        fo.write("\t/lustre/work/yongdeng/software/protokaryon/soft/Align/bowtie2-2.3.3.1/bowtie2 -p 8 -x %s -1 %s -2 %s -S %s "%(genome_prefix,r1,r2,os.path.join(od,"accepted_hits.sam")))
        fo.write("> %s 2> %s\n\n" %(os.path.join(od,"bowtie.out"),os.path.join(od,"bowtie.err")))
    fo.close()

def write_fastqc_makeflow(raw_dir,clean_dir,makeflow,samples):
    fo =open(makeflow,"a+")
    for sn in samples:
        od_raw = os.path.join(os.path.dirname(os.path.normpath(raw_dir)),"fastqcRaw",sn)
        od_clean = os.path.join(os.path.dirname(os.path.normpath(clean_dir)),"fastqcClean",sn)
        raw_r1 = os.path.join(raw_dir,sn,"R1.fastq.gz")
        raw_r2 = os.path.join(raw_dir,sn,"R2.fastq.gz")
        clean_r1 = os.path.join(clean_dir,sn,"R1.clean.fq.gz")
        clean_r2 = os.path.join(clean_dir,sn,"R2.clean.fq.gz")
        fo.write("CATEGORY=fastqc_Raw_%s\n"%sn)
        fo.write(od_raw + "/R1_fastqc.html " + od_raw + "/R2_fastqc.html : " + raw_r1 + " " + raw_r2 + "\n")
        fo.write("\t/lustre/work/yongdeng/software/protokaryon/soft/QC/FastQC/fastqc -o %s -f fastq -t 8 --extract %s %s " %(od_raw,raw_r1,raw_r2))
        fo.write("> %s 2> %s\n\n"%(od_raw + "/fastqc.out",od_raw + "/fastqc.err"))        
        fo.write("CATEGORY=fastqc_Clean_%s\n"%sn)
        fo.write(od_clean + "/R1.clean_fastqc.html " + od_clean + "/R2.clean_fastqc.html : " + clean_r1 + " " + clean_r2 + "\n")
        fo.write("\t/lustre/work/yongdeng/software/protokaryon/soft/QC/FastQC/fastqc -o %s -f fastq -t 8 --extract %s %s "%(od_clean,clean_r1,clean_r2))
        fo.write("> %s 2> %s\n\n"%(od_clean + "/fastqc.out",od_clean + "/fastqc.err"))
    fo.close()    

def write_summary_align(clean_dir,makeflow,samples,reflat,strand):
    fo =open(makeflow,"a+")
    align_dir = os.path.join(os.path.dirname(os.path.normpath(clean_dir)),"align")
    fo.write("CATEGORY=Stat_align\n")
    fo.write(os.path.join(align_dir,"align.state.txt : ") + " ".join(map(lambda x:os.path.join(align_dir,x,"bowtie.err"),samples)) + "\n")
    fo.write("\tpython /lustre/work/yongdeng/software/protokaryon/soft/Align/stat_align_result.py -a %s -s %s -o %s " %(align_dir,",".join(samples),os.path.join(align_dir,"align.state.txt")))
    fo.write("> %s 2> %s\n\n"%(os.path.join(align_dir,"stat_align_result.out"),os.path.join(align_dir,"stat_align_result.err")))     
    for sn in samples:
        picard_dir = os.path.join(os.path.dirname(os.path.normpath(clean_dir)),"picard",sn)
        fo.write("CATEGORY=picard_%s\n"%sn)
        fo.write(os.path.join(picard_dir,"rna_seq_metrics.txt : ") + os.path.join(align_dir,sn,"accepted_hits.bam %s\n" % reflat))
        fo.write("@BATCH_OPTIONS = -l h_vmem=50G\n")
        fo.write("\tjava -XX:ParallelGCThreads=12 -Xmx30g -Djava.io.tmpdir=%s -jar /lustre/work/yongdeng/software/protokaryon/soft/Align/picard-2.14.0/picard.jar CollectRnaSeqMetrics REF_FLAT=%s STRAND=%s MINIMUM_LENGTH=200 I=%s O=%s VALIDATION_STRINGENCY=LENIENT "%(os.path.join(picard_dir,"tmp"),reflat,strand,os.path.join(align_dir,sn,"accepted_hits.bam"),os.path.join(picard_dir,"rna_seq_metrics.txt")))
        fo.write("> %s 2> %s\n\n"%(os.path.join(picard_dir,"CollectRnaSeqMetrics.out"),os.path.join(picard_dir,"CollectRnaSeqMetrics.err")))
        fo.write("CATEGORY=parse_rna_seq_metrics_%s\n"%sn)
        fo.write(os.path.join(picard_dir,"element_distribution.txt ") + os.path.join(picard_dir,"uniform_distribution.txt : ") + os.path.join(picard_dir,"rna_seq_metrics.txt\n"))
        fo.write("\tpython /lustre/work/yongdeng/software/protokaryon/soft/Align/parse_rna_seq_metrics.py -i %s -u %s -e %s &> %s\n\n"%(os.path.join(picard_dir,"rna_seq_metrics.txt"),os.path.join(picard_dir,"uniform_distribution.txt"),os.path.join(picard_dir,"element_distribution.txt"),os.path.join(picard_dir,"ele_unif_distribution.log")))
        fo.write("CATEGORY=element_distribution_png_%s\n"%sn)
        fo.write(os.path.join(picard_dir,"element_distribution.png : ") + os.path.join(picard_dir,"element_distribution.txt\n"))
        fo.write("\tRscript /lustre/work/yongdeng/software/protokaryon/soft/Align/Element_dis.r %s %s &> %s\n\n"%(os.path.join(picard_dir,"element_distribution.txt"),os.path.join(picard_dir,"element_distribution.png"),os.path.join(picard_dir,"Element_dis.rscript.log")))
        fo.write("CATEGORY=uniform_distribution_png_%s\n"%sn)
        fo.write(os.path.join(picard_dir,"uniform_distribution.png : ") + os.path.join(picard_dir,"uniform_distribution.txt\n"))
        fo.write("\tRscript /lustre/work/yongdeng/software/protokaryon/soft/Align/Uniform_dis.r %s %s &> %s\n\n"%(os.path.join(picard_dir,"uniform_distribution.txt"),os.path.join(picard_dir,"uniform_distribution.png"),os.path.join(picard_dir,"Uniform_dis.rscript.log")))
        fo.write("CATEGORY=saturation_%s\n"%sn)
        fo.write(os.path.join(picard_dir,"rnaseq.saturation : ") + os.path.join(align_dir,sn,"accepted_hits.bam\n"))
        fo.write("@BATCH_OPTIONS = -l h_vmem=50G\n")
        fo.write("\tperl /lustre/work/yongdeng/software/protokaryon/soft/Align/RNA-Seq_saturation_fastV3_lessMEM.pl -s %s -sd 5 -g %s -o %s "%(os.path.join(align_dir,sn,"accepted_hits.bam"),reflat.partition(".refFlat")[0],os.path.join(picard_dir,"rnaseq.saturation")))
        fo.write("> %s 2> %s\n\n"%(os.path.join(picard_dir,"RNA-Seq_saturation.out"),os.path.join(picard_dir,"RNA-Seq_saturation.err")))
    fo.close()

def write_avgQual_makeflow(samples,raw_dir,clean_dir,makeflow):
    fo =open(makeflow,"a+")
    for sn in samples:
        for r in "12":
            fo.write("CATEGORY=avgQual_%s\n"%sn)
            fo.write(os.path.join(raw_dir,sn,"R%s_raw_avgQual.txt : "%r) + os.path.join(raw_dir,sn,"R%s.fastq.gz\n"%r))
            fo.write("\tpython /lustre/work/yongdeng/software/protokaryon/soft/QC/avgQual.py -i %s -o %s "%(os.path.join(raw_dir,sn,"R%s.fastq.gz"%r),os.path.join(raw_dir,sn,"R%s_raw_avgQual.txt"%r)))
            fo.write("&> %s\n\n"%os.path.join(raw_dir,sn,"R%s_raw_avgQual.log"%r))            
            fo.write("CATEGORY=avgQual_%s\n"%sn)
            fo.write(os.path.join(clean_dir,sn,"R%s_clean_avgQual.txt : "%r) + os.path.join(clean_dir,sn,"R%s.clean.fq.gz\n"%r))
            fo.write("\tpython /lustre/work/yongdeng/software/protokaryon/soft/QC/avgQual.py -i %s -o %s "%(os.path.join(clean_dir,sn,"R%s.clean.fq.gz"%r),os.path.join(clean_dir,sn,"R%s_clean_avgQual.txt"%r)))
            fo.write("&> %s\n\n"%os.path.join(clean_dir,sn,"R%s_clean_avgQual.log"%r))
    fo.close()

#def write_illuqc_stat_clean_makeflow(samples,raw_dir,clean_dir,makeflow):

def write_readcount_makeflow(samples,makeflow,input_d,output_d,gtf,strand):
    fo =open(makeflow,"a+")
    for sn in samples:
        od = os.path.join(output_d,sn)
        count = open(os.path.join(od,"%s_count.txt"%sn),"w")
        count.write("gene_id\t%s\n"%sn)
        count.close()
        fo.write("CATEGORY=Sam-sort_%s\n"%sn)
        fo.write("%s : %s\n" %(os.path.join(input_d,sn,"accepted_hits.bam"),os.path.join(input_d,sn,"accepted_hits.sam")))
        fo.write("@BATCH_OPTIONS = -l h_vmem=50G\n")
        fo.write("\t/lustre/work/yongdeng/software/protokaryon/soft/readcount/samtools-1.6/samtools sort -n -@ 8 -o %s -O bam -T %s %s " %(os.path.join(input_d,sn,"accepted_hits.bam"),os.path.join(input_d,sn,"temp"),os.path.join(input_d,sn,"accepted_hits.sam")))
        fo.write("> %s 2> %s && rm -fr %s\n\n" %(os.path.join(input_d,sn,"sort.out"),os.path.join(input_d,sn,"sort.err"),os.path.join(input_d,sn,"accepted_hits.sam")))        
        fo.write("CATEGORY=HTseq-count_%s\n"%sn)
        fo.write(os.path.join(od,"%s_count.txt : "%sn) + os.path.join(input_d,sn,"accepted_hits.bam\n"))
        fo.write("\t/home/yongdeng/software/python27/bin/htseq-count -s %s -f bam %s %s >> %s " %(strand,os.path.join(input_d,sn,"accepted_hits.bam"),gtf,os.path.join(od,"%s_count.txt"%sn)))
        fo.write("2> %s\n\n" %os.path.join(od,"%s_count.progress"%sn))
    fo.write("CATEGORY=get_count_matrix\n")
    fo.write(os.path.join(output_d,"gene_count_matrix.txt : ") + " ".join([os.path.join(output_d,i,i+"_count.txt") for i in samples]) + "\n")
    fo.write("\tpython /home/yongdeng/software/python_scripts/get_count_matrix.py %s %s &> %s\n\n"%(" ".join([os.path.join(output_d,i,i+"_count.txt") for i in samples]),os.path.join(output_d,"gene_count_matrix.txt"),os.path.join(output_d,"gene_count_matrix.log")))
    fo.close()

def write_deg_analysis_makeflow(makeflow,countable,samples,groups,vs_list,fc,outdir):
    fo = open(makeflow,"a+")
    fo.write("CATEGORY=get_deg_rscript\n")
    fo.write(os.path.join(outdir,"diff_gene_DESeq2.r : ") + countable + "\n")
    fo.write("\tpython /home/yongdeng/software/python_scripts/get_deg.py -c %s -fc %s -s %s -g %s -vs %s -o %s "%(countable,fc,",".join(samples),",".join(groups),",".join(vs_list),outdir))
    fo.write('> %s 2> %s\n\n'%(os.path.join(outdir,"get_deg_rscript.out"),os.path.join(outdir,"get_deg_rscript.err")))
    fo.write("CATEGORY=deg_analysis\n")
    fo.write(os.path.join(outdir,"gene.normalized.count.table ")  + os.path.join(outdir,"gene.pseudo.log2.counts.table ") + " ".join([os.path.join(outdir,i,i+".txt") for i in vs_list]) + " : " + os.path.join(outdir,"diff_gene_DESeq2.r ") + countable + "\n")
    fo.write("\tRscript %s > %s 2> %s\n\n"%(os.path.join(outdir,"diff_gene_DESeq2.r"),os.path.join(outdir,"diff_gene_DESeq2.out"),os.path.join(outdir,"diff_gene_DESeq2.err")))
    fo.close()
    
def write_stat_diff_gene_makeflow(makeflow,vs_list,outdir):
    fo = open(makeflow,"a+")
    fo.write("CATEGORY=get_scatter_valcon\n")
    fo.write(" ".join([os.path.join(outdir,i,i+j) for i in vs_list for j in [".diff_info",".down.diff",".plot_scatter.file",".plot_valcano.file",".up.diff"]]) + " : " + " ".join([os.path.join(outdir,i,i+".txt") for i in vs_list]) + "\n")
    fo.write("\tpython /home/yongdeng/software/python_scripts/pick_diffexp.py %s &> %s\n\n"%(" ".join([os.path.join(outdir,i,i+".txt") for i in vs_list]),os.path.join(outdir,"pick_diffexp.log")))
    for vs in vs_list:
        fo.write("CATEGORY=get_scatter_valcon_%s\n"%vs)
        fo.write(os.path.join(outdir,vs,vs+".plot_scatter.png : ") + os.path.join(outdir,vs,vs+".plot_scatter.file\n"))
        fo.write("\tRscript /lustre/work/yongdeng/software/protokaryon/soft/deseq/Scatter_volcano.r %s %s xy &> %s\n\n"%(os.path.join(outdir,vs,".plot_scatter.file"),os.path.join(outdir,vs,vs+".plot_scatter.png"),os.path.join(outdir,vs,"Scatter_rscript.log")))
        fo.write("CATEGORY=get_scatter_valcon_%s\n"%vs)
        fo.write(os.path.join(outdir,vs,vs+".plot_valcano.png : ") + os.path.join(outdir,vs,vs+".plot_valcano.file\n"))
        fo.write("\tRscript /lustre/work/yongdeng/software/protokaryon/soft/deseq/Scatter_volcano.r %s %s &> %s\n\n"%(os.path.join(outdir,vs,".plot_valcano.file"),os.path.join(outdir,vs,vs+".plot_valcano.png"),os.path.join(outdir,vs,"Volcano_rscript.log")))
    fo.write("CATEGORY=stat_diff\n")
    fo.write(" ".join([os.path.join(outdir,vs,vs+".diff_stat") for vs in vs_list]) + " : " + " ".join([os.path.join(outdir,i,i+".txt") for i in vs_list]) + "\n")
    fo.write("\tpython /home/yongdeng/software/python_scripts/stat_diffexp.py %s &> %s\n\n"%(" ".join([os.path.join(outdir,i,i+".txt") for i in vs_list]),os.path.join(outdir,"stat_diffexp.log")))     
    fo.close()    

def write_back_deg_list_makeflow(makeflow,vs_list,anno,outdir,dbtype,reftype):
    fo = open(makeflow,"a+")
    fo.write("CATEGORY=get_deg_back_list\n")
    fo.write(os.path.join(outdir,"back.list ") + " ".join([os.path.join(outdir,vs,j) for vs in vs_list for j in ["deg.list","deg.num.txt","deg.up.txt","deg.down.txt","deg.up_down.list"]]) + " : " + anno + " " + " ".join([os.path.join(outdir,i,i+".txt") for i in vs_list]) + "\n")
    fo.write("\tpython /home/yongdeng/software/python_scripts/call_deg_back_list.py -anno %s -t %s -r %s -diff %s &> %s\n\n"%(anno,dbtype,reftype," ".join([os.path.join(outdir,i,i+".txt") for i in vs_list]),os.path.join(outdir,"call_deg_back_list.log")))
    fo.close()
    
def write_stat_deg_num_makeflow(makeflow,outdir,vs_list):
    fo = open(makeflow,"a+")
    fo.write("CATEGORY=deg_num_stat\n")
    fo.write(os.path.join(outdir,"deg.state.txt : ") + " ".join([os.path.join(outdir,vs,"deg.num.txt") for vs in vs_list]) + "\n")
    fo.write("\tpython /home/yongdeng/software/python_scripts/stat_deg_num.py %s %s &> %s\n\n"%(" ".join([os.path.join(outdir,vs,"deg.num.txt") for vs in vs_list]),os.path.join(outdir,"deg.state.txt"),os.path.join(outdir,"stat_deg_num.log")))
    fo.write("CATEGORY=deg_stat_plot\n")
    fo.write(os.path.join(outdir,"deg.stat.png ") + os.path.join(outdir,"deg.stat.pdf : ") + os.path.join(outdir,"deg.state.txt\n"))
    fo.write("\tRscript /lustre/work/yongdeng/software/protokaryon/soft/deseq/deg_stat.r %s %s %s &> %s\n\n"%(os.path.join(outdir,"deg.state.txt"),os.path.join(outdir,"deg.stat.png"),os.path.join(outdir,"deg.stat.pdf"),os.path.join(outdir,"deg_stat.log")))
    fo.close()

def cp_genelist_cmd(genes,genes_up_down,subdir):
    l = []
    for s,d in zip(genes,subdir):
        l.append("cp " + s + " " + d + "/deg.list")
    for s,d in zip(genes_up_down,subdir):
        l.append("cp " + s + " " + d + "/deg.up_down.list")
    return " && ".join(l)
    
def write_cp_gene_list_makeflow(makeflow,outputdir,back,gene,genes_up_down,subdir):
    fo = open(makeflow,"a+")
    fo.write("CATEGORY=cp_back_deg_list\n")
    fo.write(os.path.join(outputdir,"back","back.list ") + " ".join([os.path.join(outputdir,vs,i) for vs in subdir for i in ["deg.list","deg.up_down.list"]]) + " : " + back  + " " + " ".join(gene) + " " + " ".join(genes_up_down) + "\n")
    cp_list = ["cp " + back + " " + os.path.join(outputdir,"back","back.list")]
    for i in range(len(gene)):
        cp_list.append("cp " + gene[i] + " " + subdir[i])
    for i in range(len(genes_up_down)):
        cp_list.append("cp " + genes_up_down[i] + " " + subdir[i])
    fo.write("\t%s &> %s\n\n"%(" && ".join(cp_list),os.path.join(outputdir,"cp_deg_back_list.log")))
    fo.close()
    
def write_kegg_anno_makeflow(makeflow,annolist,subdir,idtype,abbr):
    fo = open(makeflow,"a+")
    out = os.path.join(subdir,os.path.splitext(os.path.basename(annolist))[0])
    fo.write("CATEGORY=annotate_%s\n"%os.path.basename(subdir))
    fo.write( out + ".ko : " + annolist + "\n")
    fo.write("\tpython /lustre/software/target/kobas-2.0/src/annotate.py -i %s -t %s -s %s -o %s "%(annolist,idtype,abbr,out + ".ko"))
    fo.write("> %s 2> %s\n\n"%(os.path.join(out +".annotate.out"),os.path.join(out +".annotate.err")))
    fo.write("CATEGORY=annotate_parse_%s\n"%os.path.basename(subdir))
    fo.write(" ".join([os.path.join(subdir) + "/KEGG_PATHWAY" + i for i in [".first.class",".second.class",".third.class",".first.list",".second.list",".third.list : "]]) + out + ".ko\n")
    fo.write("\tperl /lustre/work/yongdeng/software/protokaryon/soft/gokegg/KobasAnnoForClassify.pl -KobasAnno %s -KEGGClassify %s -od %s "%(out+".ko","/lustre/work/yongdeng/software/protokaryon/soft/gokegg/KEGG.pathway.classify.list",subdir))
    fo.write("&> %s\n\n"%os.path.join(subdir,"kegg.classify.log"))
    fo.write("CATEGORY=kegg_plot_%s\n"%os.path.basename(subdir))
    fo.write(os.path.join(subdir,"KEGG_PATHWAY.second.class.png ") + os.path.join(subdir,"KEGG_PATHWAY.second.class.pdf : ") + os.path.join(subdir,"KEGG_PATHWAY.second.class\n"))
    fo.write("\tRscript /lustre/work/yongdeng/software/protokaryon/soft/gokegg/kegg-plot.R %s %s %s "%(os.path.join(subdir,"KEGG_PATHWAY.second.class"),os.path.join(subdir,"KEGG_PATHWAY.second.class.png"),os.path.join(subdir,"KEGG_PATHWAY.second.class.pdf")))
    fo.write("&> %s\n\n"%os.path.join(subdir,"KEGG_PATHWAY.second.class.plot.log"))
    fo.close()
    
def write_kegg_ident_makeflow(makeflow,subdir,dbtype,fdr="BH"):
    dbtype = dbtype.replace("/G","")
    fo = open(makeflow,"a+")
    ko = os.path.join(subdir,"deg.ko")
    back_ko = os.path.join( os.path.dirname(subdir),"back","back.ko")
    fo.write("CATEGORY=kegg_identity_%s\n"%os.path.basename(subdir))
    fo.write(os.path.join(subdir,"deg.pathway.enrichment : ") + ko + " " + back_ko + "\n")
    fo.write("\tpython /lustre/software/target/kobas-2.0/bin/identify.py -c 1 -f %s -b %s -d %s -o %s -n %s "%(ko,back_ko,dbtype,os.path.join(subdir,"deg.pathway.enrichment"),fdr))
    fo.write("> %s 2> %s\n\n"%(os.path.join(subdir,"identify.out"),os.path.join(subdir,"identify.err")))
    fo.close()
    
def write_go_enrichment_makeflow(makeflow,annolist,subdir,output_dir,dbtype,abbr,idtype):
    fo = open(makeflow,"a+")
    back_dir = os.path.join(output_dir,"back")
    go_info = os.popen("grep '^%s' /lustre/work/yongdeng/software/protokaryon/soft/gokegg/specise_v1.1 | awk '{print $NF}'"%abbr).read().strip()
    if go_info == "kobas_nogo" and idtype != "id:ncbigi":
        fo.write("CATEGORY=prepare_go\n")
        fo.write(os.path.join(back_dir,"back.ref.idmapping : ") + os.path.join(back_dir,"back.list\n"))
        fo.write("\tperl /lustre/work/yongdeng/software/protokaryon/soft/gokegg/RNA-Seq_Anno1-3_id_mapping_byfile.pl -r -m /lustre/work/yongdeng/software/protokaryon/soft/gokegg/%s.dbfile -p back. -o %s %s "%(idtype,back_dir,os.path.join(back_dir,"back.list")))
        fo.write("> %s 2> %s\n\n"%(os.path.join(back_dir,"RNA-Seq_Anno1-3_id_mapping_byfile.out"),os.path.join(back_dir,"RNA-Seq_Anno1-3_id_mapping_byfile.err")))
        fo.write("CATEGORY=prepare_go\n")
        fo.write(os.path.join(back_dir,"back_go.idmapping : ") + os.path.join(back_dir,"back.ref.idmapping\n"))
        fo.write("\tperl /lustre/work/yongdeng/software/protokaryon/soft/gokegg/idmapping_remove_ko.pl -in %s -od %s &> %s\n\n"%(os.path.join(back_dir,"back.ref.idmapping"),back_dir,os.path.join(back_dir,"idmapping_remove_ko.log")))
        fo.write("CATEGORY=prepare_go\n")
        fo.write(os.path.join(back_dir,"back.go.normalization : ") + os.path.join(back_dir,"back_go.idmapping ") + "/lustre/work/yongdeng/software/protokaryon/soft/gokegg/go-basic.obo\n")
        fo.write("\tperl /lustre/work/yongdeng/software/protokaryon/soft/gokegg/go_normalization.pl /lustre/work/yongdeng/software/protokaryon/soft/gokegg/go-basic.obo %s %s "%(os.path.join(back_dir,"back_go.idmapping"),os.path.join(back_dir,"back.go.normalization")))
        fo.write("> %s\n\n"%os.path.join(back_dir,"back.go_normalization.log"))
    elif go_info == "kobas_havego":    
        fo.write("CATEGORY=prepare_go\n")
        fo.write(os.path.join(back_dir,"back.go.normalization : ") + os.path.join(back_dir,"back.tsv ") + os.path.join(back_dir,"back.go ") + "/lustre/work/yongdeng/software/protokaryon/soft/gokegg/go-basic.obo\n")
        fo.write("\tperl /lustre/work/yongdeng/software/protokaryon/soft/gokegg/go_normalization.pl /lustre/work/yongdeng/software/protokaryon/soft/gokegg/go-basic.obo %s %s "%(os.path.join(back_dir,"back.go"),os.path.join(back_dir,"back.go.normalization")))
        fo.write("&> %s\n\n"%os.path.join(back_dir,"back.go_normalization.log"))
    else:
        print "GO enrichment can not use %s type"%idtype
        sys.exit()
    fo.write("CATEGORY=prepare_go\n")
    fo.write(os.path.join(back_dir,"back.tsv ") + os.path.join(back_dir,"back.go : ") + os.path.join(back_dir,"back.ko ") + "/lustre/work/yongdeng/software/protokaryon/soft/gokegg/go-basic.obo\n")
    fo.write("\tperl /lustre/work/yongdeng/software/protokaryon/soft/gokegg/parse_kobas_annotate.pl %s %s %s "%(os.path.join(back_dir,"back.ko"),"/lustre/work/yongdeng/software/protokaryon/soft/gokegg/go-basic.obo",os.path.join(back_dir,"back")))
    fo.write("&> %s\n\n"%os.path.join(back_dir,"parse_kobas_annotate.log"))        
    fo.write("CATEGORY=back_anno\n")
    fo.write(os.path.join(back_dir,"back.anno.xls : ") + os.path.join(back_dir,"back.tsv ") + os.path.join(back_dir,"back.go.normalization ") + "/lustre/work/yongdeng/software/protokaryon/soft/gokegg/go-basic.obo\n")
    fo.write("\tperl /lustre/work/yongdeng/software/protokaryon/soft/gokegg/ForKobasTsv.pl -tsv %s -od %s -database %s -go %s -obo /lustre/work/yongdeng/software/protokaryon/soft/gokegg/go-basic.obo "%(os.path.join(back_dir,"back.tsv"),back_dir,dbtype,os.path.join(back_dir,"back.go.normalization")))
    fo.write("&> %s\n\n"%os.path.join(back_dir,"back.anno.log"))
    
    for deg,dir in zip(annolist,subdir):   ### annolist为各样本分组的deg.up_down.list文件，subdir为group1-VS-group2目录
        dir_name = os.path.basename(dir)
        out = os.path.join(output_dir,dir)
        fo.write("CATEGORY=go_enrichment_%s\n"%dir_name)
        fo.write(os.path.join(out,"deg.go.class ") + os.path.join(out,"deg.go.enrichment : ") + "/lustre/work/yongdeng/software/protokaryon/soft/gokegg/go-basic.obo " + deg + " " + os.path.join(back_dir,"back.go.normalization\n"))
        fo.write("\tperl /lustre/work/yongdeng/software/protokaryon/soft/gokegg/RNA-Seq_Anno3_enrichmentV1.pl -g /lustre/work/yongdeng/software/protokaryon/soft/gokegg/go-basic.obo -o %s -q deg. %s %s "%(out,os.path.join(back_dir,"back.go.normalization"),deg))
        fo.write(" > %s 2> %s\n\n"%(os.path.join(out,"RNA-Seq_Anno3_enrichmentV1.out"),os.path.join(out,"RNA-Seq_Anno3_enrichmentV1.err")))
    fo.close()
    
def write_go_pathway_summary_makeflow(makeflow,subdir):
    fo = open(makeflow,"a+")
    for vs in subdir:
        vs_name = os.path.basename(vs)
        go_enrich = os.path.join(vs,"deg.go.enrichment")
        path_enrich = os.path.join(vs,"deg.pathway.enrichment")
        fo.write("CATEGORY=ForKobasDot_%s\n"%vs_name)
        fo.write(os.path.join(vs,"enrichment.xls : ") + go_enrich + " " + path_enrich + "\n")
        fo.write("\tperl /lustre/work/yongdeng/software/protokaryon/soft/gokegg/ForKobasDot_v1.1.pl -ge %s -ke %s -o %s &> %s\n\n"%(go_enrich,path_enrich,vs,os.path.join(vs,"ForKobasDot.log")))
        fo.write("CATEGORY=parse_identify_%s\n"%vs_name)
        fo.write(" ".join([os.path.join(vs,"enrichment.") + i for i in ["GO.hierarchy.molecular_function.dot","GO.hierarchy.cellular_component.dot","GO.hierarchy.biological_process.dot","GO.txt","Pathway.txt","Disease.txt : "]]) + os.path.join(vs,"enrichment.xls\n"))
        fo.write("\tperl /lustre/work/yongdeng/software/protokaryon/soft/gokegg/parse_kobas_identify.pl %s /lustre/work/yongdeng/software/protokaryon/soft/gokegg/go-basic.obo 30 %s "%(os.path.join(vs,"enrichment.xls"),os.path.join(vs,"enrichment")))
        fo.write("> %s 2> %s\n\n"%(os.path.join(vs,"parse_kobas_identify.out"),os.path.join(vs,"parse_kobas_identify.err")))
        for level in ["biological_process","molecular_function","cellular_component"]:
            for p in ["png","pdf"]:
                fo.write("CATEGORY=enrichment.GO_%s_%s\n"%(level,p))
                fo.write(os.path.join(vs,"enrichment.GO.hierarchy.") + level + ".%s : "%p + os.path.join(vs,"enrichment.GO.hierarchy.") + level + ".dot\n")
                fo.write("\tdot -T" + p + " -o " + os.path.join(vs,"enrichment.GO.hierarchy.") + level + ".%s "%p + os.path.join(vs,"enrichment.GO.hierarchy.") + level + ".dot &> %s\n\n"%os.path.join(vs,"enrichment.GO.hierarchy." + level + "." + p + ".log"))
        fo.write("CATEGORY=GO_pvalue_%s\n"%vs_name)
        fo.write(os.path.join(vs,"enrichment.GO.p_value.top.pdf ") + os.path.join(vs,"enrichment.GO.p_value.top.png : ") + os.path.join(vs,"enrichment.GO.txt\n"))
        fo.write("\tRscript /lustre/work/yongdeng/software/protokaryon/soft/gokegg/draw_barplot_pvalue.r GO 30 %s %s %s &> %s\n\n"%(os.path.join(vs,"enrichment.GO.txt"),os.path.join(vs,"enrichment.GO.p_value"),os.path.join(vs,"enrichment.GO.p_value"),os.path.join(vs,"enrichment.GO.p_value.log")))
        fo.write("CATEGORY=GO_qvalue_%s\n"%vs_name)
        fo.write(os.path.join(vs,"enrichment.GO.q_value.top.pdf ") + os.path.join(vs,"enrichment.GO.q_value.top.png : ") + os.path.join(vs,"enrichment.GO.txt\n"))
        fo.write("\tRscript /lustre/work/yongdeng/software/protokaryon/soft/gokegg/draw_barplot.r GO 30 %s %s %s &> %s\n\n"%(os.path.join(vs,"enrichment.GO.txt"),os.path.join(vs,"enrichment.GO.q_value"),os.path.join(vs,"enrichment.GO.q_value"),os.path.join(vs,"enrichment.GO.q_value.log")))
    fo.close()
    
def write_kobas_color_makeflow(makeflow,subdir,anno,idtype):
    fo = open(makeflow,"a+")
    for vs in subdir:
        vs_name = os.path.basename(vs)
        dp_down = os.path.join(vs,"deg.up_down.list")
        for e in ["GO","Pathway"]:            
            fo.write("CATEGORY=kobas_color_%s_%s\n"%(e,vs_name))
            fo.write(os.path.join(vs,"%s.temp.xls : "%e) + os.path.join(vs,"enrichment.%s.txt\n"%e))
            fo.write("\tperl /lustre/work/yongdeng/software/protokaryon/soft/gokegg/kobas_colour1.6.pl -deg %s -k_anno %s -k_iden %s -inf %s "%(dp_down,os.path.join(vs,"deg.ko"),os.path.join(vs,"enrichment.%s.txt"%e),anno))
            fo.write("-S symbol -Y %s -od %s -p %s.temp -gc 1 -ud 3 > %s 2> %s\n\n" %(idtype,vs,e,os.path.join(vs,"enrichment.%s.change.out"%e),os.path.join(vs,"enrichment.%s.change.err"%e)))
            fo.write("CATEGORY=enrichment_%s_xls_%s\n"%(e,vs_name))
            fo.write(os.path.join(vs,"enrichment.%s.xls : "%e) + os.path.join(vs,"%s.temp.xls\n"%e))
            fo.write("\tcut -f 1-13 %s > %s\n\n"%(os.path.join(vs,"%s.temp.xls"%e),os.path.join(vs,"enrichment.%s.xls"%e)))
    fo.close()
    
        