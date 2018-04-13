#!/usr/bin/env python
#coding:utf-8
### cd /lustre/work/yongdeng/Project && python ~/summary_data.py -raw 160761_test/rawdata/ -clean 160761_test/QC/ -reflat 160761_test/genome/GCF_000934625.1_ASM93462v1_genomic.gtf.new.gtf.refFlat -m 160761_test/summary_data.mf
import argparse,os,re,sys,glob,commands
import RNA

parser = argparse.ArgumentParser(description="This Script is used to summarize raw clean align results.")
parser.add_argument("-raw","--RawDir",type=str,help="The input rawdata directory")
#parser.add_argument("-raw","--RawDir",type=str,help="The input rawdata directory",required = True)
parser.add_argument("-align","--AlignDir",type=str,help="The input align directory")
#parser.add_argument("-align","--AlignDir",type=str,help="The input align directory",required = True)
parser.add_argument("-qc_m","--qc_methods",type=str,help="The quality control methods you have done, default:'IlluQC'",required = False, default = "IlluQC",choices=["IlluQC","trimmomatic"])
parser.add_argument("-align_m","--align_methods",type=str,help="The align methods you have done, default:'hisat'",required = False, default = "hisat2",choices=["hisat2","bowtie2"])
parser.add_argument("-clean","--CleanDir",type=str,help="The input clean directory")
#parser.add_argument("-clean","--CleanDir",type=str,help="The input clean directory",required = True)
#parser.add_argument("-reflat","--RefFlat",type=str,help="The RefFlat file for picard",required = True)
parser.add_argument("-reflat","--RefFlat",type=str,help="The RefFlat file for picard")
parser.add_argument("-s","--Strand",type=str,choices=["yes","no","reverse"],default="no",help="whether the data is from a strand-specific assay. Specify 'yes', 'no', or 'reverse' (default: no). 'reverse' means 'yes' with reversed strand interpretation",required = False)
parser.add_argument("-v","--version",action="version",version='%(prog)s 1.0')
parser.add_argument("-p","--parallel",type=int,help="the cpu number used",required = False, default = 8)
parser.add_argument("-o","--output_dir",type=str,help="The output summary directory, three dirs 'fastqcRaw', 'fastqcClean', 'picard' will be creat under this directory. default: the parent directory of rawdata or cleandata or aligndir")
parser.add_argument("-mf","--makeflow",type=str,help="the out makeflow file",required = True)
args = parser.parse_args()

def mkdir(dir,sampls):
    for n in sampls:
        d = os.path.join(dir,n)
        if not os.path.exists(d):
            os.makedirs(d)
        else:
            continue

def write_fastqc_makeflow(raw_dir,clean_dir,outdir,makeflow,samples,thread):
    fastqc = RNA.get_soft_abspath("fastqc")
    fo =open(makeflow,"a+")
    for sn in samples:
        od_raw = os.path.join(outdir,"fastqcRaw",sn)
        od_clean = os.path.join(outdir,"fastqcClean",sn)
        raw_r1 = os.path.join(raw_dir,sn,"R1.fq.gz")
        raw_r2 = os.path.join(raw_dir,sn,"R2.fq.gz")
        clean_r1 = os.path.join(clean_dir,sn,"R1.clean.fq.gz")
        clean_r2 = os.path.join(clean_dir,sn,"R2.clean.fq.gz")
        if raw_dir:
            fo.write("CATEGORY=fastqc_Raw_%s\n"%sn)
            fo.write(od_raw + "/R1_fastqc.html " + od_raw + "/R2_fastqc.html : " + raw_r1 + " " + raw_r2 + "\n")
            fo.write("\t%s -o %s -f fastq -t %d --extract %s %s " %(fastqc,od_raw,thread,raw_r1,raw_r2))
            fo.write("> %s 2> %s\n\n"%(od_raw + "/fastqc.out",od_raw + "/fastqc.err"))        
        if clean_dir:
            fo.write("CATEGORY=fastqc_Clean_%s\n"%sn)
            fo.write(od_clean + "/R1.clean_fastqc.html " + od_clean + "/R2.clean_fastqc.html : " + clean_r1 + " " + clean_r2 + "\n")
            fo.write("\t%s -o %s -f fastq -t %d --extract %s %s "%(fastqc,od_clean,thread,clean_r1,clean_r2))
            fo.write("> %s 2> %s\n\n"%(od_clean + "/fastqc.out",od_clean + "/fastqc.err"))
    fo.close()    

def write_summary_align(clean_dir,align_dir,outdir,makeflow,samples,reflat,strand,thread,align_err_file):
    stat_align_result = RNA.get_bin_abspath("stat_align_result.py")
    picard = RNA.get_soft_abspath("picard")
    parse_rna_seq_metrics = RNA.get_bin_abspath("parse_rna_seq_metrics.py")
    Element_dis = RNA.get_bin_abspath("Element_dis.r")
    Uniform_dis = RNA.get_bin_abspath("Uniform_dis.r")
    saturation = RNA.get_bin_abspath("RNA-Seq_saturation_fastV3_lessMEM.pl")
    Saturation_png = RNA.get_bin_abspath("Saturation.r")
    fo =open(makeflow,"a+")
    fo.write("CATEGORY=Stat_align\n")
    fo.write(os.path.join(align_dir,"align.state.txt : ") + " ".join(map(lambda x:os.path.join(align_dir,x,align_err_file),samples)) + "\n")
    fo.write("\tpython %s -a %s -s %s -o %s " %(stat_align_result,align_dir,",".join(samples),os.path.join(align_dir,"align.state.txt")))
    fo.write("> %s 2> %s\n\n"%(os.path.join(align_dir,"stat_align_result.out"),os.path.join(align_dir,"stat_align_result.err")))     
    for sn in samples:
        if os.path.exists(os.path.join(align_dir,sn,"accepted_hits.bam")):
            bam = os.path.join(align_dir,sn,"accepted_hits.bam")
        elif not os.path.exists(os.path.join(align_dir,sn,"accepted_hits.bam")) and os.path.exists(os.path.join(align_dir,sn,"accepted_hits.sam")):
            bam = os.path.join(align_dir,sn,"accepted_hits.sam")
        else:
            bam = os.path.join(align_dir,sn,"accepted_hits.bam")
        picard_dir = os.path.join(outdir,"picard",sn)
        fo.write("CATEGORY=picard_%s\n"%sn)
        fo.write(os.path.join(picard_dir,"rna_seq_metrics.txt : ") + bam + " " + reflat + "\n")
        fo.write("@BATCH_OPTIONS = -l h_vmem=50G\n")
        fo.write("\tjava -XX:ParallelGCThreads=%d -Xmx30g -Djava.io.tmpdir=%s -jar %s CollectRnaSeqMetrics REF_FLAT=%s STRAND=%s MINIMUM_LENGTH=200 I=%s O=%s VALIDATION_STRINGENCY=LENIENT "%(thread,os.path.join(picard_dir,"tmp"),picard,reflat,strand,bam,os.path.join(picard_dir,"rna_seq_metrics.txt")))
        fo.write("> %s 2> %s\n\n"%(os.path.join(picard_dir,"CollectRnaSeqMetrics.out"),os.path.join(picard_dir,"CollectRnaSeqMetrics.err")))
        fo.write("CATEGORY=parse_rna_seq_metrics_%s\n"%sn)
        fo.write(os.path.join(picard_dir,"element_distribution.txt ") + os.path.join(picard_dir,"uniform_distribution.txt : ") + os.path.join(picard_dir,"rna_seq_metrics.txt\n"))
        fo.write("\tpython %s -i %s -u %s -e %s &> %s\n\n"%(parse_rna_seq_metrics,os.path.join(picard_dir,"rna_seq_metrics.txt"),os.path.join(picard_dir,"uniform_distribution.txt"),os.path.join(picard_dir,"element_distribution.txt"),os.path.join(picard_dir,"ele_unif_distribution.log")))
        fo.write("CATEGORY=element_distribution_png_%s\n"%sn)
        fo.write(os.path.join(picard_dir,"element_distribution.png : ") + os.path.join(picard_dir,"element_distribution.txt\n"))
        fo.write("\tRscript %s %s %s &> %s\n\n"%(Element_dis,os.path.join(picard_dir,"element_distribution.txt"),os.path.join(picard_dir,"element_distribution.png"),os.path.join(picard_dir,"Element_dis.rscript.log")))
        fo.write("CATEGORY=uniform_distribution_png_%s\n"%sn)
        fo.write(os.path.join(picard_dir,"uniform_distribution.png : ") + os.path.join(picard_dir,"uniform_distribution.txt\n"))
        fo.write("\tRscript %s %s %s &> %s\n\n"%(Uniform_dis,os.path.join(picard_dir,"uniform_distribution.txt"),os.path.join(picard_dir,"uniform_distribution.png"),os.path.join(picard_dir,"Uniform_dis.rscript.log")))
        fo.write("CATEGORY=saturation_%s\n"%sn)
        fo.write(os.path.join(picard_dir,"rnaseq.saturation : ") + bam + "\n")
        fo.write("@BATCH_OPTIONS = -l h_vmem=50G\n")
        fo.write("\tperl %s -s %s -sd 5 -g %s -o %s "%(saturation,bam,reflat.partition(".refFlat")[0],picard_dir))
        fo.write("> %s 2> %s\n\n"%(os.path.join(picard_dir,"RNA-Seq_saturation.out"),os.path.join(picard_dir,"RNA-Seq_saturation.err")))
        fo.write("CATEGORY=saturation_png_%s\n"%sn)
        fo.write(os.path.join(picard_dir,"saturation.png : ") + os.path.join(picard_dir,"rnaseq.saturation\n"))
        fo.write("\tRscript %s %s %s &> %s\n\n"%(Saturation_png,os.path.join(picard_dir,"rnaseq.saturation"),os.path.join(picard_dir,"saturation.png"),os.path.join(picard_dir,"saturation_png.log")))
    fo.close()

def write_avgQual_makeflow(samples,raw_dir,clean_dir,makeflow):
    avgQual = RNA.get_bin_abspath("avgQual.py")
    fo =open(makeflow,"a+")
    for sn in samples:
        for r in "12":
            fo.write("CATEGORY=avgQual_raw_R%s_%s\n"%(r,sn))
            fo.write(os.path.join(raw_dir,sn,"R%s_raw_avgQual.txt : "%r) + os.path.join(raw_dir,sn,"R%s.fq.gz\n"%r))
            fo.write("\tpython %s -i %s -o %s "%(avgQual,os.path.join(raw_dir,sn,"R%s.fq.gz"%r),os.path.join(raw_dir,sn,"R%s_raw_avgQual.txt"%r)))
            fo.write("&> %s\n\n"%os.path.join(raw_dir,sn,"R%s_raw_avgQual.log"%r))            
            fo.write("CATEGORY=avgQual_clean_R%s_%s\n"%(r,sn))
            fo.write(os.path.join(clean_dir,sn,"R%s_clean_avgQual.txt : "%r) + os.path.join(clean_dir,sn,"R%s.clean.fq.gz\n"%r))
            fo.write("\tpython %s -i %s -o %s "%(avgQual,os.path.join(clean_dir,sn,"R%s.clean.fq.gz"%r),os.path.join(clean_dir,sn,"R%s_clean_avgQual.txt"%r)))
            fo.write("&> %s\n\n"%os.path.join(clean_dir,sn,"R%s_clean_avgQual.log"%r))
    fo.close()

def write_illuqc_stat_clean_makeflow(makeflow,samples,clean_dir):
    ngsStat = RNA.get_bin_abspath("ngsStat.pl")
    fo = open(makeflow,"a+")
    fo.write("CATEGORY=stat_clean\n")
    fo.write(os.path.join(clean_dir,"clean_data.state.txt : ") + " ".join([os.path.join(clean_dir,sn,"Filter_stat") for sn in samples]) + "\n")
    fo.write("\tperl %s -dir %s -o %s &> %s\n\n"%(ngsStat,clean_dir,clean_dir,os.path.join(clean_dir,"stat_clean_data.log")))
    fo.close()
    
def write_trimmomatic_stat_clean_makeflow(makeflow,samples,clean_dir,raw_dir):
    stat_cleandata = RNA.get_bin_abspath("stat_cleandata.py")
    fo = open(makeflow,"a+")
    fo.write("CATEGORY=stat_clean\n")
    avgQual_filelist = [os.path.join(raw_dir,sn,"R"+n+"_raw_avgQual.txt") for sn in samples for n in "12"] + [os.path.join(clean_dir,sn,"R"+n+"_clean_avgQual.txt") for sn in samples for n in "12"]
    fo.write(os.path.join(clean_dir,"clean_data.state.txt : ") + " ".join(avgQual_filelist)  + "\n")
    fo.write("\tpython %s %s &> %s\n\n"%(stat_cleandata," ".join(avgQual_filelist),os.path.join(clean_dir,"stat_clean.log")))
    fo.close()

def main():
    reflat = os.path.abspath(args.RefFlat) if args.RefFlat else ""
    rawdir = os.path.abspath(args.RawDir) if args.RawDir else ""
    aligndir = os.path.abspath(args.AlignDir) if args.AlignDir else ""
    cleandir = os.path.abspath(args.CleanDir) if args.CleanDir else ""
    if not any([rawdir,aligndir,cleandir]):
        print "Error:rawdir or aligndir or cleandir must be given"
        sys.exit(1)
    if not args.output_dir:
        outdir  = os.path.dirname(os.path.normpath(rawdir))
    else:
        outdir = os.path.abspath(args.output_dir)
    if rawdir:
        samples = sorted(map(os.path.basename,filter(os.path.isdir,[os.path.join(rawdir,x) for x in os.listdir(rawdir)])))
        mkdir(outdir,[i + "/" + j for i in ["fastqcRaw"] for j in samples])
    if cleandir:
        samples = sorted(map(os.path.basename,filter(os.path.isdir,[os.path.join(cleandir,x) for x in os.listdir(cleandir)])))
        mkdir(outdir,[i + "/" + j for i in ["fastqcClean"] for j in samples])
    if aligndir:
        samples = sorted(map(os.path.basename,filter(os.path.isdir,[os.path.join(aligndir,x) for x in os.listdir(aligndir)])))
        mkdir(outdir,[i + "/" + j for i in ["picard"] for j in samples])
    write_fastqc_makeflow(rawdir,cleandir,outdir,args.makeflow,samples,args.parallel)
    if cleandir and args.qc_methods == "trimmomatic":
        write_avgQual_makeflow(samples,rawdir,cleandir,args.makeflow)
        write_trimmomatic_stat_clean_makeflow(args.makeflow,samples,cleandir,rawdir)
    elif cleandir and args.qc_methods == "IlluQC":
        write_illuqc_stat_clean_makeflow(args.makeflow,samples,cleandir)
    strand = {"yes":"FIRST_READ_TRANSCRIPTION_STRAND","no":"NONE","reverse":"SECOND_READ_TRANSCRIPTION_STRAND"}
    if aligndir and args.align_methods == "bowtie2":
        align_err_file = 'bowtie.err'
    elif aligndir and args.align_methods == "hisat2":
        align_err_file = "hisat2.err"
    if aligndir:
        if reflat:
            write_summary_align(cleandir,aligndir,outdir,args.makeflow,samples,reflat,strand[args .Strand],args.parallel,align_err_file)
        else:
            print "Error:No refFlat file, please give"
            sys.exit(1)
    
if __name__=="__main__":
    main()        

    
    
    
