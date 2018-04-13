#!/usr/bin/env python
#coding:utf-8
### 质控步骤
## cd /lustre/work/yongdeng/Project && python ~/Htseq_count.py -i  160761_test/align -g 160761_test/genome/GCF_000934625.1_ASM93462v1_genomic.gtf.new.gtf -s no -o 160761_test/Htseq -mf 160761_test/Htseq_count.mf
import argparse,os,sys
import RNA

parser = argparse.ArgumentParser(description="This Script is used to count the reads of gene and trans in all sample by HTSeq-count program.")
parser.add_argument("-i","--input_dir",type=str,help="The input aligned directory",required = True)
parser.add_argument("-g","--gtf",type=str,help="The reference gtf file",required = True)
parser.add_argument("-s","--stranded",type=str,choices=["yes","no","reverse"],default="no",help="whether the data is from a strand-specific assay. Specify 'yes', 'no', or 'reverse' (default: no). 'reverse' means 'yes' with reversed strand interpretation",required = False)
parser.add_argument("-o","--output_dir",type=str,help="The output read-count directory",required = True)
parser.add_argument("-v","--version",action="version",version='%(prog)s 1.0')
parser.add_argument("-p","--parallel",type=int,help="the cpu number used",required = False, default = 8)
parser.add_argument("-mf","--makeflow",type=str,help="the out makeflow file",required = True)
args = parser.parse_args()

def mkdir(dir,sampls):
    for n in sampls:
        d = os.path.join(dir,n)
        if not os.path.exists(d):
            os.makedirs(d)
        else:
            continue

def write_readcount_makeflow(samples,makeflow,input_d,output_d,gtf,strand,thread):
    samtools = RNA.get_bin_abspath("samtools")
    read_count = RNA.get_bin_abspath("read_count.py")
    get_count_matrix = RNA.get_bin_abspath("get_count_matrix.py")
    fo =open(makeflow,"a+")
    for sn in samples:
        od = os.path.join(output_d,sn)
        if not os.path.exists(os.path.join(input_d,sn,"accepted_hits.bam")):
            fo.write("CATEGORY=Sam-sort_%s\n"%sn)
            fo.write("%s : %s\n" %(os.path.join(input_d,sn,"accepted_hits.bam"),os.path.join(input_d,sn,"accepted_hits.sam")))
            fo.write("@BATCH_OPTIONS = -l h_vmem=50G\n")
            fo.write("\t%s sort -n -@ %d -o %s -O bam -T %s %s " %(samtools,thread,os.path.join(input_d,sn,"accepted_hits.bam"),os.path.join(input_d,sn,"temp"),os.path.join(input_d,sn,"accepted_hits.sam")))
            fo.write("> %s 2> %s && rm -fr %s\n\n" %(os.path.join(input_d,sn,"sort.out"),os.path.join(input_d,sn,"sort.err"),os.path.join(input_d,sn,"accepted_hits.sam")))    
        fo.write("CATEGORY=HTseq_gene_count_%s\n"%sn)
        fo.write(os.path.join(od,"%s_gene_count.txt : "%sn) + os.path.join(input_d,sn,"accepted_hits.bam\n"))
        fo.write("\tpython %s -s %s -f bam %s %s > %s 2> %s\n\n" %(read_count,strand,os.path.join(input_d,sn,"accepted_hits.bam"),gtf,os.path.join(od,"%s_gene_count.txt"%sn),os.path.join(od,"%s_gene_count.progress"%sn)))
        fo.write("CATEGORY=HTseq_trans_count_%s\n"%sn)
        fo.write(os.path.join(od,"%s_trans_count.txt : "%sn) + os.path.join(input_d,sn,"accepted_hits.bam\n"))
        fo.write("\tpython %s -s %s -i transcript_id -f bam %s %s > %s 2> %s\n\n" %(read_count,strand,os.path.join(input_d,sn,"accepted_hits.bam"),gtf,os.path.join(od,"%s_trans_count.txt"%sn),os.path.join(od,"%s_trans_count.progress"%sn)))
    fo.write("CATEGORY=get_count_matrix\n")
    fo.write(os.path.join(output_d,"gene_count_matrix.txt : ") + " ".join([os.path.join(output_d,i,i+"_gene_count.txt") for i in samples]) + "\n")
    fo.write("\tpython %s %s %s &> %s\n\n"%(get_count_matrix," ".join([os.path.join(output_d,i,i+"_gene_count.txt") for i in samples]),os.path.join(output_d,"gene_count_matrix.txt"),os.path.join(output_d,"gene_count_matrix.log")))
    fo.write("CATEGORY=get_trans_matrix\n")
    fo.write(os.path.join(output_d,"trans_count_matrix.txt : ") + " ".join([os.path.join(output_d,i,i+"_trans_count.txt") for i in samples]) + "\n")
    fo.write("\tpython %s %s %s &> %s\n\n"%(get_count_matrix," ".join([os.path.join(output_d,i,i+"_trans_count.txt") for i in samples]),os.path.join(output_d,"trans_count_matrix.txt"),os.path.join(output_d,"trans_count_matrix.log")))
    fo.close()

def main():
    input_dir = os.path.abspath(args.input_dir)
    output_dir = os.path.abspath(args.output_dir)
    gtf = os.path.abspath(args.gtf)
    if not gtf:
        print "Error: No %d file found!" %gtf
        sys.exit(1)
    strand = args.stranded
    samples = sorted(map(os.path.basename,filter(os.path.isdir,[os.path.join(input_dir,x) for x in os.listdir(input_dir)])))
    mkdir(output_dir,samples)
    write_readcount_makeflow(samples,args.makeflow,input_dir,output_dir,gtf,strand,args.parallel)

if __name__=="__main__":
    main()
    
        
