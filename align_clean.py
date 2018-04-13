#!/usr/bin/env python
#coding:utf-8
### 比对步骤
## cd /lustre/work/yongdeng/Project && python ~/align_clean.py -i 160761_test/QC/ -g 160761_test/genome/ -mf 160761_test/align_clean.mf -o 160761_test/align
import argparse,os,sys,glob
import RNA

parser = argparse.ArgumentParser(description="This Script is used to align cleandata to reference genome by bowtie2 program.")
parser.add_argument("-i","--input_dir",type=str,help="The input cleandata directory",required = True)
parser.add_argument("-g","--genome_dir",type=str,help="the genome reference dir, 'fa' and 'gtf' file must under this dir, and align index will be generate in this dir",required=True)
parser.add_argument("-p","--parallel",type=int,help="the cpu number used",required = False, default = 8)
parser.add_argument("-o","--output_dir",type=str,help="The output ailgned directory",required = True)
parser.add_argument("-s","--strand",action = "store_true",default= False,help="whether the data is from a strand-specific assay, default: False if no set")
parser.add_argument("-v","--version",action="version",version='%(prog)s 1.0')
parser.add_argument("-mf","--makeflow",type=str,help="the out makeflow file",required = True)
args = parser.parse_args()

def mkdir(dir,sampls):
    for n in sampls:
        d = os.path.join(dir,n)
        if not os.path.exists(d):
            os.makedirs(d)
        else:
            continue           
            
def write_align_makeflow(samples,genome_dir,makeflow,input_d,output_d,thread,strand):
    gff2gtf = RNA.get_bin_abspath("gff2gtf.for.NCBIgff.pl")
    build_index = RNA.get_soft_abspath("bowtie2-build")
    bowtie = RNA.get_soft_abspath("bowtie2")
    reflat = RNA.get_bin_abspath("create_refFlat.pl")
    if not os.path.exists(genome_dir):
        print "%s dir not exists, please check!" %genome_dir
        sys.exit(0)
    gtf = glob.glob(os.path.join(genome_dir,"*.gtf"))
    gff = glob.glob(os.path.join(genome_dir,"*.gff"))
    fa = glob.glob(os.path.join(genome_dir,"*.fa"))
    if gtf and fa:
        fa = fa[0]
        gtf = gtf[0]
    elif not gtf and gff and fa:
        fa = fa[0]
        gff = gff[0]
        gtf = os.path.splitext(gff) + ".gtf"
        fo.write('CATEGORY=gff2gtf\n')
        fo.write(os.path.splitext(gff) + ".gtf : " + gff + "\n")
        fo.write("\tperl %s %s %s &> %s\n\n"%(gff2gtf,gff,gtf,os.path.join(os.path.dirname(gff) + "gff2gtf.log")))
    else:
        print "Error: No gtf or gff file in %s dir"%genome_dir
        sys.exit(1)

    if not fa:
        print "please check the fa and gtf file in %s" %genome_dir
        sys.exists(1)
        
    genome_prefix = os.path.splitext(fa)[0] + '.bowtie'
    fo =open(makeflow,"a+")
    if not os.path.exists(genome_dir+"/bowtie_build-genome.log"):
        fo.write("CATEGORY=Bowtie2-build\n")
        fo.write("%s : %s\n"%(genome_dir+"/bowtie_build-genome.log",fa))
        fo.write("\t%s --threads %d %s %s &> %s\n\n" %(build_index,thread,fa,genome_prefix,genome_dir+"/bowtie_build-genome.log"))
    if not os.path.exists(gtf + ".refFlat"):
        fo.write("CATEGORY=create_refFlat\n")
        fo.write(gtf + ".refFlat : " + gtf + "\n")
        fo.write("\tperl %s %s %s &> %s\n\n"%(reflat,gtf,gtf + ".refFlat",genome_dir+"/refflat.err"))
    for sn in samples:
        r1 = os.path.join(input_d,sn,"R1.clean.fq.gz")
        r2 = os.path.join(input_d,sn,"R2.clean.fq.gz")
        od = os.path.join(output_d,sn)
        fo.write("CATEGORY=Bowtie2-align_%s\n"%sn)
        fo.write(os.path.join(od,"bowtie.err ") + os.path.join(od,"accepted_hits.sam : ") + r1 + " " + r2 + " " + genome_dir + "/bowtie_build-genome.log\n")
        fo.write("@BATCH_OPTIONS = -l h_vmem=50G\n")
        if strand:
            fo.write("\t%s --ff -p %d -x %s -1 %s -2 %s -S %s "%(bowtie,thread,genome_prefix,r1,r2,os.path.join(od,"accepted_hits.sam")))
        else:
            fo.write("\t%s --fr -p %d -x %s -1 %s -2 %s -S %s "%(bowtie,thread,genome_prefix,r1,r2,os.path.join(od,"accepted_hits.sam")))
        fo.write("> %s 2> %s\n\n" %(os.path.join(od,"bowtie.out"),os.path.join(od,"bowtie.err")))
    fo.close()

def main():
    input_dir = os.path.abspath(args.input_dir)
    output_dir = os.path.abspath(args.output_dir)
    genome_dir = os.path.abspath(args.genome_dir)
    samples = sorted(map(os.path.basename,filter(os.path.isdir,[os.path.join(input_dir,x) for x in os.listdir(input_dir)])))
    mkdir(output_dir,samples)
    write_align_makeflow(samples,genome_dir,args.makeflow,input_dir,output_dir,args.parallel,args.strand)

if __name__=="__main__":
    main()
    
