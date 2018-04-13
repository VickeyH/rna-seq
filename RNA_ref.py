#!/usr/bin/env python
#coding:utf-8
import argparse,os,sys,datetime,commands,glob,re
from ConfigParser import ConfigParser

def mkdir(dir,sampls):
    for n in sampls:
        d = os.path.join(dir,n)
        if not os.path.exists(d):
            os.makedirs(d)
        else:
            continue
            
def get_sample_dict(sampleinfo):
    sampleDict = {}
    for i in re.findall('(\(.+?\))',sampleinfo):
        l = i.strip("()").split()
        for s in l[2:]:
            sampleDict[s] = l[0]
    return sampleDict
            
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
           
parser = argparse.ArgumentParser(description="This Script is used to creat a makeflow for RNA_ref.")
parser.add_argument("-config","--config",type=str,help="The config file of this project",required = True)
parser.add_argument("-o","--output_dir",type=str,help="the output dir of all output",required = True)
parser.add_argument("-mf","--makeflow",type=str,help="the output makeflow file name",required = True)
args = parser.parse_args()

config = ConfigParser(defaults={"threads":"8","fdr":"BH","Strand":"no","fold_change":"2","outdir":os.path.abspath(args.output_dir),"method":"trimmomatic"})
config.read(args.config)
         
flow_name = "flow/Project_" + config.get("global","project_name") + "_workflow_" + datetime.datetime.now().strftime("%Y_%m_%d_%H_%M_%S")
od =  config.get("DEFAULT","outdir")
makeflow = os.path.join(od,args.makeflow)   
if os.path.exists(makeflow):
    print "%s makeflow file exists, Please change another name or remove this file"%makeflow
    sys.exit()  
mkdir(od,[flow_name,"analysis"])
sh = open(os.path.join(od,flow_name,"flow.sh"),"w")
log = open(os.path.join(od,flow_name,"creat_all_flow.log"),"w")

modules = config.get("global","moduels").split()
thread = config.getint("global","threads") if config.get("global","threads") else config.getint("global","threads")
sampleinfo = config.get("global","sample_info")
sd = get_sample_dict(sampleinfo)

if "rawdata" in modules:
    inputdir = os.path.abspath(config.get("rawdata","raw_dir"))
    cat_raw_output_dir = os.path.join(od,"analysis","rawdata")
    samples = ",".join(sorted(config.get("rawdata","sample_name").split()))
    mf = os.path.join(od,flow_name,"cat_rawdata.mf")
    cmd = 'python /lustre/work/yongdeng/software/protokaryon/flow/cat_rawdata.py -i %s -s %s -mf %s -o %s' %(inputdir,samples,mf,cat_raw_output_dir)
    sh.write("\t"+cmd+"\n\n")
    status,logout = commands.getstatusoutput(cmd);log.write(logout)
    os.popen('cat %s >> %s'%(mf,makeflow))
    if status != 0:
        print "Error in cat_rawdata.py, please check the %s file" %os.path.join(od,flow_name,+"creat_all_flow.log")
        sh.close();log.close();sys.exit(1)
    
if "QC" in modules:
    inputdir = os.path.join(od,"analysis","rawdata")
    qc_outputdir = os.path.join(od,"analysis","cleandata")
    mf = os.path.join(od,flow_name,"qc_rawdata.mf")
    methods = config.get("QC","method") if config.get("QC","method") else config.get("DEFAULT","method") 
    cmd = 'python /lustre/work/yongdeng/software/protokaryon/flow/QC_rawdata.py -p %d -i %s -mf %s -o %s -m %s' %(thread,inputdir,mf,qc_outputdir,methods)
    sh.write("\t"+cmd+"\n\n")
    status,logout = commands.getstatusoutput(cmd);log.write(logout)
    os.popen('cat %s >> %s'%(mf,makeflow))
    if status != 0:
        print "Error in QC_rawdata.py, please check the %s file" %os.path.join(od,flow_name,"creat_all_flow.log")
        sh.close();log.close();sys.exit(1)

if "align" in modules:
    inputdir = os.path.join(od,"analysis","cleandata")
    genome_dir = config.get("align","genome_dir")
    align_output_dir = os.path.join(od,"analysis","align")
    mf = os.path.join(od,flow_name,"align_clean.mf")
    cmd = 'python /lustre/work/yongdeng/software/protokaryon/flow/align_clean.py -p %d -i %s -mf %s -o %s -g %s' %(thread,inputdir,mf,align_output_dir,genome_dir)
    sh.write("\t"+cmd+"\n\n")
    status,logout = commands.getstatusoutput(cmd);log.write(logout)
    os.popen('cat %s >> %s'%(mf,makeflow))
    if status != 0:
        print "Error in align_clean.py, please check the %s file" %os.path.join(od,flow_name,"creat_all_flow.log")
        sh.close();log.close();sys.exit(1)
     
if "summary_data" in modules:
    RawDir = os.path.join(od,"analysis","rawdata")
    CleanDir = os.path.join(od,"analysis","cleandata")
    genome_dir = config.get("align","genome_dir")
    methods = config.get("QC","method") if config.get("QC","method") else config.get("DEFAULT","method")
    reflat = glob.glob(os.path.join(genome_dir,"*.gtf"))[0] + ".refFlat"
    strand = config.get("DEFAULT","Strand") if not config.get("summary_data","Strand") else config.get("summary_data","Strand")
    mf = os.path.join(od,flow_name,"summary_data.mf")
    cmd = 'python /lustre/work/yongdeng/software/protokaryon/flow/summary_data.py -p %d -raw %s -clean %s -reflat %s -m %s -mf %s' %(thread,RawDir,CleanDir,reflat,methods,mf)
    sh.write("\t"+cmd+"\n\n")
    status,logout = commands.getstatusoutput(cmd);log.write(logout)
    os.popen('cat %s >> %s'%(mf,makeflow))
    if status != 0:
        print "Error in summary_data.py, please check the %s file" %os.path.join(od,flow_name,"creat_all_flow.log")
        sh.close();log.close();sys.exit(1)

if "read_count" in modules:
    strand = config.get("DEFAULT","Strand") if not config.get("read_count","Strand") else config.get("read_count","Strand")
    genome_dir = config.get("align","genome_dir")
    gtf = glob.glob(os.path.join(genome_dir,"*.gtf"))[0]
    inputdir = os.path.join(od,"analysis","align")
    read_count_output_dir = os.path.join(od,"analysis","Htseq")
    mf = os.path.join(od,flow_name,"Htseq_count.mf")
    cmd = 'python /lustre/work/yongdeng/software/protokaryon/flow/Htseq_count.py -p %d -i %s -g %s -s %s -o %s -mf %s'%(thread,inputdir,gtf,strand,read_count_output_dir,mf)
    sh.write("\t"+cmd+"\n\n")
    status,logout = commands.getstatusoutput(cmd);log.write(logout)
    os.popen('cat %s >> %s'%(mf,makeflow))
    if status != 0:
        print "Error in Htseq_count.py, please check the %s file" %os.path.join(od,flow_name,"creat_all_flow.log")
        sh.close();log.close();sys.exit(1)

if "diff_gene" in modules:
    countable = os.path.join(od,"analysis","Htseq","gene_count_matrix.txt")
    anno = config.get("diff_gene","anno")
    fc = config.get("DEFAULT","fold_change") if not config.get("diff_gene","fold_change") else config.get("diff_gene","fold_change")
    samples = ",".join([i[0] for i in sorted(sd.items(),key=lambda x:x[0])])
    group = ",".join([i[1] for i in sorted(sd.items(),key=lambda x:x[0])])
    diff_out_dir = os.path.join(od,"analysis","DESeq2")
    vs = ",".join(config.get("diff_gene","vs").split())
    dbtype = config.get("diff_gene","dbtype")
    reftype = config.get("diff_gene","reftype")
    mf = os.path.join(od,flow_name,"diff_gene.mf")
    cmd = 'python /lustre/work/yongdeng/software/protokaryon/flow/diff_gene.py -c %s -anno %s -fc %s -s %s -g %s -o %s -vs %s -mf %s -t %s -r %s'%(countable,anno,fc,samples,group,diff_out_dir,vs,mf,dbtype,reftype)
    sh.write("\t"+cmd+"\n\n")
    status,logout = commands.getstatusoutput(cmd);log.write(logout)
    os.popen('cat %s >> %s'%(mf,makeflow))
    if status != 0:
        print "Error in diff_gene.py, please check the %s file" %os.path.join(od,flow_name,"creat_all_flow.log")
        sh.close();log.close();sys.exit(1)
        
if "go_kegg" in modules and "go_kegg_blast" not in modules:
    fdr = config.get("go_kegg","fdr") if config.get("go_kegg","fdr") else config.get("DEFAULT","fdr")
    bg = os.path.join(od,"analysis","DESeq2","back.list")
    vs_list = config.get("diff_gene","vs").split()
    gene = ",".join([os.path.join(od,"analysis","DESeq2",i,"deg.list") for i in vs_list])
    updown = ",".join([os.path.join(od,"analysis","DESeq2",i,"deg.up_down.list") for i in vs_list])
    vs = ",".join(vs_list)
    specise = config.get("go_kegg","specise_abbr")
    dbtype = config.get("diff_gene","dbtype") if not config.get("go_kegg","dbtype") else config.get("go_kegg","dbtype")
    kobasdb = config.get("go_kegg","kobasdb")
    anno  = config.get("diff_gene","anno")
    gokegg_output_dir = os.path.join(od,"analysis","enrichment")
    mf = os.path.join(od,flow_name,"go_kegg.mf")
    cmd = 'python /lustre/work/yongdeng/software/protokaryon/flow/gokegg.py -fdr %s -gene %s -updown %s -vs %s -specise %s -t %s -d %s -mf %s -bg %s -anno %s -o %s'%(fdr,gene,updown,vs,specise,dbtype,kobasdb,mf,bg,anno,gokegg_output_dir)
    sh.write("\t"+cmd+"\n\n")
    status,logout = commands.getstatusoutput(cmd);log.write(logout)
    os.popen('cat %s >> %s'%(mf,makeflow))
    if status != 0:
        print "Error in go_kegg.py, please check the %s file" %os.path.join(od,flow_name,"creat_all_flow.log")
        sh.close();log.close();sys.exit(1)
 
if "go_kegg_blast" in modules and "go_kegg" not in modules:
    fdr = config.get("go_kegg_blast","fdr") if config.get("go_kegg_blast","fdr") else config.get("DEFAULT","fdr")
    genome_dir = config.get("align","genome_dir")
    gtf = glob.glob(os.path.join(genome_dir,"*.gtf"))[0]
    vs_list = config.get("diff_gene","vs").split()
    updown = ",".join([os.path.join(od,"analysis","DESeq2",i,"deg.up_down.list") for i in vs_list])
    vs = ",".join(vs_list)
    dbtype = config.get("go_kegg_blast","dbtype")
    ref_fasta = glob.glob(os.path.join(genome_dir,"*.fa"))[0]
    gokegg_blast_output_dir = os.path.join(od,"analysis","enrichment")
    mf = os.path.join(od,flow_name,"go_kegg_blast.mf")
    cmd = 'python /lustre/work/yongdeng/software/protokaryon/flow/gokegg_blast.py -gtf %s -updown %s -r %s -o %s -vs %s -d %s -mf %s -fdr %s'%(gtf,update,ref_fasta,gokegg_blast_output_dir,vs,dbtype,mf,fdr)
    sh.write("\t"+cmd+"\n\n")
    status,logout = commands.getstatusoutput(cmd);log.write(logout)
    os.popen('cat %s >> %s'%(mf,makeflow))
    if status != 0:
        print "Error in go_kegg_blast.py, please check the %s file" %os.path.join(od,flow_name,"creat_all_flow.log")
        sh.close();log.close();sys.exit(1)

bool_value,message = check_makeflow(makeflow)
sys.stdout.write(message)
log.write(message)
sh.close()
log.close()



