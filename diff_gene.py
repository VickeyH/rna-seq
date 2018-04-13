#!/usr/bin/env python
#coding:utf-8
## cd /lustre/work/yongdeng/Project/160761_test/deseq2 && python ~/diff_gene.py -c gene_count_matrix.txt -anno /lustre/work/yongdeng/GenomeRef/Lactobacillus_acidophilus/anno4.txt -fc 2 -s 2641,2642,2643,2720IB,2821,3136IB,3328NT,3404IB,3551 -g IR,IR,NT,IB,IR,IB,NT,IB,NT -o deg -vs IB-VS-NT,IR-VS-NT -t id:ncbigene -r ensembl -mf diff_gene.mf
import argparse,os,sys

parser = argparse.ArgumentParser(description="This Script is used to do DEG analysis")
parser.add_argument("-c","--count_tab",type=str,help="The read-count table",required = True)
parser.add_argument("-fc","--fold_change",type=float,help="The fc Threshold for diff_gene, 2 is default",default = 2)
parser.add_argument("-anno","--anno",type=str,help="the anno file of ref species",required = True)
parser.add_argument("-s","--sample",type=str,help="The sample name,the order must be same as group name",required = True,nargs="*")
parser.add_argument("-g","--group",type=str,help="The group name,the order must be same as sample name",required = True,nargs="*")
parser.add_argument("-o","--output_dir",type=str,help="The output directory",required = True)
parser.add_argument("-vs","--vs",type=str,help="the diff groups, like CasegroupA-VS-ControlgroupB ...", required = True,nargs="*")
parser.add_argument("-t","--dbtype",type=str,help="the type of gene_id used for kobsa, default=id:ncbigene",default = "id:ncbigene",choices=["id:ensembl","id:ncbigene","id:uniprot"])
parser.add_argument("-r","--reftype",type=str,help="the type of gene_id in your diff file(first colnum), default=ensembl",default="ensembl",choices=["ensembl","entrez","uniprot"])
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
            
def write_deg_analysis_makeflow(makeflow,countable,samples,groups,vs_list,outdir):
    fo = open(makeflow,"a+")
    fo.write("CATEGORY=get_deg_rscript\n")
    fo.write(os.path.join(outdir,"diff_gene_DESeq2.r : ") + countable + "\n")
    fo.write("\tpython /lustre/work/yongdeng/software/protokaryon/flow/get_deg2.py -c %s -s %s -g %s -vs %s -o %s "%(countable,",".join(samples),",".join(groups),",".join(vs_list),outdir))
    fo.write('> %s 2> %s\n\n'%(os.path.join(outdir,"get_deg_rscript.out"),os.path.join(outdir,"get_deg_rscript.err")))
    fo.write("CATEGORY=deg_analysis\n")
    #fo.write(os.path.join(outdir,"gene.rpm.count.table ")  + os.path.join(outdir,"gene.rlog2.counts.table ") + " ".join([os.path.join(outdir,i,i+".txt") for i in vs_list]) + " : " + os.path.join(outdir,"diff_gene_DESeq2.r ") + countable + "\n")
    fo.write(os.path.join(outdir,"gene_count_matrix.normalized.txt ") + " ".join([os.path.join(outdir,i,i+".txt") for i in vs_list]) + " : " + os.path.join(outdir,"diff_gene_DESeq2.r ") + countable + "\n")
    fo.write("\tRscript %s > %s 2> %s\n\n"%(os.path.join(outdir,"diff_gene_DESeq2.r"),os.path.join(outdir,"diff_gene_DESeq2.out"),os.path.join(outdir,"diff_gene_DESeq2.err")))
    fo.close()
    
def write_stat_diff_gene_makeflow(makeflow,vs_list,outdir):
    fo = open(makeflow,"a+")
    fo.write("CATEGORY=get_scatter_valcon\n")
    fo.write(" ".join([os.path.join(outdir,i,i+j) for i in vs_list for j in [".diff_info",".down.diff",".plot_scatter.file",".plot_valcano.file",".up.diff"]]) + " : " + " ".join([os.path.join(outdir,i,i+".txt") for i in vs_list]) + "\n")
    fo.write("\tpython /lustre/work/yongdeng/software/protokaryon/flow/pick_diffexp.py %s &> %s\n\n"%(" ".join([os.path.join(outdir,i,i+".txt") for i in vs_list]),os.path.join(outdir,"pick_diffexp.log")))
    for vs in vs_list:
        fo.write("CATEGORY=get_scatter_valcon_%s\n"%vs)
        fo.write(os.path.join(outdir,vs,vs+".plot_scatter.png : ") + os.path.join(outdir,vs,vs+".plot_scatter.file\n"))
        fo.write("\tRscript /lustre/work/yongdeng/software/protokaryon/soft/deseq/Scatter_volcano.r %s %s xy &> %s\n\n"%(os.path.join(outdir,vs,vs+".plot_scatter.file"),os.path.join(outdir,vs,vs+".plot_scatter.png"),os.path.join(outdir,vs,"Scatter_rscript.log")))
        fo.write("CATEGORY=get_scatter_valcon_%s\n"%vs)
        fo.write(os.path.join(outdir,vs,vs+".plot_valcano.png : ") + os.path.join(outdir,vs,vs+".plot_valcano.file\n"))
        fo.write("\tRscript /lustre/work/yongdeng/software/protokaryon/soft/deseq/Scatter_volcano.r %s %s &> %s\n\n"%(os.path.join(outdir,vs,vs+".plot_valcano.file"),os.path.join(outdir,vs,vs+".plot_valcano.png"),os.path.join(outdir,vs,"Volcano_rscript.log")))
    fo.write("CATEGORY=stat_diff\n")
    fo.write(" ".join([os.path.join(outdir,vs,vs+".diff_stat") for vs in vs_list]) + " : " + " ".join([os.path.join(outdir,i,i+".txt") for i in vs_list]) + "\n")
    fo.write("\tpython /lustre/work/yongdeng/software/protokaryon/flow/stat_diffexp.py %s &> %s\n\n"%(" ".join([os.path.join(outdir,i,i+".txt") for i in vs_list]),os.path.join(outdir,"stat_diffexp.log")))     
    fo.close()    

def write_back_deg_list_makeflow(makeflow,vs_list,anno,outdir,dbtype,reftype):
    fo = open(makeflow,"a+")
    fo.write("CATEGORY=get_deg_back_list\n")
    fo.write(os.path.join(outdir,"back.list ") + " ".join([os.path.join(outdir,vs,j) for vs in vs_list for j in ["deg.list","deg.num.txt","deg.up.txt","deg.down.txt","deg.up_down.list"]]) + " : " + anno + " " + " ".join([os.path.join(outdir,i,i+".txt") for i in vs_list]) + "\n")
    fo.write("\tpython /lustre/work/yongdeng/software/protokaryon/flow/call_deg_back_list.py -anno %s -t %s -r %s -diff %s &> %s\n\n"%(anno,dbtype,reftype," ".join([os.path.join(outdir,i,i+".txt") for i in vs_list]),os.path.join(outdir,"call_deg_back_list.log")))
    fo.close()
    
def write_stat_deg_num_makeflow(makeflow,outdir,vs_list):
    fo = open(makeflow,"a+")
    fo.write("CATEGORY=deg_num_stat\n")
    fo.write(os.path.join(outdir,"deg.state.txt : ") + " ".join([os.path.join(outdir,vs,"deg.num.txt") for vs in vs_list]) + "\n")
    fo.write("\tpython /lustre/work/yongdeng/software/protokaryon/flow/stat_deg_num.py %s %s &> %s\n\n"%(" ".join([os.path.join(outdir,vs,"deg.num.txt") for vs in vs_list]),os.path.join(outdir,"deg.state.txt"),os.path.join(outdir,"stat_deg_num.log")))
    fo.write("CATEGORY=deg_stat_plot\n")
    fo.write(os.path.join(outdir,"deg.stat.png ") + os.path.join(outdir,"deg.stat.pdf : ") + os.path.join(outdir,"deg.state.txt\n"))
    fo.write("\tRscript /lustre/work/yongdeng/software/protokaryon/soft/deseq/deg_stat.r %s %s %s &> %s\n\n"%(os.path.join(outdir,"deg.state.txt"),os.path.join(outdir,"deg.stat.png"),os.path.join(outdir,"deg.stat.pdf"),os.path.join(outdir,"deg_stat.log")))
    fo.close()
    
if __name__ == "__main__":
    makeflow = args.makeflow
    outdir = os.path.abspath(args.output_dir)
    countable = os.path.abspath(args.count_tab)
    anno = os.path.abspath(args.anno)
    mkdir(outdir,args.vs)
    write_deg_analysis_makeflow(makeflow,countable,args.sample,args.group,args.vs,outdir)
    write_stat_diff_gene_makeflow(makeflow,args.vs,outdir)
    write_back_deg_list_makeflow(makeflow,args.vs,anno,outdir,args.dbtype,args.reftype)
    write_stat_deg_num_makeflow(makeflow,outdir,args.vs)
    
    
    
