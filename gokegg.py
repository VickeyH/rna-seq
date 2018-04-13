#!/usr/bin/env python
#coding:utf-8
## cd /lustre/work/yongdeng/Project/160761_test && python  ~/gokegg.py -gene /lustre/work/yongdeng/Project/161015_mmu/analysis/enrichment/ICA69-VS-WT/deg.list -updown /lustre/work/yongdeng/Project/161015_mmu/analysis/enrichment/ICA69-VS-WT/deg.up_down.list -vs ICA69-VS-WT -specise mmu -t id:ncbigene -d /K/R/B/p/G -mf gokegg_161015.mf -bg /lustre/work/yongdeng/Project/161015_mmu/analysis/enrichment/back.list -anno /lustre/work/zhonghuali/database/genomeInfo/Mus_musculus/38.84/fa.gtf/Mus_musculus.GRCm38.84.gtf.new.gtf.gene_id.list.ENSEMBL.anno -o enrichment

import argparse,os,re,sys
import RNA

parser = argparse.ArgumentParser(description="This Script is used to annotate and  enrichment go and kegg for gene list.")
parser.add_argument("-gene","--gene_list",type=str,help="The input deg.list file. If more than two file, all files must be separate by comma without any space, the order of file must be same as '--updown'",required = True)
parser.add_argument("-updown","--updown",type=str,help="The input deg.up_down.list file. if more than two file, all files must be separate by comma without any space, the order of file must be same as '--gene_list'",required = True)
parser.add_argument("-o","--output_dir",type=str,help="The output directory",required=True)
parser.add_argument("-v","--version",action="version",version='%(prog)s 1.0')
parser.add_argument("-vs","--vs",type=str,help="the diff groups, like CasegroupA-VS-ControlgroupB,CasegroupC-VS-ControlgroupD..., '-' symbol must not in group name, the order must be same as '--gene_list' and '--updown'",required = True)
parser.add_argument("-specise","---specise",type=str,help="the specise abbr.",required = True)
parser.add_argument("-t","---idtype",type=str,choices=["id:ncbigi","id:uniprot", "id:ensembl","id:ncbigene"],help="the input idtype of deg genes.",default = "id:ncbigene",required = False)
parser.add_argument("-d","---dbtype",type=str,help="the dbtype of kobas.",required = True)
parser.add_argument("-mf","--makeflow",type=str,help="the out makeflow file",required = True)
parser.add_argument("-bg","--backlist",type=str,help="the backgroud gene list",required = True)
parser.add_argument("-fdr","--fdr",type=str,help="the method of adjust p_value,'BH' for Benjamini and Hochberg, 'BY' for Benjamini and Yekutieli, 'BH' is default",choices=["BH","BY"],default = "BH")
parser.add_argument("-anno","--anno",type=str,required=True,help="The anno file of reference gtf")
args = parser.parse_args()

def mkdir(dir,sampls):
    for n in sampls:
        d = os.path.join(dir,n)
        if not os.path.exists(d):
            os.makedirs(d)
        else:
            continue
            
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
    fo.write(os.path.join(outputdir,"back","back.list ") + " ".join([os.path.join(outputdir,vs,i) for vs in subdir for i in ["deg.list","deg.up_down.list"]]) + " : " + os.path.abspath(back) + " " + " ".join(gene) + " " + " ".join(genes_up_down) + "\n")
    cp_list = ["cp " + back + " " + os.path.join(outputdir,"back","back.list")]
    for i in range(len(gene)):
        cp_list.append("cp " + gene[i] + " " + subdir[i])
    for i in range(len(genes_up_down)):
        cp_list.append("cp " + genes_up_down[i] + " " + subdir[i])
    fo.write("\t%s \n\n"%(" && ".join(cp_list)))
    fo.close()
    
def write_kegg_anno_makeflow(makeflow,annolist,subdir,idtype,abbr):
    annotate = RNA.get_soft_abspath("annotate")
    KobasAnnoForClassify = RNA.get_bin_abspath("KobasAnnoForClassify.pl")
    kegg_plot = RNA.get_bin_abspath("kegg-plot.R")
    fo = open(makeflow,"a+")
    out = os.path.join(subdir,os.path.splitext(os.path.basename(annolist))[0])
    fo.write("CATEGORY=annotate_%s\n"%os.path.basename(subdir))
    fo.write( out + ".ko : " + annolist + "\n")
    fo.write("\tpython %s -i %s -t %s -s %s -o %s "%(annotate,annolist,idtype,abbr,out + ".ko"))
    fo.write("> %s 2> %s\n\n"%(os.path.join(out +".annotate.out"),os.path.join(out +".annotate.err")))
    fo.write("CATEGORY=annotate_parse_%s\n"%os.path.basename(subdir))
    fo.write(" ".join([os.path.join(subdir) + "/KEGG_PATHWAY" + i for i in [".first.class",".second.class",".third.class",".first.list",".second.list",".third.list : "]]) + out + ".ko\n")
    fo.write("\tperl %s -KobasAnno %s -KEGGClassify %s -od %s "%(KobasAnnoForClassify,out+".ko","/lustre/work/yongdeng/software/protokaryon/soft/gokegg/KEGG.pathway.classify.list",subdir))
    fo.write("&> %s\n\n"%os.path.join(subdir,"kegg.classify.log"))
    fo.write("CATEGORY=kegg_plot_%s\n"%os.path.basename(subdir))
    fo.write(os.path.join(subdir,"KEGG_PATHWAY.second.class.png ") + os.path.join(subdir,"KEGG_PATHWAY.second.class.pdf : ") + os.path.join(subdir,"KEGG_PATHWAY.second.class\n"))
    fo.write("\tRscript %s %s %s %s "%(kegg_plot,os.path.join(subdir,"KEGG_PATHWAY.second.class"),os.path.join(subdir,"KEGG_PATHWAY.second.class.png"),os.path.join(subdir,"KEGG_PATHWAY.second.class.pdf")))
    fo.write("&> %s\n\n"%os.path.join(subdir,"KEGG_PATHWAY.second.class.plot.log"))
    fo.close()
    
def write_kegg_ident_makeflow(makeflow,subdir,dbtype,fdr="BH"):
    identify = RNA.get_soft_abspath("identify")
    dbtype = dbtype.replace("/G","")
    fo = open(makeflow,"a+")
    ko = os.path.join(subdir,"deg.ko")
    back_ko = os.path.join( os.path.dirname(subdir),"back","back.ko")
    fo.write("CATEGORY=kegg_identity_%s\n"%os.path.basename(subdir))
    fo.write(os.path.join(subdir,"deg.pathway.enrichment : ") + ko + " " + back_ko + "\n")
    fo.write("\tpython %s -c 1 -f %s -b %s -d %s -o %s -n %s "%(identify,ko,back_ko,dbtype,os.path.join(subdir,"deg.pathway.enrichment"),fdr))
    fo.write("> %s 2> %s\n\n"%(os.path.join(subdir,"identify.out"),os.path.join(subdir,"identify.err")))
    fo.close()
    
def write_go_enrichment_makeflow(makeflow,annolist,subdir,output_dir,dbtype,abbr,idtype):
    id_mapping = RNA.get_bin_abspath("RNA-Seq_Anno1-3_id_mapping_byfile.pl")
    idmapping_remove_ko = RNA.get_bin_abspath("idmapping_remove_ko.pl")
    go_normalization = RNA.get_bin_abspath("go_normalization.pl")
    parse_kobas_annotate = RNA.get_bin_abspath("parse_kobas_annotate.pl")
    ForKobasTsv = RNA.get_bin_abspath("ForKobasTsv.pl")
    RNA_Seq_Anno3_enrichmentV1 = RNA.get_bin_abspath("RNA-Seq_Anno3_enrichmentV1.pl")
    fo = open(makeflow,"a+")
    back_dir = os.path.join(output_dir,"back")
    go_info = os.popen("grep '^%s' /lustre/work/yongdeng/software/protokaryon/soft/gokegg/specise_v1.1 | awk '{print $NF}'"%abbr).read().strip()
    if go_info == "kobas_nogo" and idtype != "id:ncbigi":
        fo.write("CATEGORY=prepare_go\n")
        fo.write(os.path.join(back_dir,"back.ref.idmapping : ") + os.path.join(back_dir,"back.list\n"))
        fo.write("\tperl %s -r -m /lustre/work/yongdeng/software/protokaryon/soft/gokegg/%s.dbfile -p back. -o %s %s "%(id_mapping,idtype,back_dir,os.path.join(back_dir,"back.list")))
        fo.write("> %s 2> %s\n\n"%(os.path.join(back_dir,"RNA-Seq_Anno1-3_id_mapping_byfile.out"),os.path.join(back_dir,"RNA-Seq_Anno1-3_id_mapping_byfile.err")))
        fo.write("CATEGORY=prepare_go\n")
        fo.write(os.path.join(back_dir,"back_go.idmapping : ") + os.path.join(back_dir,"back.ref.idmapping\n"))
        fo.write("\tperl %s -in %s -od %s &> %s\n\n"%(idmapping_remove_ko,os.path.join(back_dir,"back.ref.idmapping"),back_dir,os.path.join(back_dir,"idmapping_remove_ko.log")))
        fo.write("CATEGORY=prepare_go\n")
        fo.write(os.path.join(back_dir,"back.go.normalization : ") + os.path.join(back_dir,"back_go.idmapping ") + "/lustre/work/yongdeng/software/protokaryon/soft/gokegg/go-basic.obo\n")
        fo.write("\tperl %s /lustre/work/yongdeng/software/protokaryon/soft/gokegg/go-basic.obo %s %s "%(go_normalization,os.path.join(back_dir,"back_go.idmapping"),os.path.join(back_dir,"back.go.normalization")))
        fo.write("> %s\n\n"%os.path.join(back_dir,"back.go_normalization.log"))
    elif go_info == "kobas_havego":    
        fo.write("CATEGORY=prepare_go\n")
        fo.write(os.path.join(back_dir,"back.go.normalization : ") + os.path.join(back_dir,"back.tsv ") + os.path.join(back_dir,"back.go ") + "/lustre/work/yongdeng/software/protokaryon/soft/gokegg/go-basic.obo\n")
        fo.write("\tperl %s /lustre/work/yongdeng/software/protokaryon/soft/gokegg/go-basic.obo %s %s "%(go_normalization,os.path.join(back_dir,"back.go"),os.path.join(back_dir,"back.go.normalization")))
        fo.write("&> %s\n\n"%os.path.join(back_dir,"back.go_normalization.log"))
    else:
        print "GO enrichment can not use %s type"%idtype
        sys.exit()
    fo.write("CATEGORY=prepare_go\n")
    fo.write(os.path.join(back_dir,"back.tsv ") + os.path.join(back_dir,"back.go : ") + os.path.join(back_dir,"back.ko ") + "/lustre/work/yongdeng/software/protokaryon/soft/gokegg/go-basic.obo\n")
    fo.write("\tperl %s %s %s %s "%(parse_kobas_annotate,os.path.join(back_dir,"back.ko"),"/lustre/work/yongdeng/software/protokaryon/soft/gokegg/go-basic.obo",os.path.join(back_dir,"back")))
    fo.write("&> %s\n\n"%os.path.join(back_dir,"parse_kobas_annotate.log"))        
    fo.write("CATEGORY=back_anno\n")
    fo.write(os.path.join(back_dir,"back.anno.xls : ") + os.path.join(back_dir,"back.tsv ") + os.path.join(back_dir,"back.go.normalization ") + "/lustre/work/yongdeng/software/protokaryon/soft/gokegg/go-basic.obo\n")
    fo.write("\tperl %s -tsv %s -od %s -database %s -go %s -obo /lustre/work/yongdeng/software/protokaryon/soft/gokegg/go-basic.obo "%(ForKobasTsv,os.path.join(back_dir,"back.tsv"),back_dir,dbtype,os.path.join(back_dir,"back.go.normalization")))
    fo.write("&> %s\n\n"%os.path.join(back_dir,"back.anno.log"))
    
    for deg,dir in zip(annolist,subdir):   ### annolist为各样本分组的deg.up_down.list文件，subdir为group1-VS-group2目录
        dir_name = os.path.basename(dir)
        out = os.path.join(output_dir,dir)
        fo.write("CATEGORY=go_enrichment_%s\n"%dir_name)
        fo.write(os.path.join(out,"deg.go.class ") + os.path.join(out,"deg.go.enrichment : ") + "/lustre/work/yongdeng/software/protokaryon/soft/gokegg/go-basic.obo " + deg + " " + os.path.join(back_dir,"back.go.normalization\n"))
        fo.write("\tperl %s -g /lustre/work/yongdeng/software/protokaryon/soft/gokegg/go-basic.obo -o %s -q deg. %s %s "%(RNA_Seq_Anno3_enrichmentV1,out,os.path.join(back_dir,"back.go.normalization"),deg))
        fo.write(" > %s 2> %s\n\n"%(os.path.join(out,"RNA-Seq_Anno3_enrichmentV1.out"),os.path.join(out,"RNA-Seq_Anno3_enrichmentV1.err")))
    fo.close()
    
def write_go_pathway_summary_makeflow(makeflow,subdir):
    ForKobasDot = RNA.get_bin_abspath("ForKobasDot_v1.1.pl")
    parse_kobas_identify = RNA.get_bin_abspath("parse_kobas_identify.pl")
    draw_barplot_pvalue = RNA.get_bin_abspath("draw_barplot_pvalue.r")
    draw_barplot = RNA.get_bin_abspath("draw_barplot.r")
    fo = open(makeflow,"a+")
    for vs in subdir:
        vs_name = os.path.basename(vs)
        go_enrich = os.path.join(vs,"deg.go.enrichment")
        path_enrich = os.path.join(vs,"deg.pathway.enrichment")
        fo.write("CATEGORY=ForKobasDot_%s\n"%vs_name)
        fo.write(os.path.join(vs,"enrichment.xls : ") + go_enrich + " " + path_enrich + "\n")
        fo.write("\tperl %s -ge %s -ke %s -o %s &> %s\n\n"%(ForKobasDot,go_enrich,path_enrich,vs,os.path.join(vs,"ForKobasDot.log")))
        fo.write("CATEGORY=parse_identify_%s\n"%vs_name)
        fo.write(" ".join([os.path.join(vs,"enrichment.") + i for i in ["GO.hierarchy.molecular_function.dot","GO.hierarchy.cellular_component.dot","GO.hierarchy.biological_process.dot","GO.txt","Pathway.txt","Disease.txt : "]]) + os.path.join(vs,"enrichment.xls\n"))
        fo.write("\tperl %s %s /lustre/work/yongdeng/software/protokaryon/soft/gokegg/go-basic.obo 30 %s "%(parse_kobas_identify,os.path.join(vs,"enrichment.xls"),os.path.join(vs,"enrichment")))
        fo.write("> %s 2> %s\n\n"%(os.path.join(vs,"parse_kobas_identify.out"),os.path.join(vs,"parse_kobas_identify.err")))
        for level in ["biological_process","molecular_function","cellular_component"]:
            for p in ["png","pdf"]:
                fo.write("CATEGORY=enrichment.GO_%s_%s\n"%(level,p))
                fo.write(os.path.join(vs,"enrichment.GO.hierarchy.") + level + ".%s : "%p + os.path.join(vs,"enrichment.GO.hierarchy.") + level + ".dot\n")
                fo.write("\tdot -T" + p + " -o " + os.path.join(vs,"enrichment.GO.hierarchy.") + level + ".%s "%p + os.path.join(vs,"enrichment.GO.hierarchy.") + level + ".dot &> %s\n\n"%os.path.join(vs,"enrichment.GO.hierarchy." + level + "." + p + ".log"))
        fo.write("CATEGORY=GO_pvalue_%s\n"%vs_name)
        fo.write(os.path.join(vs,"enrichment.GO.p_value.top.pdf ") + os.path.join(vs,"enrichment.GO.p_value.top.png : ") + os.path.join(vs,"enrichment.GO.txt\n"))
        fo.write("\tRscript %s GO 30 %s %s %s &> %s\n\n"%(draw_barplot_pvalue,os.path.join(vs,"enrichment.GO.txt"),os.path.join(vs,"enrichment.GO.p_value"),os.path.join(vs,"enrichment.GO.p_value"),os.path.join(vs,"enrichment.GO.p_value.log")))
        fo.write("CATEGORY=GO_qvalue_%s\n"%vs_name)
        fo.write(os.path.join(vs,"enrichment.GO.q_value.top.pdf ") + os.path.join(vs,"enrichment.GO.q_value.top.png : ") + os.path.join(vs,"enrichment.GO.txt\n"))
        fo.write("\tRscript %s GO 30 %s %s %s &> %s\n\n"%(draw_barplot,os.path.join(vs,"enrichment.GO.txt"),os.path.join(vs,"enrichment.GO.q_value"),os.path.join(vs,"enrichment.GO.q_value"),os.path.join(vs,"enrichment.GO.q_value.log")))
    fo.close()
    
def write_kobas_color_makeflow(makeflow,subdir,anno,idtype):
    kobas_colour = RNA.get_bin_abspath("kobas_colour1.6.pl")
    fo = open(makeflow,"a+")
    for vs in subdir:
        vs_name = os.path.basename(vs)
        dp_down = os.path.join(vs,"deg.up_down.list")
        for e in ["GO","Pathway"]:            
            fo.write("CATEGORY=kobas_color_%s_%s\n"%(e,vs_name))
            fo.write(os.path.join(vs,"%s.temp.xls : "%e) + os.path.join(vs,"enrichment.%s.txt\n"%e))
            fo.write("\tperl %s -deg %s -k_anno %s -k_iden %s -inf %s "%(kobas_colour,dp_down,os.path.join(vs,"deg.ko"),os.path.join(vs,"enrichment.%s.txt"%e),anno))
            fo.write("-S symbol -Y %s -od %s -p %s.temp -gc 1 -ud 3 > %s 2> %s\n\n" %(idtype,vs,e,os.path.join(vs,"enrichment.%s.change.out"%e),os.path.join(vs,"enrichment.%s.change.err"%e)))
            fo.write("CATEGORY=enrichment_%s_xls_%s\n"%(e,vs_name))
            fo.write(os.path.join(vs,"enrichment.%s.xls : "%e) + os.path.join(vs,"%s.temp.xls\n"%e))
            fo.write("\tcut -f 1-13 %s > %s\n\n"%(os.path.join(vs,"%s.temp.xls"%e),os.path.join(vs,"enrichment.%s.xls"%e)))
    fo.close()
            
def main():
    outputdir = os.path.abspath(args.output_dir)
    genes = [os.path.abspath(i) for i in args.gene_list.strip().split(",")]
    genes_up_down = [os.path.abspath(i) for i in args.updown.strip().split(",")]
    vs = args.vs.strip().split(",")
    subdir = [os.path.join(outputdir,v) for v in vs]
    mkdir(outputdir,["back"] + vs)
    # for i in genes + genes_up_down + [args.backlist]:
        # if not os.path.exists(i):
            # print "%s file not exists,please check"%i
            # sys.exit()
    # os.system("cp %s %s"%(args.backlist,os.path.join(outputdir,"back","back.list")))
    # os.system(cp_genelist_cmd(genes,genes_up_down,subdir))
    write_cp_gene_list_makeflow(args.makeflow,outputdir,os.path.abspath(args.backlist),genes,genes_up_down,subdir)
    deglist = [os.path.join(i,"deg.list") for i in subdir]
    write_kegg_anno_makeflow(args.makeflow,os.path.join(outputdir,"back","back.list"),os.path.join(outputdir,"back"),args.idtype,args.specise)
    for d in subdir:
        write_kegg_anno_makeflow(args.makeflow,os.path.join(d,"deg.list"),d,args.idtype,args.specise)
    for x in subdir:
        write_kegg_ident_makeflow(args.makeflow,x,args.dbtype,args.fdr)        
    write_go_enrichment_makeflow(args.makeflow,genes_up_down,subdir,outputdir,args.dbtype,args.specise,args.idtype)
    write_go_pathway_summary_makeflow(args.makeflow,subdir)    
    write_kobas_color_makeflow(args.makeflow,subdir,args.anno,args.idtype)
    
if __name__ == "__main__"   :
    main()
    
    
### 多个分组的deg.list和deg.up_down.list, 对应多个分组的vs参数, 所有的分组共用一个back.list, 该back.list在输出目录的子目录之下    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
