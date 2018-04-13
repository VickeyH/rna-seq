#!/usr/bin/env python
#coding:utf-8
import argparse,os,re,sys,glob

parser = argparse.ArgumentParser(description="This Script is used to annotate and  enrichment go and kegg for gene list.")
parser.add_argument("-gtf","--gtf",type=str,help="The gtf file",required = True)
parser.add_argument("-updown","--updown",type=str,help="The input deg.up_down.list file. if more than two file, all files must be separate by comma without any space.",required = True)
parser.add_argument("-r","--ref_fasta",type=str,help="the reference fasta file",required = True)
parser.add_argument("-o","--output_dir",type=str,help="The output directory",required=True)
parser.add_argument("-v","--version",action="version",version='%(prog)s 1.0')
parser.add_argument("-vs","--vs",type=str,help="the diff groups, like CasegroupA-VS-ControlgroupB,CasegroupC-VS-ControlgroupD..., '-' symbol must not in group name, the order must be same as '--updown'",required = True)
parser.add_argument("-d","---dbtype",type=str,choices=["animals","bacteria","fungi","plants","protists"],help="the dbtype for blast to kobas and uniprot , you can see it in '/lustre/work/yongdeng/software/protokaryon/soft/gokegg/pep_fa' or '/lustre/work/yongdeng/software/protokaryon/soft/kobas/kobas_database'directory.",required = True)
parser.add_argument("-mf","--makeflow",type=str,help="the out makeflow file",required = True)
parser.add_argument("-fdr","--fdr",type=str,help="the method of adjust p_value,'BH' for Benjamini and Hochberg, 'BY' for Benjamini and Yekutieli, default fdr='BH'.",choices=["BH","BY"],default = "BH")
args = parser.parse_args()

def mkdir(dir,sampls):
    for n in sampls:
        d = os.path.join(dir,n)
        if not os.path.exists(d):
            os.makedirs(d)
        else:
            continue
            
def cp_genelist_cmd(genes_up_down,subdir):
    l = []
    for s,d in zip(genes_up_down,subdir):
        l.append("cp " + s + " " + d + "/deg.up_down.list")
    return " && ".join(l)
            
def write_preapre_blast_makeflow(makeflow,fa,gtf,outdir):
    fo = open(makeflow,"a+")
    fo.write("CATEGORY=extract_fa\n")
    fo.write(os.path.join(outdir,"blast","allgeneseq.fa : ")+gtf + " " + fa + "\n")
    fo.write("\tperl /lustre/work/yongdeng/software/protokaryon/soft/gokegg/extract_geneseq_bylongesttrans_inGTF.pl -Gtf %s -refFa %s -od %s &> %s\n\n"%(gtf,fa,os.path.join(outdir,"blast"),os.path.join(outdir,"blast/extract_geneseq.log")))
    fo.write("CATEGORY=split_fasta\n")
    fo.write(" ".join([os.path.join(outdir,"blast","allgeneseq.fa.") + str(i) for i in range(1,11)]) + " : " + os.path.join(outdir,"blast","allgeneseq.fa\n"))
    fo.write("\tpython /lustre/work/yongdeng/software/protokaryon/soft/gokegg/split_fastafile.py -i %s -n 10 -o %s &> %s\n\n"%(os.path.join(outdir,"blast","allgeneseq.fa"),os.path.join(outdir,"blast"),os.path.join(outdir,"blast/split_fa.log")))
    fo.close()
    
def write_blastx_makeflow(makeflow,fa_list,dbtype):
    fo = open(makeflow,"a+")
    for fa in fa_list:
        fo.write("CATEGORY=blastx_%s_to_kobas\n"%os.path.basename(fa))
        fo.write(fa+".kobas.blast : " + fa + "\n@BATCH_OPTIONS = -l h_vmem=10G\n")
        fo.write("\tblastx -query %s -db /lustre/work/yongdeng/software/protokaryon/soft/kobas/kobas_database/%s.fa -evalue 1e-5 -outfmt 6 -num_threads 8 -max_target_seqs 10 -out %s "%(fa,dbtype,fa+".kobas.blast"))
        fo.write("> %s 2> %s\n\n"%(fa + ".kobas.blast.out",fa + ".kobas.blast.err"))
        fo.write("CATEGORY=blastx_%s_to_uniprot\n"%os.path.basename(fa))
        fo.write(fa+".uniprot.blast : " + fa + "\n@BATCH_OPTIONS = -l h_vmem=10G\n")
        fo.write("\tblastx -query %s -db /lustre/work/yongdeng/software/protokaryon/soft/gokegg/pep_fa/uniprot.%s.pep.fa -evalue 1e-5 -outfmt 6 -num_threads 8 -max_target_seqs 10 -out %s "%(fa,dbtype,fa+".uniprot.blast"))
        fo.write("> %s 2> %s\n\n"%(fa + ".uniprot.blast.out",fa + ".uniprot.blast.err"))
    fo.write("CATEGORY=cat_kobas_blast\n")
    fo.write(os.path.join(os.path.dirname(fa),"kobas.blastx.outfmt6 : ") + " ".join([i + ".kobas.blast" for i in fa_list]) +"\n")
    fo.write("\tcat " + " ".join([i + ".kobas.blast" for i in fa_list]) + " >> " +  os.path.join(os.path.dirname(fa),"kobas.blastx.outfmt6\n\n"))
    fo.write("CATEGORY=cat_uniprot_blast\n")
    fo.write(os.path.join(os.path.dirname(fa),"uniprot.blastx.outfmt6 : ") + " ".join([i + ".uniprot.blast" for i in fa_list]) +"\n")
    fo.write("\tcat " + " ".join([i + ".uniprot.blast" for i in fa_list]) + " >> " +  os.path.join(os.path.dirname(fa),"uniprot.blastx.outfmt6\n\n"))
    fo.close()
    
def write_kegg_anno_makeflow(makeflow,annolist,subdir,outko):  ## annolist = kobas.blastx.outfmt6,.... 
    fo = open(makeflow,"a+")
    out = os.path.join(subdir,outko)   ### 输出文件路径和名字
    fo.write("CATEGORY=annotate_%s\n"%os.path.basename(subdir))
    fo.write( out + " : " + annolist + "\n")
    fo.write("\tpython /lustre/software/target/kobas-2.0/src/annotate.py -i %s -t blastout:tab -s ko -o %s "%(annolist,out))
    fo.write("> %s 2> %s\n\n"%(os.path.join(out +".annotate.out"),os.path.join(out +".annotate.err")))
    fo.write("CATEGORY=annotate_parse_%s\n"%os.path.basename(subdir))
    fo.write(" ".join([os.path.join(subdir) + "/KEGG_PATHWAY" + i for i in [".first.class",".second.class",".third.class",".first.list",".second.list",".third.list : "]]) + out + "\n")
    fo.write("\tperl /lustre/work/yongdeng/software/protokaryon/soft/gokegg/KobasAnnoForClassify.pl -KobasAnno %s -KEGGClassify %s -od %s "%(out,"/lustre/work/yongdeng/software/protokaryon/soft/gokegg/KEGG.pathway.classify.list",subdir))
    fo.write("&> %s\n\n"%os.path.join(subdir,"kegg.classify.log"))
    fo.write("CATEGORY=kegg_plot_%s\n"%os.path.basename(subdir))
    fo.write(os.path.join(subdir,"KEGG_PATHWAY.second.class.png ") + os.path.join(subdir,"KEGG_PATHWAY.second.class.pdf : ") + os.path.join(subdir,"KEGG_PATHWAY.second.class\n"))
    fo.write("\tRscript /lustre/work/yongdeng/software/protokaryon/soft/gokegg/kegg-plot.R %s %s %s "%(os.path.join(subdir,"KEGG_PATHWAY.second.class"),os.path.join(subdir,"KEGG_PATHWAY.second.class.png"),os.path.join(subdir,"KEGG_PATHWAY.second.class.pdf")))
    fo.write("&> %s\n\n"%os.path.join(subdir,"KEGG_PATHWAY.second.class.plot.log"))
    fo.close()
        
def write_kegg_ident_makeflow(makeflow,subdir,fdr="BH"):
    fo = open(makeflow,"a+")
    ko = os.path.join(subdir,"deg.ko")
    back_ko = os.path.join( os.path.dirname(subdir),"back","back.ko")
    fo.write("CATEGORY=kegg_identity_%s\n"%os.path.basename(subdir))
    fo.write(os.path.join(subdir,"deg.pathway.enrichment : ") + ko + " " + back_ko + "\n")
    fo.write("\tpython /lustre/software/target/kobas-2.0/bin/identify.py -c 1 -f %s -b %s -d K -o %s -n %s "%(ko,back_ko,os.path.join(subdir,"deg.pathway.enrichment"),fdr))
    fo.write("> %s 2> %s\n\n"%(os.path.join(subdir,"identify.out"),os.path.join(subdir,"identify.err")))
    fo.close()        
    
def write_ForGoKo_enrichmnet_makeflow(makeflow,annolist,subdir,output_dir):
    fo = open(makeflow,"a+")
    blast_dir = os.path.join(output_dir,"blast")
    for deg,dir in zip(annolist,subdir):
        dir_name = os.path.basename(dir)
        out = os.path.join(output_dir,dir)
        fo.write("CATEGORY=ForGoKo_enrichmnet_pre_%s\n"%dir_name)
        fo.write(os.path.join(out,"deg.uniprot.blastx ") + os.path.join(out,"deg.ko.blastx ") + os.path.join(out,"deg.go.list : ") + deg + " " + os.path.join(blast_dir,"uniprot.blastx.outfmt6 ") + os.path.join(blast_dir,"kobas.blastx.outfmt6\n"))
        fo.write("\tperl /lustre/work/yongdeng/software/protokaryon/soft/gokegg/ForGoKo_enrichmnet_pre.pl -bk %s -bn %s -od %s -deg %s -p deg &> %s\n\n"%(os.path.join(blast_dir,"kobas.blastx.outfmt6"),os.path.join(blast_dir,"uniprot.blastx.outfmt6"),out,deg,os.path.join(out,"deg.up_down.list.KoGo.log")))
    fo.close()

def write_go_enrichment_makeflow(makeflow,annolist,subdir,output_dir):
    fo = open(makeflow,"a+")
    blast_dir = os.path.join(output_dir,"blast")
    back_dir = os.path.join(output_dir,"back")
    fo.write("CATEGORY=prepare_go_back\n")
    fo.write(os.path.join(back_dir,"back.tsv : ") + os.path.join(back_dir,"back.ko ") + "/lustre/work/yongdeng/software/protokaryon/soft/gokegg/go-basic.obo\n")
    fo.write("\tperl /lustre/work/yongdeng/software/protokaryon/soft/gokegg/parse_kobas_annotate_blastx.pl %s %s %s "%(os.path.join(back_dir,"back.ko"),"/lustre/work/yongdeng/software/protokaryon/soft/gokegg/go-basic.obo",os.path.join(back_dir,"back")))
    fo.write("&> %s\n\n"%os.path.join(back_dir,"parse_kobas_annotate.log"))
    fo.write("CATEGORY=prepare_go_back\n")
    fo.write(os.path.join(back_dir,"back.go : ") + os.path.join(blast_dir,"uniprot.blastx.outfmt6 ") + os.path.join(back_dir,"back.tsv\n"))
    fo.write("\tperl /lustre/work/yongdeng/software/protokaryon/soft/gokegg/ForGoIdmapping_denovo3.pl -bn %s -od %s "%(os.path.join(blast_dir,"uniprot.blastx.outfmt6"),back_dir))
    fo.write("&> %s\n\n"%os.path.join(back_dir,"back.go.log"))
    fo.write("CATEGORY=back_go_idmapping\n")
    fo.write(os.path.join(back_dir,"back.ref.idmapping : ")+ os.path.join(back_dir,"back.go\n"))
    fo.write("\tperl /lustre/work/yongdeng/software/protokaryon/soft/gokegg/RNA-Seq_Anno1-3_id_mapping_byfile.pl -r -m /lustre/work/yongdeng/software/protokaryon/soft/gokegg/id:uniprot.dbfile -p back. -o %s %s "%(back_dir,os.path.join(back_dir,"back.go")))
    fo.write("&> %s\n\n"%os.path.join(back_dir,"RNA-Seq_Anno1-3_id_mapping_byfile.log"))
    fo.write("CATEGORY=back_go_idmapping\n")
    fo.write(os.path.join(back_dir,"back_go.idmapping : ")+ os.path.join(back_dir,"back.ref.idmapping\n"))
    fo.write("\tperl /lustre/work/yongdeng/software/protokaryon/soft/gokegg/idmapping_remove_ko.pl -in %s -od %s &> %s\n\n"%(os.path.join(back_dir,"back.ref.idmapping"),back_dir,os.path.join(back_dir,"idmapping_remove_ko.log")))
    fo.write("CATEGORY=back_go_normalization\n")
    fo.write(os.path.join(back_dir,"back.go.normalization_temp : ") + os.path.join(back_dir,"back_go.idmapping\n"))
    fo.write("\tperl /lustre/work/yongdeng/software/protokaryon/soft/gokegg/go_normalization_blast.pl /lustre/work/yongdeng/software/protokaryon/soft/gokegg/go-basic.obo %s %s &> %s\n\n"%(os.path.join(back_dir,"back_go.idmapping"),os.path.join(back_dir,"back.go.normalization_temp"),os.path.join(back_dir,"back.go_normalization_temp.log")))
    fo.write("CATEGORY=back_go_normalization\n")
    fo.write(os.path.join(back_dir,"back.go.normalization : ") + os.path.join(back_dir,"back.go.normalization_temp ") + os.path.join(blast_dir,"uniprot.blastx.outfmt6\n"))
    fo.write("\tperl /lustre/work/yongdeng/software/protokaryon/soft/gokegg/go_anno_denovo3.pl -nr %s -go %s -od %s &> %s && cp %s %s\n\n"%(os.path.join(blast_dir,"uniprot.blastx.outfmt6"),os.path.join(back_dir,"back.go.normalization_temp"),back_dir,os.path.join(back_dir,"go.final.idmapping.log"),os.path.join(back_dir,"go.final.idmapping"),os.path.join(back_dir,"back.go.normalization")))
    fo.write("CATEGORY=back_go_normalization\n")
    fo.write(os.path.join(back_dir,"go.annotation.list : ") + os.path.join(back_dir,"back.go.normalization\n"))
    fo.write("\tperl /lustre/work/yongdeng/software/protokaryon/soft/gokegg/obo_GO_anno.pl -obo /lustre/work/yongdeng/software/protokaryon/soft/gokegg/go-basic.obo -od %s -in %s &> %s\n\n"%(back_dir,os.path.join(back_dir,"back.go.normalization"),os.path.join(back_dir,"obo_GO_anno.log")))
    fo.write("CATEGORY=go_list_anno\n")
    fo.write(os.path.join(back_dir,"go.list : ") + os.path.join(back_dir,"go.annotation.list\n"))
    fo.write("\tperl /lustre/work/yongdeng/software/protokaryon/soft/gokegg/go_list_anno.pl -in %s -od %s &> %s\n\n"%(os.path.join(back_dir,"go.annotation.list"),back_dir,os.path.join(back_dir,"go.list")))
    fo.write("CATEGORY=back_anno\n")
    fo.write(os.path.join(back_dir,"back.anno.xls : ") + os.path.join(back_dir,"back.tsv ") + os.path.join(back_dir,"back.go.normalization ") + "/lustre/work/yongdeng/software/protokaryon/soft/gokegg/go-basic.obo\n")
    fo.write("\tperl /lustre/work/yongdeng/software/protokaryon/soft/gokegg/ForKobasTsv_blast.pl -tsv %s -od %s -database K/G -go %s -obo /lustre/work/yongdeng/software/protokaryon/soft/gokegg/go-basic.obo "%(os.path.join(back_dir,"back.tsv"),back_dir,os.path.join(back_dir,"back.go.normalization")))
    fo.write("&> %s\n\n"%os.path.join(back_dir,"back.anno.log"))   
    for deg,dir in zip(annolist,subdir):
        dir_name = os.path.basename(dir)
        out = os.path.join(output_dir,dir)
        fo.write("CATEGORY=go_enrichment_%s\n"%dir_name)
        fo.write(os.path.join(out,"deg.go.class ") + os.path.join(out,"deg.go.enrichment : ") + "/lustre/work/yongdeng/software/protokaryon/soft/gokegg/go-basic.obo " + deg + " " + os.path.join(back_dir,"back.go.normalization\n"))
        fo.write("\tperl /lustre/work/yongdeng/software/protokaryon/soft/gokegg/RNA-Seq_Anno3_enrichmentV1_blast.pl -g /lustre/work/yongdeng/software/protokaryon/soft/gokegg/go-basic.obo -o %s -q deg. %s %s "%(out,os.path.join(back_dir,"back.go.normalization"),deg))
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
    
def write_kobas_color_makeflow(makeflow,subdir):
    fo = open(makeflow,"a+")
    for vs in subdir:
        vs_name = os.path.basename(vs)
        dp_down = os.path.join(vs,"deg.up_down.list")
        for e in ["GO","Pathway"]:            
            fo.write("CATEGORY=kobas_color_%s_%s\n"%(e,vs_name))
            fo.write(os.path.join(vs,"%s.temp.xls : "%e) + os.path.join(vs,"enrichment.%s.txt\n"%e))
            fo.write("\tperl /lustre/work/yongdeng/software/protokaryon/soft/gokegg/kobas_colour_v2.1.pl -deg %s -k_anno %s -k_iden %s "%(dp_down,os.path.join(vs,"deg.ko"),os.path.join(vs,"enrichment.%s.txt"%e)))
            fo.write("-od %s -p %s.temp -gc 1 -ud 3 > %s 2> %s\n\n" %(vs,e,os.path.join(vs,"enrichment.%s.change.out"%e),os.path.join(vs,"enrichment.%s.change.err"%e)))
            fo.write("CATEGORY=enrichment_%s_xls_%s\n"%(e,vs_name))
            fo.write(os.path.join(vs,"enrichment.%s.xls : "%e) + os.path.join(vs,"%s.temp.xls\n"%e))
            fo.write("\tcut -f 1-11,14,15 %s > %s\n\n"%(os.path.join(vs,"%s.temp.xls"%e),os.path.join(vs,"enrichment.%s.xls"%e)))
    fo.close()

def main():
    outputdir = os.path.abspath(args.output_dir)
    makeflow = args.makeflow
    genes_up_down = [os.path.abspath(i) for i in args.updown.strip().split(",")]
    vs = args.vs.strip().split(",")
    subdir = [os.path.join(outputdir,v) for v in vs]
    fa = os.path.abspath(args.ref_fasta)
    gtf = os.path.abspath(args.gtf)
    mkdir(outputdir,["back","blast"] + vs)
    for i in genes_up_down:
      if not os.path.exists(i):
        print "%s file not exists,please check"%i
        sys.exit()
    os.system(cp_genelist_cmd(genes_up_down,subdir))
    write_preapre_blast_makeflow(makeflow,fa,gtf,outputdir)
    fa_list = [os.path.join(outputdir,"blast","allgeneseq.fa.") + str(i) for i in range(1,11)]
    write_blastx_makeflow(makeflow,fa_list,args.dbtype)
    write_ForGoKo_enrichmnet_makeflow(makeflow,genes_up_down,subdir,outputdir)
    write_kegg_anno_makeflow(makeflow,os.path.join(outputdir,"blast","kobas.blastx.outfmt6"),os.path.join(outputdir,"back"),"back.ko")
    for d in subdir:
        write_kegg_anno_makeflow(makeflow,os.path.join(d,"deg.ko.blastx"),d,"deg.ko")
    for x in subdir:
        write_kegg_ident_makeflow(makeflow,x,args.fdr)
    write_go_enrichment_makeflow(makeflow,genes_up_down,subdir,outputdir)
    write_go_pathway_summary_makeflow(makeflow,subdir)
    write_kobas_color_makeflow(makeflow,subdir)

if __name__ == "__main__"   :
    main()    
        
        
        
        
        
