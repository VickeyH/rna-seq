#!/usr/bin/env python
#coding:utf-8

#python script.py gtf abbr uniport.tab anno

### uniprot_id先从kobsa库里面查找，没有的再从uniport文件中按照symbol查找
import sys,os,re

def get_gtf_gen_group(gtfile):
    gtf_group = []
    with open(gtfile) as gtf:
        a = []
        for line in gtf:
            if len(a) == 0:
                a.append(line)
                continue
            if "\tgene\t" in line:
                gtf_group.append(a)
                a=[line,]
            else:
                a.append(line)
        gtf_group.append(a)
    return gtf_group
    
def get_uniport_dict(abbr):
    d = {}
    with open("/lustre/work/yongdeng/software/protokaryon/soft/kobas/kobas_database/sqlit3/"+abbr + "/GeneUniprotkbAcs.txt") as udb:
        for line in udb:
            k = line.split(":")[-1].split()[0]
            v = line.split(":")[-1].split()[-1].strip()
            d[k] = v
    return d
    
def get_uniport_tab_dict(tab):
    d = {}
    fu = open(tab)
    head = fu.next().split("\t")
    index = head.index(head[-1]) + 1
    s = os.popen("cut -f%d %s | sort |uniq -u" %(index,tab)).read().strip().split("\n")
    for line in fu:
        if line.strip("\n").split("\t")[-1] in s or not line.strip("\n").split("\t")[-1]:
            continue
        else:
            d[(line.strip("\n").split("\t")[-1]).split()[0]] = line.strip("\n").split("\t")[0]
    return d
        
    
def main():
    with open(sys.argv[4],"w") as anno:
        anno.write("#Gene_id\tEntrez_id\tUniprot_id\tChr\tStart\tEnd\tSymbol\tGene_type\tDescription\n")
        ud = get_uniport_dict(sys.argv[2])
        ut = get_uniport_tab_dict(sys.argv[3])
        for gene in get_gtf_gen_group(sys.argv[1]):
            ref_id = re.search('gene_id "(.+?)";',str(gene)).group(1)
            entrez_id  = re.search('Entrez_id "(.+?)";',str(gene)).group(1) if re.search('Entrez_id "(.+?)";',str(gene)) else "-"
            uniport_id = ud[entrez_id] if ud.has_key(entrez_id) else "-"
            chr = gene[0].split("\t")[0]
            s,e = gene[0].split("\t")[3:5]
            symbol = re.split("_\d+?",re.search('gene_id "(.+?)";',str(gene)).group(1))[0] if re.search('gene_id "(.+?)";',str(gene)) else "-"
            uniport_id = ut[symbol] if uniport_id == "-" and ut.has_key(symbol) else uniport_id
            Gene_type = re.search('gene_biotype "(.+?)";',str(gene)).group(1) if re.search('gene_biotype "(.+?)";',str(gene)) else "-"
            Description = re.search('product "(.+?)";',str(gene)).group(1) if re.search('product "(.+?)";',str(gene)) else "-"
            Description = re.search('Note "(.+?)";',str(gene)).group(1) if Description == "-"  and re.search('Note "(.+?)";',str(gene)) else Description
            anno.write("\t".join([ref_id,entrez_id,uniport_id,chr,s,e,symbol,Gene_type,Description]) + "\n")
            
if __name__ == "__main__":
    main()
            
