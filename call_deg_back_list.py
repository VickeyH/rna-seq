#!/usr/bin/env python
#coding:utf-8

import argparse,os,sys
parser = argparse.ArgumentParser(description="This Script is used to call deg.list, back.list, deg.up_down.list for kobas. 'back.list' will be generate in the parent directory of 'deg.*' file")
parser.add_argument("-diff","--diff_info",type=str,help="The diff info file, like group1-VS-group2.txt",required = True,nargs='*')
parser.add_argument("-anno","--anno",type=str,help="the anno file of ref species",required = True)
parser.add_argument("-t","--dbtype",type=str,help="the type of gene_id used for kobsa, default: id:ncbigene",default="id:ncbigene",choices=["id:ncbigene","id:ensembl","id:uniprot"])
parser.add_argument("-r","--reftype",type=str,help="the type of gene_id in your diff file(first colnum), default: ensembl",default="ensembl",choices=["ensembl","entrez","uniprot"])
args = parser.parse_args()

def get_id_dict(anno,kid,rid):
    kd = {"id:ensembl":"Ensembl_id","id:ncbigene":"Entrez_id","id:uniprot":"Uniprot_id"}
    rd = {"ensembl":"Ensembl_id","entrez":"Entrez_id","id:uniprot":"Uniprot_id"}
    kid = kd[kid]
    rid = rd[rid]
    all_d,id_dict = {},{}
    all_d["Ensembl_id"] = os.popen('cut -f1 %s'%anno).read().strip().split("\n")[1:]
    all_d["Entrez_id"] = os.popen('cut -f2 %s'%anno).read().strip().split("\n")[1:]
    all_d["Uniprot_id"] = os.popen('cut -f3 %s'%anno).read().strip().split("\n")[1:]
    id_pair_list = zip(all_d[rid],all_d[kid])
    id_pair_new = [i for i in id_pair_list if i[0] != "-"]   ## 去除[("-","ENSMUSG00000028185"),("-","ENSMUSG00000028186")..] 这种类似元素
    for i in id_pair_new:
        if "," in i[0]:
            id_pair_new.remove(i)
            for g in i[0].split(","):
                id_pair_new.append((g,i[1]))
    for i in id_pair_new:
        id_dict.setdefault(i[0],[]).append(i[1])     ### 合并重复元素     
    return {k:",".join(v) for k,v in id_dict.iteritems()}  ## 将值转换为字符串
 
def get_anno_dict(anno):
    ann_dict = {}
    with open(anno) as ano:
        for line in ano:
            line_list = line.strip().split("\t")
            ann_dict.setdefault(line_list[0],{})["entrez_id"] = line_list[1]
            ann_dict.setdefault(line_list[0],{})["uniprot_id"] = line_list[2]
            ann_dict.setdefault(line_list[0],{})["ref_pos"] = line_list[3] + ":" + line_list[4] + "-" + line_list[5]
            ann_dict.setdefault(line_list[0],{})["symbol"] = line_list[6]
            ann_dict.setdefault(line_list[0],{})["gene_type"] = line_list[7]
            ann_dict.setdefault(line_list[0],{})["description"] = line_list[8]
    return ann_dict    

def call_deg_list(difffile,anno,outputdir,kid,rid):  ### 生成deg.list deg.up_down.list文件
    id_dict = get_id_dict(anno,kid,rid)
    deg_up = os.popen('grep -E "\tup$" %s | wc -l'%difffile).read().strip()
    deg_down = os.popen('grep -E "\tdown$" %s | wc -l'%difffile).read().strip()
    with open(difffile) as df,open(os.path.join(outputdir,"deg.up_down.list"),"w") as up_down:
        header = df.next()
        for line in df:
            line_list = line.split("\t")
            if not line.endswith("\tno\n") and id_dict.has_key(line_list[0]):
                up_down.writelines([i + "\t" + line_list[-2] + "\t" + line_list[-1] for i in id_dict[line_list[0]].split(",") if id_dict[line_list[0]] != "-"])
    os.system('cut -f1 %s > %s'%(os.path.join(outputdir,"deg.up_down.list"),os.path.join(outputdir,"deg.list")))
    os.system("echo -ne 'Sample\tTotal_DEG\tUp\tDown\n%s\t%d\t%d\t%d\n' > %s"%(os.path.splitext(os.path.basename(difffile))[0],int(deg_up)+int(deg_down),int(deg_up),int(deg_down),os.path.join(outputdir,"deg.num.txt")))
    
def call_back_list(difffile,anno,outputdir,kid,rid):
    id_dict = get_id_dict(anno,kid,rid)
    with open(difffile) as df,open(os.path.join(outputdir,"back.list"),"w") as back:
        header = df.next()
        for line in df:
            if float(line.split("\t")[1]) > 0 and id_dict.has_key(line.split("\t")[0]):
                back.writelines([i + "\n" for i in id_dict[line.split("\t")[0]].split(",") if id_dict[line.split("\t")[0]] != "-"])
    # os.system('sort -u %s -o %s'%(os.path.join(outputdir,"back.list"),os.path.join(outputdir,"back.list")))  ## back.list去重
                
def call_deg_up_down_anno(anno,difffile,outputdir):   ### anno文件添加
    anno_dict = get_anno_dict(anno)
    na_dict = dict.fromkeys(["entrez_id","uniprot_id","ref_pos","symbol","gene_type","description"],"-")
    with open(difffile) as df, open(os.path.join(outputdir,"deg.up.txt"),"w") as du, open(os.path.join(outputdir,"deg.down.txt"),"w") as dd:
        header = df.next()
        new_header = header.split("\t")[0] + "\t" + "\t".join(header.split("\t")[2:-1]) + "\t" + "\t".join(["track","entrez_id","uniprot_id","ref_pos","symbol","gene_type","description"]) + "\n"
        du.write(new_header)
        dd.write(new_header)
        for line in df:
            line_list = line.strip().split("\t")
            if line.endswith("\tdown\n"):
                dd.write(line_list[0] + "\t" + "\t".join(line_list[2:-1]) + "\tdown\t" + anno_dict.get(line_list[0],na_dict)["entrez_id"] + "\t" + anno_dict.get(line_list[0],na_dict)["uniprot_id"] + "\t" + anno_dict.get(line_list[0],na_dict)["ref_pos"] + "\t" + anno_dict.get(line_list[0],na_dict)["symbol"] + "\t" + anno_dict.get(line_list[0],na_dict)["gene_type"] + "\t" + anno_dict.get(line_list[0],na_dict)["description"] + "\n")
            elif line.endswith("\tup\n"):
                du.write(line_list[0] + "\t" + "\t".join(line_list[2:-1]) + "\tup\t" + anno_dict.get(line_list[0],na_dict)["entrez_id"] + "\t" + anno_dict.get(line_list[0],na_dict)["uniprot_id"] + "\t" + anno_dict.get(line_list[0],na_dict)["ref_pos"] + "\t" + anno_dict.get(line_list[0],na_dict)["symbol"] + "\t" + anno_dict.get(line_list[0],na_dict)["gene_type"] + "\t" + anno_dict.get(line_list[0],na_dict)["description"] + "\n")
            else:
                continue

if __name__ == "__main__":
    anno = args.anno
    if not os.path.exists(anno):
        print "%s file not exists!" %anno
        sys.exit(0)
    kid = args.dbtype
    rid = args.reftype
    difffiles = args.diff_info
    for f in difffiles:
        if not os.path.exists(f):
            print "%s file not exists!" %f
            sys.exit(0)
    outdir = os.path.dirname(os.path.dirname(os.path.abspath(difffiles[0])))
    call_back_list(difffiles[0],anno,outdir,kid,rid)
    for diff in difffiles:
        call_deg_list(diff,anno,os.path.dirname(diff),kid,rid)
        call_deg_up_down_anno(anno,diff,os.path.dirname(diff))
        
