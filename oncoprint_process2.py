# -*- coding: utf-8 -*-
"""
Created on Fri Jun  5 09:58:45 2020

@author: jxf43
"""

import os
import pandas as pd
from itertools import chain

os.chdir("I:\\xiaojun\\shiguanai\\oncoprint_table")
os.listdir(".")

mutation_oncoprint=pd.read_csv("mutation_oncoprint_pre.csv",index_col=False)
gene_pathway=pd.read_csv("oncogenic_pathway.csv",index_col=False)
cnv_oncoprint_pre=pd.read_csv("gene_cnv_oncoprint_pre.csv",index_col=False)

oncoprint_table=mutation_oncoprint
pathway_table=gene_pathway
def  join_mutation(a_list):
    a_list=a_list[a_list.notnull()]
    a_list=list(set(a_list))
    a_list=["SNV" for item in a_list]
    a_list=list(set(a_list))
    join_result=(",").join(a_list)
    return (join_result) 

def  join_cnv(a_list):
    a_list=a_list[a_list.notnull()]
    a_list=list(set(a_list))
    a_list=["SCNA" for item in a_list]
    a_list=list(set(a_list))
    join_result=(",").join(a_list)
    return (join_result)

oncoprint_type=["oncogenic","total"]
oncoprint_type="oncogenic"
mutation_cnv="mutation"
keep_oncogenic=2
mutation_cnv=["mutation","cnv]
keep_pathway=5
def gene_pathway_match(oncoprint_table, pathway_table,mutation_cnv,oncoprint_type="oncogenic",keep_oncogenic=10,keep_pathway=5):
    return_result={}
    oncoprint_table=pd.merge(oncoprint_table,pathway_table,on="Gene",how="left")
    oncoprint_table["Pathway"]=oncoprint_table["Pathway"].fillna("Others")
    #oncoprint_table["]
    return_result["original_table"]=oncoprint_table
    if oncoprint_type=="oncogenic":
        keep_gene=list(oncoprint_table[(oncoprint_table["count"]>=0) & (oncoprint_table["Pathway"] !="Others")]["Gene"])
        keep_gene=keep_gene[0:keep_oncogenic]
    if oncoprint_type=="total":
        keep_gene=list(oncoprint_table[oncoprint_table["count"]>=0]["Gene"])
        keep_gene=keep_gene[0:keep_oncogenic]
    #keep_gene=list(set(keep_gene+keep_gene2))
    #keep_gene3=list(oncoprint_table[(oncoprint_table["count"]>=2) & (oncoprint_table["Pathway"] !="Others")]["Gene"])
    oncoprint_table_keep_gene=oncoprint_table[oncoprint_table["Gene"].isin(keep_gene)]
    oncoprint_table_keep_gene=oncoprint_table_keep_gene.loc[:,['Gene', 'anchuanguo', 'chenguoquan', 'chenxingming', 'chenzongye', 'lanjingxian', 'taoxiaokun', 'yangbaoxian', 'zhanggaohua', 'zhanqingping', 'zhengyoufu']]
    pathway_oncoprint_table=oncoprint_table.loc[:,['Pathway', 'anchuanguo', 'chenguoquan', 'chenxingming', 'chenzongye', 'lanjingxian', 'taoxiaokun', 'yangbaoxian', 'zhanggaohua', 'zhanqingping', 'zhengyoufu']]
    function_dict={}
    if mutation_cnv=="mutation":
        for patient in pathway_oncoprint_table.columns[1:11]:
            function_dict[patient]=join_mutation
    if mutation_cnv=="cnv":
        for patient in pathway_oncoprint_table.columns[1:11]:
            function_dict[patient]=join_cnv        
    pathway_oncoprint_table_grouped = pathway_oncoprint_table.groupby("Pathway")     
    pathway_oncoprint_table_grouped = pathway_oncoprint_table_grouped.aggregate(function_dict)
    pathway_oncoprint_table_grouped=pathway_oncoprint_table_grouped.reset_index()
    count=[]
    if mutation_cnv=="mutation":
        for i in range(pathway_oncoprint_table_grouped.shape[0]):
            num=(pathway_oncoprint_table_grouped.iloc[i,1:11]=="SNV").sum()
            count.append(num)
    if mutation_cnv=="cnv":
        for i in range(pathway_oncoprint_table_grouped.shape[0]):
            num=(pathway_oncoprint_table_grouped.iloc[i,1:11]=="SCNA").sum()
            count.append(num)        
    pathway_oncoprint_table_grouped["count"]= count   
    pathway_oncoprint_table_grouped =pathway_oncoprint_table_grouped.sort_values(by=['count'], ascending=False)
    return_result["original_pathway_table"]=pathway_oncoprint_table_grouped
    keep_pathway_table=pathway_oncoprint_table_grouped[(pathway_oncoprint_table_grouped["count"]>=0) & (pathway_oncoprint_table_grouped["Pathway"]!="Others")]
    keep_pathway_table=keep_pathway_table.iloc[0:keep_pathway,]
    keep_pathway_table=keep_pathway_table.iloc[:,0:11]    
    keep_pathway_table["Pathway"]=keep_pathway_table["Pathway"]+"_p"
    keep_pathway_table=keep_pathway_table.rename(columns={"Pathway": "Gene"})
    gene_pathway_table=pd.concat([oncoprint_table_keep_gene,keep_pathway_table])
    #gene_pathway_table_oncogenic=gene_pathway_table[gene_pathway_table["Gene"].isin(keep_gene3)]
    return_result["gene_pathway_table"]=gene_pathway_table
    #return_result["gene_pathway_table_oncogenic"]=gene_pathway_table_oncogenic
    return(return_result)
    
mutation_oncoprint_pre=return_result

cnv_oncoprint_pre_table=gene_pathway_match(oncoprint_table=cnv_oncoprint_pre, pathway_table=gene_pathway,keep_pathway=4)

a_data_frame=gene_pathway_table.drop(["Gene"], axis=1)

mutation_oncoprint_pre_neo=gene_pathway_match(oncoprint_table=mutation_oncoprint, pathway_table=gene_pathway,mutation_cnv="mutation",oncoprint_type="oncogenic",keep_oncogenic=10,keep_pathway=5)    
cnv_oncoprint_pre_neo=gene_pathway_match(oncoprint_table=cnv_oncoprint_pre, pathway_table=gene_pathway,mutation_cnv="cnv",oncoprint_type="oncogenic",keep_oncogenic=5,keep_pathway=3) 

pre_oncoprint=mutation_oncoprint_pre_neo["gene_pathway_table"]
pre_oncoprint_cnv=cnv_oncoprint_pre_neo["gene_pathway_table"]
#pre_oncoprint_cnv["Gene"]=pre_oncoprint_cnv["Gene"]+"_C"
pre_oncoprint_total=pd.concat([pre_oncoprint,pre_oncoprint_cnv]) 

arm_cnv_pre=pd.read_csv("arm_cnv_oncoprint_pre.csv",index_col=False) 
arm_cnv_pre=arm_cnv_pre.iloc[0:5,] 
arm_cnv_pre=arm_cnv_pre.iloc[:,0:11]
arm_cnv_pre=arm_cnv_pre.rename(columns={"Arm": "Gene"})
pre_oncoprint_total=pd.concat([pre_oncoprint_total,arm_cnv_pre])

pre_oncoprint_total.to_csv("pre_oncoprint_total.csv",index=False)


mutation_oncoprint_post=pd.read_csv("mutation_oncoprint_post.csv",index_col=False)
#gene_pathway=pd.read_csv("oncogenic_pathway.csv",index_col=False)
cnv_oncoprint_post=pd.read_csv("gene_cnv_oncoprint_post.csv",index_col=False)

mutation_oncoprint_post_neo=gene_pathway_match(oncoprint_table=mutation_oncoprint_post, pathway_table=gene_pathway,mutation_cnv="mutation",oncoprint_type="oncogenic",keep_oncogenic=10,keep_pathway=5)    
cnv_oncoprint_post_neo=gene_pathway_match(oncoprint_table=cnv_oncoprint_post, pathway_table=gene_pathway,mutation_cnv="cnv",oncoprint_type="oncogenic",keep_oncogenic=5,keep_pathway=3) 

post_oncoprint=mutation_oncoprint_post_neo["gene_pathway_table"]
post_oncoprint_cnv=cnv_oncoprint_post_neo["gene_pathway_table"]
#pre_oncoprint_cnv["Gene"]=pre_oncoprint_cnv["Gene"]+"_C"
post_oncoprint_total=pd.concat([post_oncoprint,post_oncoprint_cnv]) 

arm_cnv_post=pd.read_csv("arm_cnv_oncoprint_post.csv",index_col=False) 
arm_cnv_post=arm_cnv_post.iloc[0:5,] 
arm_cnv_post=arm_cnv_post.iloc[:,0:11]
arm_cnv_post=arm_cnv_post.rename(columns={"Arm": "Gene"})
post_oncoprint_total=pd.concat([post_oncoprint_total,arm_cnv_post])
post_oncoprint_total.to_csv("post_oncoprint_total.csv",index=False)


mutation_oncoprint_pre_neo["original_table"].to_csv("pre_gene_mutation_table.csv",index=False)
mutation_oncoprint_pre_neo["original_pathway_table"].to_csv("pre_pathway_table.csv",index=False)

cnv_oncoprint_pre_neo["original_table"].to_csv("pre_gene_cnv_table.csv",index=False)
cnv_oncoprint_pre_neo["original_pathway_table"].to_csv("pre_cnv_pathway_table.csv",index=False)

mutation_oncoprint_post_neo["original_table"].to_csv("post_gene_mutation_table.csv",index=False)
mutation_oncoprint_post_neo["original_pathway_table"].to_csv("post_pathway_table.csv",index=False)

cnv_oncoprint_post_neo["original_table"].to_csv("post_gene_cnv_table.csv",index=False)
cnv_oncoprint_post_neo["original_pathway_table"].to_csv("post_cnv_pathway_table.csv",index=False)

pre_table=cnv_oncoprint_pre_neo["original_pathway_table"]
post_table=cnv_oncoprint_post_neo["original_pathway_table"]
import scipy.stats as stats

def fisher_test(pre_table,post_table):
    genes=[]
    p_values=[]
    gene_list=list(set(list(pre_table["Pathway"])+list(post_table["Pathway"])))
    for gene in gene_list:
        if gene in list(pre_table["Pathway"]):
            gene_pre_count=list(pre_table[pre_table["Pathway"]==gene].loc[:,"count"])[0]
        else:
            gene_pre_count=0
        if gene in list(post_table["Pathway"]):
            gene_post_count=list(post_table[post_table["Pathway"]==gene].loc[:,"count"])[0]
        else:
            gene_post_count=0
        if gene_pre_count+gene_post_count>=4:
            pvalue = stats.fisher_exact([[10-gene_pre_count, gene_pre_count], [10-gene_post_count, gene_post_count]])[1]
            genes.append(gene)
            p_values.append(pvalue)

def get_unique_value(a_data_frame):
    values=[]
    for column in a_data_frame.columns:
        value=list(a_data_frame[column].unique())
        values=values+value
    values=list(set(values))
    values = [x for x in values if str(x) != 'nan']
    values= [x for x in values if str(x) != '']
    values_1=[item.split(",") for item in values if "," in item]
    values_2=[item for item in values if "," not in item]
    output=[]
    for v in values_1:
        for j in v:
            output.append(j)
    values_total=output+values_2
    #result=reemoNestings(values_total)
    result=list(set(values_total))
    return (result)
        
alterations= get_unique_value(pre_oncoprint_total.iloc[:,1:11])  

colorcode=pd.read_csv("colorcode.csv")

color=colorcode
variations=alterations
start=0
step=15
variation=alterations[0]
background_color='#CCCCCC'



def onprint_color(variations,color=colorcode,start=0,step=15):
    first="col=c("
    for variation in variations:
        variation_color="{0}='{1}'".format(variation,list(colorcode["colorcode"])[start])
        first=first+variation_color+","
        start=start+step
    return(first)    

def oncoprint_color(variations,background_color='#CCCCCC',color=colorcode,start=0,step=15):
    background="alter_fun = list(background = alter_graphic('rect', fill = '{0}'),".format(background_color)
    for variation in variations:
        variation_color="{0} = alter_graphic('rect', fill = col['{1}']),".format(variation,variation)
        background=background+variation_color
        start=start+step
    return(background)  

    
 gene_pathway_table.to_csv("gene_pathway_table_oncoprint.csv",index=False)       
 
def oncoprint_color_2(variations,background_color='#CCCCCC',color=colorcode,start=0,step=15):
    background="alter_fun = list(background = function(x, y, w, h) {grid.rect(x, y, w-unit(2, 'pt'), h-unit(2, 'pt'), gp = gpar(fill = '%s'))}," %background_color
    for variation in variations:
        variation_color="%s" %variation +"="+"function(x, y, w, h) {grid.rect(x, y, w-unit(2, 'pt'), h-unit(2, 'pt'), gp = gpar(fill = '%s'))}," %(list(color["colorcode"])[start])
        background=background+variation_color
        start=start+step
    return(background)  
    