# -*- coding: utf-8 -*-
"""
Created on Sun May 31 13:51:17 2020

@author: jxf43
"""

import os
import pandas as pd
from itertools import compress

input_dir="I:\\xiaojun\\shiguanai"
patient_information=pd.read_csv(os.path.join(input_dir,"sample_information.csv"))
input_dir="I:\\xiaojun\\shiguanai\\mutation_table_correct"
files=os.listdir(input_dir)
patientname=[item.split(".")[0] for item in files]
patientname=list(set(patientname))

pre_sample=[item.split("-")[0] for item in patient_information["pre"]]
post_sample=[item.split("-")[0] for item in patient_information["post"]]

#for name in patientname:
    #os.makedirs("I:\\xiaojun\\shiguanai\\\\pyclone\\"+name)
#os.mkdir("I:\\xiaojun\\shiguanai\\sigprofiler")
file='lanjingxian.GetBasecounts.result.txt' 
file_content=pd.read_csv(os.path.join(input_dir,file),sep="\t",index_col=False)
snv_table=file_content
sample_name=file_content["Tumor_Sample_Barcode"]    
sample_name=[item.split("-")[0] for item in sample_name]
snv_table["sample_name"]=sample_name
unique_sample_name=list(set(sample_name))

silent_mutation=["5'UTR",'Silent',"3'UTR","3'Flank","5'Flank",'Intron','IGR','RNA']  

def mutation_compare(snv_table,patientname):
    sample_name=snv_table["Tumor_Sample_Barcode"]
    sample_name=[item.split("-")[0] for item in sample_name]
    snv_table["sample_name"]=sample_name
    unique_sample_name=list(set(sample_name))
    for name in unique_sample_name:
        if name in pre_sample:
            pre_snv_table=snv_table[snv_table["sample_name"]==name]
            pre_snv_table=pre_snv_table.reset_index()
        if name in post_sample:
            post_snv_table=snv_table[snv_table["sample_name"]==name]
            post_snv_table=post_snv_table.reset_index() 
    if pre_snv_table.shape[0]!=post_snv_table.shape[0]  :
        print (patientname+":"+"pre and post have different rows")
    elif pre_snv_table.shape[0]==post_snv_table.shape[0] :
        pre_snv_table_reduce=pre_snv_table.loc[:,["Hugo_Symbol","Chromosome","Start_Position","Reference_Allele","Tumor_Seq_Allele1"]]
        post_snv_table_reduce=post_snv_table.loc[:,["Hugo_Symbol","Chromosome","Start_Position","Reference_Allele","Tumor_Seq_Allele1"]]
        tag=[]
        wrong_row=[]
        for i in range(pre_snv_table.shape[0]):
            if (list(pre_snv_table_reduce.loc[i,:])==list(post_snv_table_reduce.loc[i,:])):
                tag.append(1)
            else:
                tag.append(0)
                wrong_row.append(pre_snv_table_reduce.loc[i,:])
        if sum(tag)==pre_snv_table_reduce.shape[0]:
           print (patientname+":"+"success")
        else:
            result={"patientname":patientname,"wrong_rows":wrong_row}
            return(result)
        
check_result=[]           
for file in files:
    if file.endswith("result.txt"):
       patient_name=file.split(".") [0]
       fc= pd.read_csv(os.path.join(input_dir,file),sep="\t",index_col=False)
       result=mutation_compare(snv_table=fc,patientname=patient_name) 
       check_result.append(result)

def  join_mutation(a_list):
    join_result=(",").join(a_list)
    return (join_result)            

def mutation_oncoprint_extraction(snv_table,patientname,minumum_alt_read=1,silent_list=silent_mutation):
    final_result=[]
    sample_name=snv_table["Tumor_Sample_Barcode"]
    sample_name=[item.split("-")[0] for item in sample_name]
    snv_table["sample_name"]=sample_name
    unique_sample_name=list(set(sample_name))
    for name in unique_sample_name:
        if name in pre_sample:
            pre_snv_table=snv_table[snv_table["sample_name"]==name]
            pre_snv_table=pre_snv_table.reset_index()
        if name in post_sample:
            post_snv_table=snv_table[snv_table["sample_name"]==name]
            post_snv_table=post_snv_table.reset_index() 
    if pre_snv_table.shape[0]!=post_snv_table.shape[0]  :
        print (patientname+":"+"pre and post have different rows")
    elif pre_snv_table.shape[0]==post_snv_table.shape[0] :
        pre_snv_table_reduce=pre_snv_table[pre_snv_table["t_alt_count"]>=minumum_alt_read]
        pre_snv_table_reduce=pre_snv_table_reduce.loc[:,["Hugo_Symbol","Variant_Classification"]]
        pre_snv_table_reduce=pre_snv_table_reduce[~pre_snv_table_reduce["Variant_Classification"].isin(silent_list)]
        pre_snv_table_reduce=pre_snv_table_reduce.reset_index()
        pre_snv_table_reduce=pre_snv_table_reduce.loc[:,["Hugo_Symbol","Variant_Classification"]]
        pre_snv_table_reduce=pre_snv_table_reduce.drop_duplicates()
        pre_snv_grouped = pre_snv_table_reduce.groupby('Hugo_Symbol')
        #print (pre_snv_grouped.groups)
        #aggregate:predefined function (heres is joint_mutation) acted on the selected column (here is Variant_Classification) of grouped dataframe
        # that is, the parameter in the function is a column (serial)
        # apply: predefined function acts on every value of selected columns
        # that is, the parameter in the function is a value
        pre_snv_oncoprint=pre_snv_grouped.aggregate({'Variant_Classification':join_mutation})
        pre_snv_oncoprint=pre_snv_oncoprint.reset_index()
        pre_snv_oncoprint.columns=["Gene",patientname]
        
        post_snv_table_reduce=post_snv_table[post_snv_table["t_alt_count"]>=minumum_alt_read]
        post_snv_table_reduce=post_snv_table_reduce.loc[:,["Hugo_Symbol","Variant_Classification"]]
        post_snv_table_reduce=post_snv_table_reduce[~post_snv_table_reduce["Variant_Classification"].isin(silent_list)]
        post_snv_table_reduce=post_snv_table_reduce.reset_index()
        post_snv_table_reduce=post_snv_table_reduce.loc[:,["Hugo_Symbol","Variant_Classification"]]
        post_snv_table_reduce=post_snv_table_reduce.drop_duplicates()
        post_snv_grouped = post_snv_table_reduce.groupby('Hugo_Symbol')
        #print (post_snv_grouped.groups)
        #aggregate:postdefined function (heres is joint_mutation) acted on the selected column (here is Variant_Classification) of grouped dataframe
        # that is, the parameter in the function is a column (serial)
        # apply: postdefined function acts on every value of selected columns
        # that is, the parameter in the function is a value
        post_snv_oncoprint=post_snv_grouped.aggregate({'Variant_Classification':join_mutation})
        post_snv_oncoprint=post_snv_oncoprint.reset_index()
        post_snv_oncoprint.columns=["Gene",patientname]
        
        final_result.append(pre_snv_oncoprint)
        final_result.append(post_snv_oncoprint)
        return(final_result)
        
pre_oncoprint_total=[] 
post_oncoprint_total=[]          
for file in files:
    if file.endswith("result.txt"):
       patient_name=file.split(".") [0]
       fc= pd.read_csv(os.path.join(input_dir,file),sep="\t",index_col=False)
       result=mutation_oncoprint_extraction(snv_table=fc,patientname=patient_name) 
       pre_oncoprint_total.append(result[0])
       post_oncoprint_total.append(result[1])

pre_oncoprint_combine=pre_oncoprint_total[0]
for i in range(1,10) :
    pre_oncoprint_combine=pd.merge(pre_oncoprint_combine,pre_oncoprint_total[i],on="Gene",how="outer")  

#pre_oncoprint_combine.to_csv("pre_oncoprint_combine.csv")     
count=[]
for i in range(pre_oncoprint_combine.shape[0]):
    num=10-pre_oncoprint_combine.iloc[i,1:11].isna().sum()
    count.append(num)
pre_oncoprint_combine["count"] = count   
pre_oncoprint = pre_oncoprint_combine.sort_values(by=['count'], ascending=False)
pre_oncoprint.to_csv("mutation_oncoprint_pre.csv")

post_oncoprint_combine=post_oncoprint_total[0]
for i in range(1,10) :
    post_oncoprint_combine=pd.merge(post_oncoprint_combine,post_oncoprint_total[i],on="Gene",how="outer")  

#pre_oncoprint_combine.to_csv("pre_oncoprint_combine.csv")     
count=[]
for i in range(post_oncoprint_combine.shape[0]):
    num=10-post_oncoprint_combine.iloc[i,1:11].isna().sum()
    count.append(num)
post_oncoprint_combine["count"] = count   
post_oncoprint = post_oncoprint_combine.sort_values(by=['count'], ascending=False)
post_oncoprint.to_csv("mutation_oncoprint_post.csv")

all_patient_mutation=[]
for file in files:
    pc=pd.read_csv(os.path.join(input_dir,file),sep="\t",index_col=False)
    tumor_sample=[item.split("-")[0] for item in pc['Tumor_Sample_Barcode']]
    pc["sample_name"]=tumor_sample
    all_patient_mutation.append(pc)
all_patient_mutation=pd.concat(all_patient_mutation)   
all_patient_mutation=pd.merge(all_patient_mutation,patient_sample_time,on="sample_name")
all_patient_mutation.to_csv("all_patient_mutation.csv")

patient_sample_table=[]
for patient in patientname:
    for file in files:
        if patient in file: 
            pc=pd.read_csv(os.path.join(input_dir,file),sep="\t",index_col=False)
            samplename=[item.split("-")[0] for item in pc["Tumor_Sample_Barcode"]]
            samplename=list(set(samplename))
            sample_time=[1,2]
            if samplename[0] in pre_sample:
                sample_time[0]="pre_sample"
            else:
                sample_time[0]="post_sample" 
            if samplename[1] in pre_sample:
                sample_time[1]="pre_sample"
            else:
                sample_time[1]="post_sample" 
            result={"Patientname":[patient,patient],"sample_name":samplename,"sample_time":sample_time} 
            result=pd.DataFrame.from_dict(result)
            patient_sample_table.append(result)
patient_sample_time = pd.concat(patient_sample_table)

patient_sample_time.to_csv("patient_sample_time.csv")

gene_cnv_dir="I:\\xiaojun\\shiguanai\\K964_facets\\gene"
gene_cnv_files=os.listdir(gene_cnv_dir)
gene_cnv_combine=[]
for file in gene_cnv_files:
    sam_name=file.split("-")[0]
    sam_content=pd.read_csv(os.path.join(gene_cnv_dir,file),sep="\t",index_col=False)
    sam_content["sample_name"]=sam_name
    gene_cnv_combine.append(sam_content)

gene_cnv_combine=pd.concat(gene_cnv_combine)    
gene_cnv_combine_annotation=pd.merge(gene_cnv_combine,patient_sample_time,on="sample_name")

gene_cnv_pre=gene_cnv_combine_annotation[gene_cnv_combine_annotation["sample_time"]=="pre_sample"]

gene_cnv_oncoprint=gene_cnv_pre.pivot(index='GeneName', columns='Patientname', values='CNVType')
gene_cnv_oncoprint["chenzongye"]=None
count=[]
for i in range(gene_cnv_oncoprint.shape[0]):
    num=10-gene_cnv_oncoprint.iloc[i,0:10].isna().sum()
    count.append(num)
gene_cnv_oncoprint["count"] = count
gene_cnv_oncoprint_pre = gene_cnv_oncoprint.sort_values(by=['count'], ascending=False)


gene_cnv_post=gene_cnv_combine_annotation[gene_cnv_combine_annotation["sample_time"]=="post_sample"]
gene_cnv_post_oncoprint=gene_cnv_post.pivot(index='GeneName', columns='Patientname', values='CNVType')

count=[]
for i in range(gene_cnv_post_oncoprint.shape[0]):
    num=10-gene_cnv_post_oncoprint.iloc[i,0:10].isna().sum()
    count.append(num)
gene_cnv_post_oncoprint["count"] = count
gene_cnv_post_oncoprint = gene_cnv_post_oncoprint.sort_values(by=['count'], ascending=False)

gene_cnv_oncoprint_pre.to_csv("gene_cnv_oncoprint_pre.csv")
gene_cnv_post_oncoprint.to_csv("gene_cnv_oncoprint_post.csv")
gene_cnv_combine_annotation.to_csv("gene_cnv_combine_annotation.csv")




arm_cnv_dir="I:\\xiaojun\\shiguanai\\K964_facets\\arm"
arm_cnv_files=os.listdir(arm_cnv_dir)
arm_cnv_combine=[]
for file in arm_cnv_files:
    sam_name=file.split("-")[0]
    sam_content=pd.read_csv(os.path.join(arm_cnv_dir,file),sep="\t",index_col=False)
    sam_content["sample_name"]=sam_name
    arm_cnv_combine.append(sam_content)

arm_cnv_combine=pd.concat(arm_cnv_combine)    
arm_cnv_combine_annotation=pd.merge(arm_cnv_combine,patient_sample_time,on="sample_name")

arm_cnv_pre=arm_cnv_combine_annotation[arm_cnv_combine_annotation["sample_time"]=="pre_sample"]

arm_cnv_oncoprint=arm_cnv_pre.pivot(index='Arm', columns='Patientname', values='Type')

count=[]
for i in range(arm_cnv_oncoprint.shape[0]):
    num=10-arm_cnv_oncoprint.iloc[i,0:10].isna().sum()
    count.append(num)
arm_cnv_oncoprint["count"] = count
arm_cnv_oncoprint_pre = arm_cnv_oncoprint.sort_values(by=['count'], ascending=False)


arm_cnv_post=arm_cnv_combine_annotation[arm_cnv_combine_annotation["sample_time"]=="post_sample"]
arm_cnv_post_oncoprint=arm_cnv_post.pivot(index='Arm', columns='Patientname', values='Type')

count=[]
for i in range(arm_cnv_post_oncoprint.shape[0]):
    num=10-arm_cnv_post_oncoprint.iloc[i,0:10].isna().sum()
    count.append(num)
arm_cnv_post_oncoprint["count"] = count
arm_cnv_post_oncoprint = arm_cnv_post_oncoprint.sort_values(by=['count'], ascending=False)

arm_cnv_oncoprint_pre.to_csv("arm_cnv_oncoprint_pre.csv")
arm_cnv_post_oncoprint.to_csv("arm_cnv_oncoprint_post.csv")
arm_cnv_combine_annotation.to_csv("arm_cnv_combine_annotation.csv")

input_dir="I:\\xiaojun\\shiguanai\\oncoprint_table"
patient_sample_time=pd.read_csv(os.path.join(input_dir,"patient_sample_time.csv"),index_col=False)

input_dir="I:\\xiaojun\\shiguanai\\K964_facets\\cis"
os.listdir(input_dir)
cis_total=[]
samplenames=[]
for file in os.listdir(input_dir):
    pc=pd.read_csv(os.path.join(input_dir,file),sep="\t")
    cis=float(pc.columns[0].strip("%"))/100
    cis_total.append(cis)
    sample_name=file.split("-")[0]
    samplenames.append(sample_name)
cis_dic={"sample_name":samplenames,"cis":cis_total}    
cis_all=pd.DataFrame(cis_dic)
patient_sample_time=pd.merge(patient_sample_time,cis_all,on="sample_name")
pre_cis=patient_sample_time[patient_sample_time["sample_time"]=="pre_sample"]["cis"]
post_cis=patient_sample_time[patient_sample_time["sample_time"]=="post_sample"]["cis"]
from scipy.stats import wilcoxon
wilcoxon(pre_cis,post_cis)

patient_sample_time.to_csv("patient_sample_time_cis.csv")

input_dir="I:\\xiaojun\\shiguanai\\pyclone\\pyclone_input_coding_mutation"
patient=os.listdir(input_dir)

output_dir="I:\\xiaojun\\shiguanai\\donosGP"
for p in patient:
    p_dir=input_dir+"\\"+p
    p_file=os.listdir(p_dir)
    for pf in p_file:
        if "before" in pf:
            before=pd.read_csv(os.path.join(p_dir,pf),sep="\t",index_col=False)
            before["TIME"]=0
            before["R"]=before["t_ref_count"]+before["t_alt_count"]
            before["r"]=before["t_alt_count"]
            before["MUTID"]=before["gene_id"]
            before["SAMPLEID"]=before["sample_name"]
            before["CNn"]=before["normal_cn"]
            before["CNt"]=before['minor_cn']+before['major_cn']
        if "post" in pf:    
            after=pd.read_csv(os.path.join(p_dir,pf),sep="\t",index_col=False)
            after["TIME"]=1
            after["R"]=after["t_ref_count"]+after["t_alt_count"]
            after["r"]=after["t_alt_count"]
            after["MUTID"]=after["gene_id"]
            after["SAMPLEID"]=after["sample_name"]
            after["CNn"]=after["normal_cn"]
            after["CNt"]=after['minor_cn']+after['major_cn']
    for pf in p_file:
        if "purity" in pf:
            purity=pd.read_csv(os.path.join(p_dir,pf),sep="\t",index_col=False)
            pre_purity=purity[purity["sample"]=="pre"]["purity"][0]
            before["PURITY"]=pre_purity
            after_purity=purity[purity["sample"]=="post"]["purity"][1]
            after["PURITY"]=after_purity   
    final=pd.concat([before,after])
    final["patientname"] = p
    final.to_csv(os.path.join(output_dir,p+".csv"),index=False)        