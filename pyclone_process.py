# -*- coding: utf-8 -*-
"""
Created on Fri May 29 15:14:15 2020

@author: jxf43
"""

import os
import pandas as pd
from itertools import compress

#input_dir="I:\\xiaojun\\shiguanai"
patient_dir="I:\\xiaojun\\shiguanai"
patient_information=pd.read_csv(os.path.join(patient_dir,"sample_information.csv"))
input_dir="I:\\xiaojun\\shiguanai\\GetBasecounts_result"
files=os.listdir(input_dir)
patientname=[item.split(".")[0] for item in files]
patientname=list(set(patientname))

pre_sample=[item.split("-")[0] for item in patient_information["pre"]]
post_sample=[item.split("-")[0] for item in patient_information["post"]]

for name in patientname:
    os.makedirs("I:\\xiaojun\\shiguanai\\\\pyclone\\"+name)
#os.mkdir("I:\\xiaojun\\shiguanai\\sigprofiler")
file='lanjingxian.GetBasecounts.result.txt' 
file_content=pd.read_csv(os.path.join(input_dir,file),sep="\t",index_col=False)   
#total_mutation_type=[] 
#for file in files:
    #if file.endswith("result.txt"):
        #file_content=pd.read_csv(os.path.join(input_dir,file),sep="\t",index_col=False)
        #for mutation_type in file_content["Variant_Classification"]:
            #total_mutation_type.append(mutation_type)
#total_mutation_type=list(set(total_mutation_type))
silent_mutation=["5'UTR",'Silent',"3'UTR","3'Flank","5'Flank",'Intron','IGR','RNA']    
silent_list=silent_mutation
snv_table=file_content
snv_table=file_content[(file_content["Chromosome"]!="X") | (file_content["Chromosome"]!="Y") ]
cnv_dir="I:\\xiaojun\\shiguanai\\K964_facets\seg"
purity_dir="I:\\xiaojun\\shiguanai\\K964_facets\\purity"


def snv_cnv_match(snv_table,cnv_table,major_copy=2):
    snv_table["mutation_id"]=snv_table["mutation_site"]
    snv_table["ref_counts"]=snv_table["t_ref_count"]
    snv_table["var_counts"]=snv_table["t_alt_count"]
    snv_table["normal_cn"]=2
    major_copy_number=[]
    minor_copy_number=[]
    segment=[]
    cnv_tag=[]
    for i in range(snv_table.shape[0]):
        chromosome=int(snv_table["Chromosome"].iloc[i])
        start=int(snv_table["Start_Position"].iloc[i])
        segment_chromosome=cnv_table[cnv_table["chrom"]==chromosome]
        segment_chromosome=segment_chromosome.reset_index()
        findit=False
        for j in range(segment_chromosome.shape[0]):
            seg_start=int(segment_chromosome.loc[j,"start"])
            seg_end=int(segment_chromosome.loc[j,"end"])
            total_cn=segment_chromosome.loc[j,"tcn.em"]
            minor_cn=segment_chromosome.loc[j,"lcn.em"]
            seg=segment_chromosome.loc[j,"seg"]
            if start>=seg_start and start<=seg_end:
                findit=True
                if (not pd.isnull(total_cn)) and (not pd.isnull(minor_cn)):
                    major_cn=int(total_cn)-int(minor_cn)
                    major_copy_number.append(major_cn)
                    minor_copy_number.append(int(minor_cn))
                    segment.append(seg)
                    cnv_tag.append("correct")
                  
                else:
                    major_copy_number.append(total_cn)
                    minor_copy_number.append(0)
                    segment.append(seg)
                    cnv_tag.append("minor_NA")
        if findit==False and major_copy==2:
           major_copy_number.append(2)
           minor_copy_number.append(0)
           segment.append(9999)
           cnv_tag.append("not_cover_major2")
        elif findit==False and major_copy==1: 
           major_copy_number.append(1)
           minor_copy_number.append(1)
           segment.append(9999)
           cnv_tag.append("not_cover_major2")
    snv_table["minor_cn"]=minor_copy_number
    snv_table["major_cn"]=major_copy_number
    snv_table["segment"]=segment
    snv_table["cnv_tag"]=cnv_tag
    if major_copy==2:
        snv_table=snv_table.replace({'major_cn': {0: 2}})
    return(snv_table)       
                
patient_name="lanjingxian"          
test= snv_cnv_match(snv_table=post_snv_table,cnv_table=post_cnv_table,major_copy=2)  
    
output_dir="I:\\xiaojun\\shiguanai\\pyclone\\"
def pyclone_input(output_dir,patient_name,snv_table,cnv_dir,purity_dir,keep_silent=True,silent_list=silent_mutation,major_copy=2):
    snv_table=snv_table[(snv_table["Chromosome"]!="X") & (snv_table["Chromosome"]!="Y") ]
    #snv_table["gene_id"]=snv_table["Chromosome"].astype(int).astype(str)+"_"+snv_table["Start_Position"].astype(int).astype(str)
    snv_table["mutation_site"]=snv_table["Hugo_Symbol"].astype(str)+"_"+snv_table["gene_id"].astype(str)+"_"+snv_table["HGVSc"].astype(str)
    if(keep_silent==False):
        snv_table=snv_table[~snv_table["Variant_Classification"].isin(silent_list)]
    sample_name=snv_table["Tumor_Sample_Barcode"]    
    sample_name=[item.split("-")[0] for item in sample_name]
    snv_table["sample_name"]=sample_name
    unique_sample_name=list(set(sample_name))
    for name in unique_sample_name:
        if name in pre_sample:
            pre_snv_table=snv_table[snv_table["sample_name"]==name]
            pre_snv_table=pre_snv_table.reset_index()
        elif name in post_sample:
            post_snv_table=snv_table[snv_table["sample_name"]==name]
            post_snv_table=post_snv_table.reset_index()
    cnv_files=os.listdir(cnv_dir)
    purity_files=os.listdir(purity_dir)
    for name in unique_sample_name: 
        if name in pre_sample:
            name_index=[name in item for item in cnv_files]
            cnv_file=list(compress(cnv_files, name_index))
            pre_cnv_table=pd.read_csv(os.path.join(cnv_dir,cnv_file[0]))
        elif name in post_sample:
            name_index=[name in item for item in cnv_files]
            cnv_file=list(compress(cnv_files, name_index))
            post_cnv_table=pd.read_csv(os.path.join(cnv_dir,cnv_file[0]))
    pre_input=snv_cnv_match(snv_table=pre_snv_table,cnv_table=pre_cnv_table,major_copy=major_copy) 
    post_input=snv_cnv_match(snv_table=post_snv_table,cnv_table=post_cnv_table,major_copy=major_copy) 
    target_dir=output_dir+patient_name
    pre_filename=patient_name+"_before.tsv"
    pre_input.to_csv(os.path.join(target_dir,pre_filename),sep="\t",index=False)  
    post_filename=patient_name+"_post.tsv"
    post_input.to_csv(os.path.join(target_dir,post_filename),sep="\t",index=False)
    if list(pre_input["mutation_id"])==list(post_input["mutation_id"]):
        print(patient_name+": success")
    for name in unique_sample_name:
        if name in pre_sample:
            name_index=[name in item for item in purity_files]
            purity_file=list(compress(purity_files, name_index))
            pre_purity_table=pd.read_csv(os.path.join(purity_dir,purity_file[0]),sep="\t") 
            pre_purity=float(pre_purity_table.columns[1])
        if name in post_sample:
            name_index=[name in item for item in purity_files]
            purity_file=list(compress(purity_files, name_index))
            post_purity_table=pd.read_csv(os.path.join(purity_dir,purity_file[0]),sep="\t") 
            post_purity=float(post_purity_table.columns[1]) 
    purity_data={"sample":["pre","post"],"purity":[pre_purity,post_purity]}
    purity_table=pd.DataFrame.from_dict(purity_data) 
    purity_name=patient_name+"_"+"purity.tsv"
    purity_table.to_csv(os.path.join(target_dir,purity_name),sep="\t",index=False)        

            
input_dir_correct="I:\\xiaojun\\shiguanai\\mutation_table_correct"
patient="anchuanguo"
mutation_files=os.listdir(input_dir_correct)
for patient in patientname: 
   p_file=[item for item in mutation_files if patient in item]           
   pcontent=pd.read_csv(os.path.join(input_dir_correct,p_file[0]),sep="\t",index_col=False) 
   pyclone_input(output_dir=output_dir,patient_name=patient,snv_table=pcontent,cnv_dir=cnv_dir,purity_dir=purity_dir,keep_silent=True,silent_list=silent_mutation,major_copy=2)
