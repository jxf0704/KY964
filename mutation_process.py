# -*- coding: utf-8 -*-
"""
Created on Mon Jun  1 14:25:22 2020

@author: jxf43
"""

import os
import pandas as pd
from itertools import compress

input_dir="I:\\xiaojun\\shiguanai"
patient_information=pd.read_csv(os.path.join(input_dir,"sample_information.csv"))
input_dir="I:\\xiaojun\\shiguanai\\result_frank\\remove_indel_0"
files=os.listdir(input_dir)
patientname=[item.split("_")[0] for item in files]
patientname=[item.split(".")[0] for item in patientname]
patientname=list(set(patientname))

patient="lanjingxian"
patient_gene={}
for patient in patientname:
    patient_files=[item for item in files if patient in item]
    file_content=[]
    for pf in patient_files:
        fc=pd.read_csv(os.path.join(input_dir,pf),sep="\t",index_col=False)
        file_content.append(fc)
    patient_total=pd.concat(file_content)
    patient_total["gene_id"]=patient_total["Chromosome"].astype(str)+"_"+patient_total["Start_Position"].astype(int).astype(str)
    keep_mutation=list(set(list(patient_total["gene_id"])))
    patient_gene[patient]=keep_mutation
    
input_dir="I:\\xiaojun\\shiguanai\\result_frank\\all_info"
files=os.listdir(input_dir)
output_dir="I:\\xiaojun\\shiguanai\\mutation_table_correct" 

for patient in patientname: 
    patient_files=[item for item in files if patient in item]
    p_content=pd.read_csv(os.path.join(input_dir,patient_files[0]),sep="\t",index_col=False)
    p_content["gene_id"]=p_content["Chromosome"].astype(str)+"_"+p_content["Start_Position"].astype(int).astype(str)
    p_content_filter=p_content[p_content['gene_id'].isin(patient_gene[patient])]
    gene_count=p_content_filter["gene_id"].value_counts()    
    if sum(gene_count==2)==len(gene_count):
        print(patient+":"+"success")
        p_content_filter.loc[p_content_filter['HGVSc'].isnull(),'HGVSc'] = p_content_filter['Variant_Type']
        outputfile=patient+".result.txt"
        p_content_filter.to_csv(os.path.join(output_dir,outputfile),sep="\t",index=False)
        