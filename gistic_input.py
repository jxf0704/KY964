# -*- coding: utf-8 -*-
"""
Created on Sun Jun 14 10:41:41 2020

@author: jxf43
"""

import os
import pandas as pd
os.chdir("I:\\xiaojun\\shiguanai\\exome_signature")

oncogenic_gene=pd.read_csv("total_oncogenic_gene.csv")
[item.split(")")[1] for item in oncogenic_gene["change"]]

import os
import pandas as pd
from itertools import compress
import math

#input_dir="I:\\xiaojun\\shiguanai"
patient_dir="I:\\xiaojun\\shiguanai"
patient_information=pd.read_csv(os.path.join(patient_dir,"sample_information.csv"))
input_dir="I:\\xiaojun\\shiguanai\\GetBasecounts_result"
files=os.listdir(input_dir)
patientname=[item.split(".")[0] for item in files]
patientname=list(set(patientname))
pre_sample=[item.split("-")[0] for item in patient_information["pre"]]
post_sample=[item.split("-")[0] for item in patient_information["post"]]
cnv_dir="I:\\xiaojun\\shiguanai\\K964_facets\seg"

cnv_files=os.listdir(cnv_dir)

pre_gistic_input1=[]
pre_gistic_input2=[]
post_gistic_input1=[]
post_gistic_input2=[]

file=cnv_files[0]
for file in cnv_files:
    sample_name=file.split("-")[0]
    if sample_name in pre_sample:
        seg_cnv=pd.read_csv(os.path.join(cnv_dir,file))
        gistic_input1=seg_cnv[["chrom","start","end","num.mark","cnlr.median"]]
        gistic_input1["sample"]=sample_name
        gistic_input1=gistic_input1[["sample","chrom","start","end","num.mark","cnlr.median"]]
        pre_gistic_input1.append(gistic_input1)
        gistic_input2=seg_cnv[["chrom","start","end","num.mark","tcn.em"]]
        gistic_input2["sample"]=sample_name
        gistic_input2["tcn.em"]=[math.log2(item+0.01) for item in gistic_input2["tcn.em"]]
        gistic_input2=gistic_input2[["sample","chrom","start","end","num.mark","tcn.em"]]
        pre_gistic_input2.append(gistic_input2)
    if sample_name in post_sample:
        seg_cnv=pd.read_csv(os.path.join(cnv_dir,file))
        gistic_input1=seg_cnv[["chrom","start","end","num.mark","cnlr.median"]]
        gistic_input1["sample"]=sample_name
        gistic_input1=gistic_input1[["sample","chrom","start","end","num.mark","cnlr.median"]]
        post_gistic_input1.append(gistic_input1)
        gistic_input2=seg_cnv[["chrom","start","end","num.mark","tcn.em"]]
        gistic_input2["sample"]=sample_name
        gistic_input2["tcn.em"]=[math.log2(item+0.01) for item in gistic_input2["tcn.em"]]
        gistic_input2=gistic_input2[["sample","chrom","start","end","num.mark","tcn.em"]]
        post_gistic_input2.append(gistic_input2)    
        
pre_gistic_input1_dataframe=pd.concat(pre_gistic_input1)        
pre_gistic_input2_dataframe=pd.concat(pre_gistic_input2)
post_gistic_input1_dataframe=pd.concat(post_gistic_input1)
post_gistic_input2_dataframe=pd.concat(post_gistic_input2)

pre_gistic_input1_dataframe.to_csv("pre_gistic_input1_dataframe.txt",sep="\t")
pre_gistic_input2_dataframe.to_csv("pre_gistic_input2_dataframe.txt",sep="\t")

post_gistic_input1_dataframe.to_csv("post_gistic_input1_dataframe.txt",sep="\t")
post_gistic_input2_dataframe.to_csv("post_gistic_input2_dataframe.txt",sep="\t")
