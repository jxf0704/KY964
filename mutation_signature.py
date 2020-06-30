# -*- coding: utf-8 -*-
"""
Created on Mon May 25 11:30:33 2020

@author: jxf43
"""

#from SigProfilerMatrixGenerator import install as genInstall
#genInstall.install('GRCh37', rsync=False,bash=False)
import os
#os.path.dirname(os.path.abspath(__file__))
#import shutil
#shutil.which("wget") 
#os.environ()
#shutil.which("python")

import pandas as pd

input_dir="I:\\xiaojun\\shiguanai"
patient_information=pd.read_csv(os.path.join(input_dir,"sample_information.csv"))
input_dir="I:\\xiaojun\\shiguanai\\GetBasecounts\\vcf\\input"
files=os.listdir(input_dir)
patientname=[item.split("_")[0] for item in files]
patientname=list(set(patientname))

pre_sample=[item.split("-")[0] for item in patient_information["pre"]]
post_sample=[item.split("-")[0] for item in patient_information["post"]]

#os.mkdir("I:\\xiaojun\\shiguanai\\sigprofiler")
os.makedirs("I:\\xiaojun\\shiguanai\\sigprofiler\\pre_full")
os.makedirs("I:\\xiaojun\\shiguanai\\sigprofiler\\post_full")
os.makedirs("I:\\xiaojun\\shiguanai\\sigprofiler\\common")
os.makedirs("I:\\xiaojun\\shiguanai\\sigprofiler\\pre_specific")
os.makedirs("I:\\xiaojun\\shiguanai\\sigprofiler\\post_specific")
pre_full_dir="I:\\xiaojun\\shiguanai\\sigprofiler\\pre_full"
post_full_dir="I:\\xiaojun\\shiguanai\\sigprofiler\\post_full"
common_dir="I:\\xiaojun\\shiguanai\\sigprofiler\\common"
pre_specific_dir="I:\\xiaojun\\shiguanai\\sigprofiler\\pre_specific"
post_specific_dir="I:\\xiaojun\\shiguanai\\sigprofiler\\post_specific"

patient="anchuanguo"
for patient in patientname: 
    patient_file=[f for f in files if patient in f]
    for name in patient_file:
        if "common" in name:
            common_file=os.path.join(input_dir,name)
            common_table=pd.read_table(common_file,header=None)
            patient_file.remove(name)
    for name in patient_file:
        if "common" in name:
            common_file=os.path.join(input_dir,name)
            common_table=pd.read_table(common_file,header=None)
        else:
            sample_name=name.split("_")[1].split("-")[0]
            if sample_name in pre_sample:
                pre_specific_table=pd.read_table(os.path.join(input_dir,name),header=None)
                pre_full_table= pd.concat([pre_specific_table,common_table], axis = 0, ignore_index = True)
                pre_full_table.iloc[:,2]=pre_full_table.iloc[:,2][0]
            elif sample_name in post_sample:
                post_specific_table=pd.read_table(os.path.join(input_dir,name),header=None)
                post_full_table= pd.concat([post_specific_table,common_table], axis = 0, ignore_index = True)
                post_full_table.iloc[:,2]=post_full_table.iloc[:,2][0]
    common_table.to_csv(os.path.join(common_dir,patient+"_common.vcf"),sep="\t",index=False,header=False) 
    pre_specific_table.to_csv(os.path.join(pre_specific_dir,patient+"_pre_specific.vcf"),sep="\t",index=False,header=False)
    post_specific_table.to_csv(os.path.join(post_specific_dir,patient+"_post_specific.vcf"),sep="\t",index=False,header=False)
    pre_full_table.to_csv(os.path.join(pre_full_dir,patient+"_pre_full.vcf"),sep="\t",index=False,header=False)
    post_full_table.to_csv(os.path.join(post_full_dir,patient+"_post.vcf"),sep="\t",index=False,header=False)

import os       
os.chdir("I:\\xiaojun\\shiguanai")
from sigproextractor import sigpro as sig
data = "I:/xiaojun/shiguanai/sigprofiler/post_full"
sig.sigProfilerExtractor("vcf", "post_full_output", data, startProcess=1, endProcess=10, genome_build = "GRCh37", refgen="GRCh37",mtype = ["96"])

data = "I:/xiaojun/shiguanai/sigprofiler/pre_specific"
sig.sigProfilerExtractor("vcf", "pre_specific_output", data, startProcess=1, endProcess=10, genome_build = "GRCh37", refgen="GRCh37",mtype = ["96"])

data = "I:/xiaojun/shiguanai/sigprofiler/common"
sig.sigProfilerExtractor("vcf", "common_output", data, startProcess=1, endProcess=10, genome_build = "GRCh37", refgen="GRCh37",mtype = ["96"])

data = "I:/xiaojun/shiguanai/sigprofiler/post_specific"
sig.sigProfilerExtractor("vcf", "post_specific_output", data, startProcess=1, endProcess=10, genome_build = "GRCh37", refgen="GRCh37",mtype = ["96"])