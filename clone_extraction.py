# -*- coding: utf-8 -*-
"""
Created on Wed Jun 10 20:14:39 2020

@author: jxf43
"""
import os 
import pandas as pd
import itertools

tree_dataframe=tree_file
def clone_order(tree_dataframe):
    #parents=[]
    #childs=[]
    a=tree_dataframe["parent"].values.tolist()
    b=tree_dataframe["child"].values.tolist()
    a.extend(b)
    c=set(a)
    #total=[]
    parent_dict={}
    for i in c:
        if not (i in list(tree_dataframe["child"])):
            child=list(tree_dataframe["child"][tree_dataframe["parent"]==i])
            parent_dict[i]=child
    new_tree_dataframe=tree_dataframe[~tree_dataframe["parent"].isin(parent_dict.keys())]
            #parents.append(parent)
            #childs.append(child)
            #total.append(parent_dict)
            #tree_dataframe=new_tree_dataframe
            
    return(parent_dict,new_tree_dataframe) 

  
def clone_order_total(tree_dataframe):
    total=[]
    while tree_dataframe.shape[0]>0:
        result=clone_order(tree_dataframe)[0]
        total.append(result)
        tree_dataframe=clone_order(tree_dataframe)[1]
    order_list=[0]
    clone_list=[[int(k) for k in total[0]][0]]
    #order_dict={} 
    #order_dict[total[0].keys()]=0
    for i in range(len(total)):
        item=total[i]
        childs=set(list(itertools.chain.from_iterable(total[i].values())))
        for item in childs:
            #order_dict[item]=i+1
            order_list.append(i+1)
            clone_list.append(item)
    result=pd.DataFrame({"Clone":clone_list,"order":order_list})        
    return(result) 



total_file=[]
input_dir="I:\\xiaojun\\shiguanai\\pyclone_coding_mutation_20"
total_file=os.listdir(input_dir)
total_file=[item.split(".")[0] for item in total_file]
total_file=list(set(total_file)) 

patient="CXM"

ccf_final=[] 
clone_count_total=[]
clone_number_total=[]
patients=[]
  
for patient in total_file:
    #for file in total_file:
        if not patient.startswith("Rplots"):
            ccf_file=os.path.join(input_dir,patient+".result.csv")
            ccf_file=pd.read_csv(ccf_file)
            #ccf_file.drop(['Chr','Position'],axis=1,inplace=True)
            ccf_file["mutation_id"]=ccf_file["Gene"]+"_"+ccf_file["Chr"].astype(str)+"_"+ccf_file["Position"].astype(str)
            ccf_file.drop(['Chr','Position','Mutation_type','Gene'],axis=1,inplace=True)
            ccf_long=ccf_file.melt(id_vars=['Clone','mutation_id'])
            ccf_long.columns=['Clone','mutation_id', 'Sample_name','CCF_mutation']
            tree_file=os.path.join(input_dir,patient+".GA.consensusTree.new")
            tree_file=pd.read_csv(tree_file,sep="\t",header=None)
            tree_file.columns=["parent","child"]
            tree_file_rename=tree_file.rename(columns={'parent':'Clone'})
            ccf_child=ccf_long.merge(tree_file_rename,how="left")
            tree_file_rename=tree_file.rename(columns={'child':'Clone'})
            ccf_parent=ccf_child.merge(tree_file_rename,how="left")
            order_result=clone_order_total(tree_file)
            ccf_parent_order=ccf_parent.merge(order_result,how="left")
            #ccf_parent_order["Patient_Id"]=[patient[1:]]*ccf_parent_order.shape[0]#用于非独立样本，即病人只有一个树
            ccf_parent_order["Patient"]=patient#用于非独立样本，即病人只有一个树
            #ccf_final.append(ccf_parent_order)
            clone_ccf=os.path.join(input_dir,patient+".cluster.cellularity.csv")
            clone_ccf=pd.read_csv(clone_ccf)
            #clone_ccf["sampleID"]=[item.split("-")[0] for item in clone_ccf["sampleID"]]
            clone_ccf.columns=['Sample_name','Clone', 'CCF_clone','CCF_clone_sd']
            ccf_parent_order=ccf_parent_order.merge(clone_ccf,how="left",on=["Sample_name","Clone"])
            #ccf_parent_order=ccf_parent_order.rename(columns={'Sample_ID':'Sample_ID_correct'})
            ccf_parent_order["tree_id"]=[patient]*ccf_parent_order.shape[0]
            ccf_final.append(ccf_parent_order)
            clone_number=len(set(ccf_parent_order["Clone"]))
            clone_count = ccf_parent_order.groupby('Clone')
            clone_count=clone_count.aggregate({'mutation_id':len})
            clone_count=clone_count.reset_index()
            clone_count.columns=["Clone","count"]
            clone_count["count"]=clone_count["count"]/2
            clone_count["Patient"]=patient
            clone_count_total.append(clone_count)
            clone_number_total.append(clone_number)
            patients.append(patient)
            

ccf_final_all_mutation_5=pd.concat(ccf_final)
ccf_final_all_mutation_5["parameter"]="all_mutation_5"
clone_count_total_all_mutation_5=pd.concat(clone_count_total)
clone_count_total_all_mutation_5["parameter"]="all_mutation_5"
clone_number_all_mutation_5=clone_number_total


ccf_final_all_mutation_6=pd.concat(ccf_final)
ccf_final_all_mutation_6["parameter"]="all_mutation_6"
clone_count_total_all_mutation_6=pd.concat(clone_count_total)
clone_count_total_all_mutation_6["parameter"]="all_mutation_6"
clone_number_all_mutation_6=clone_number_total

ccf_final_all_mutation_7=pd.concat(ccf_final)
ccf_final_all_mutation_7["parameter"]="all_mutation_7"
clone_count_total_all_mutation_7=pd.concat(clone_count_total)
clone_count_total_all_mutation_7["parameter"]="all_mutation_7"
clone_number_all_mutation_7=clone_number_total

ccf_final_all_mutation_8=pd.concat(ccf_final)
ccf_final_all_mutation_8["parameter"]="all_mutation_8"
clone_count_total_all_mutation_8=pd.concat(clone_count_total)
clone_count_total_all_mutation_8["parameter"]="all_mutation_8"
clone_number_all_mutation_8=clone_number_total


ccf_final_all_mutation_4=pd.concat(ccf_final)
ccf_final_all_mutation_4["parameter"]="all_mutation_4"
clone_count_total_all_mutation_4=pd.concat(clone_count_total)
clone_count_total_all_mutation_4["parameter"]="all_mutation_4"
clone_number_all_mutation_4=clone_number_total

# coding mutation
ccf_final_coding_mutation_5=pd.concat(ccf_final)
ccf_final_coding_mutation_5["parameter"]="coding_mutation_5"
clone_count_total_coding_mutation_5=pd.concat(clone_count_total)
clone_count_total_coding_mutation_5["parameter"]="coding_mutation_5"
clone_number_coding_mutation_5=clone_number_total


ccf_final_coding_mutation_6=pd.concat(ccf_final)
ccf_final_coding_mutation_6["parameter"]="coding_mutation_6"
clone_count_total_coding_mutation_6=pd.concat(clone_count_total)
clone_count_total_coding_mutation_6["parameter"]="coding_mutation_6"
clone_number_coding_mutation_6=clone_number_total

ccf_final_coding_mutation_7=pd.concat(ccf_final)
ccf_final_coding_mutation_7["parameter"]="coding_mutation_7"
clone_count_total_coding_mutation_7=pd.concat(clone_count_total)
clone_count_total_coding_mutation_7["parameter"]="coding_mutation_7"
clone_number_coding_mutation_7=clone_number_total

ccf_final_coding_mutation_8=pd.concat(ccf_final)
ccf_final_coding_mutation_8["parameter"]="coding_mutation_8"
clone_count_total_coding_mutation_8=pd.concat(clone_count_total)
clone_count_total_coding_mutation_8["parameter"]="coding_mutation_8"
clone_number_coding_mutation_8=clone_number_total


ccf_final_coding_mutation_4=pd.concat(ccf_final)
ccf_final_coding_mutation_4["parameter"]="coding_mutation_4"
clone_count_total_coding_mutation_4=pd.concat(clone_count_total)
clone_count_total_coding_mutation_4["parameter"]="coding_mutation_4"
clone_number_coding_mutation_4=clone_number_total


ccf_final_coding_mutation_20=pd.concat(ccf_final)
ccf_final_coding_mutation_20["parameter"]="coding_mutation_20"
clone_count_total_coding_mutation_20=pd.concat(clone_count_total)
clone_count_total_coding_mutation_20["parameter"]="coding_mutation_20"
clone_number_coding_mutation_20=clone_number_total

clone_count_total_coding_mutation_20.to_csv("clone_count_total_coding_mutation_20.csv",index=False)


ccf_final_result=[ccf_final_all_mutation_5,ccf_final_all_mutation_6,ccf_final_all_mutation_7,ccf_final_all_mutation_8,ccf_final_coding_mutation_5,ccf_final_coding_mutation_6,ccf_final_coding_mutation_7,ccf_final_coding_mutation_8,ccf_final_coding_mutation_4]
ccf_final_result=pd.concat(ccf_final_result)
ccf_final_result.to_csv("ccf_final_result.csv")

clone_count_final_result=[clone_count_total_all_mutation_5,clone_count_total_all_mutation_6,clone_count_total_all_mutation_7,clone_count_total_all_mutation_8,clone_count_total_coding_mutation_5,clone_count_total_coding_mutation_6,clone_count_total_coding_mutation_7,clone_count_total_coding_mutation_8,clone_count_total_coding_mutation_4]
clone_count_final_result=pd.concat(clone_count_final_result)
clone_count_final_result.to_csv("clone_count_final_result.csv")


pd.set_option('display.max_columns', 500)
pd.set_option('expand_frame_repr', False)

ccf_final_result=pd.read_csv("I:\\xiaojun\\shiguanai\\ccf_final_result.csv")
ccf_final_result.shape

test=ccf_final_result_coding[(ccf_final_result_coding["parameter"]=="coding_mutation_8") & (ccf_final_result_coding["Patient"]=="ACG") & (ccf_final_result["Sample_name"]=="FB201K0346")]

ccf_table=test
clone=7
def information_extract(ccf_table):
    clone_ccf=ccf_table[["Clone","child","parent","CCF_clone"]]
    clone_ccf.drop_duplicates(inplace=True)
    parent_clone=clone_ccf[~pd.isna(clone_ccf["child"])]
    parent_clone=list(set(parent_clone["Clone"]))
    parentclone=[]
    parentccf=[]
    childtotalccf=[]
    for clone in parent_clone:
        parent_ccf=list(clone_ccf[clone_ccf["Clone"]==clone]["CCF_clone"])[0]
        child_clone=list(set(clone_ccf[clone_ccf["Clone"]==clone]["child"]))
        child_total_ccf=0
        for child in child_clone:
            child_ccf=list(clone_ccf[clone_ccf["Clone"]==child]["CCF_clone"])[0]
            child_total_ccf=child_total_ccf+child_ccf
        parentclone.append(clone) 
        parentccf.append(parent_ccf)
        childtotalccf.append(child_total_ccf)
    parent_child_ccf={"parent_clone":parentclone,"parent_ccf":parentccf,"child_total_ccf":childtotalccf}  
    parent_child_ccf=pd.DataFrame(parent_child_ccf)
    parent_child_ccf["parent_child"]=parent_child_ccf["parent_ccf"]-parent_child_ccf["child_total_ccf"]
    parent_child_dif=sum(parent_child_ccf["parent_child"])    
    parent_child_dif_count=sum(parent_child_ccf["parent_child"]<0)
    wrong_ccf=list(parent_child_ccf["parent_child"][parent_child_ccf["parent_child"]<0])
    wrong_ccf=[str(item) for item in wrong_ccf]
    wrong_ccf=",".join(wrong_ccf)
    wrong_ccf_percentage=parent_child_dif_count/parent_child_ccf.shape[0]
    gene=ccf_table["mutation_id"]
    gene=[item.split("_")[0] for item in list(gene)]
    if "TP53" in gene:
        indices = [i for i, x in enumerate(gene) if x == "TP53"]
        TP53_order=min(ccf_table.iloc[indices,]["order"])
    else:
        TP53_order=9999
    result={"parent_child_dif":[parent_child_dif],"parent_child_dif_count":[parent_child_dif_count],"wrong_ccf":wrong_ccf,"wrong_ccf_percentage":[wrong_ccf_percentage],"TP53_order":[TP53_order]}    
    result=pd.DataFrame(result)
    return(parent_child_ccf,result)

ccf_final_result_coding=ccf_final_result[ccf_final_result["parameter"].str.contains("coding") ]    
all_patient= list(set(ccf_final_result["Patient"]))  
all_parameter= list(set(ccf_final_result["parameter"])) 

total_sample_clone_ccf=[]
parent_child_total=[]
for patient in all_patient:
    for parameter in all_parameter:
        target_table=ccf_final_result[(ccf_final_result["parameter"]==parameter) & (ccf_final_result["Patient"]==patient)]
        if target_table.shape[0]>=1:
            for sample in list(set(target_table["Sample_name"])):
                sample_table=target_table[target_table["Sample_name"]==sample]
                sample_result=information_extract(sample_table)
                sample_clone_ccf=sample_result[0]
                parent_child=sample_result[1]
                sample_clone_ccf["patient"]=patient
                sample_clone_ccf["parameter"]=parameter
                sample_clone_ccf["Sample_name"]=sample
                total_sample_clone_ccf.append(sample_clone_ccf)
                parent_child["patient"]=patient
                parent_child["parameter"]=parameter
                parent_child["Sample_name"]=sample
                parent_child_total.append(parent_child)
                
parent_child_total_all=pd.concat(parent_child_total)                
parent_child_total_all.to_csv("parent_child_total_all.csv")
