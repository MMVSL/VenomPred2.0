#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 14 10:19:43 2023

@author: salvo
"""

print("""
.___  ___. .___  ___. ____    ____   _______. __      
|   \/   | |   \/   | \   \  /   /  /       ||  |     
|  \  /  | |  \  /  |  \   \/   /  |   (----`|  |     
|  |\/|  | |  |\/|  |   \      /    \   \    |  |     
|  |  |  | |  |  |  |    \    / .----)   |   |  `----.
|__|  |__| |__|  |__|     \__/  |_______/    |_______|


 _    __                           ____                __   ___    ____     
| |  / /__  ____  ____  ____ ___  / __ \________  ____/ /  |__ \  / __ \    
| | / / _ \/ __ \/ __ \/ __ `__ \/ /_/ / ___/ _ \/ __  /   __/ / / / / /    
| |/ /  __/ / / / /_/ / / / / / / ____/ /  /  __/ /_/ /   / __/_/ /_/ /     
|___/\___/_/ /_/\____/_/ /_/ /_/_/   /_/   \___/\__,_/   /____(_)____/      
                                                                            
                                                                                                  
University of Pisa - Department of Pharmacy
Authors: Salvatore Galati, Miriana Di Stefano      
""")


import os, joblib, PubChemFingerprints, argparse
import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
parser=argparse.ArgumentParser(description="VenomPred2.0 API")
parser.add_argument('-in','--input_file',required=True,type=str,help="Input csv file (It must contain a SMILES fields")
parser.add_argument('-o','--output',required=True, type=str,help="Output file")
parser.add_argument('--endpoint', choices=["Mutagenicity", "Carcinogenicity","Hepatotoxicity", "Estrogenicity", "Androgenicity", "Eye_Irritation",
                                                         "Skin_Irritation", "Acute_Oral_Toxicity", "All"], default="All",help="Endpoints to predict")

args=parser.parse_args()

output_name = args.output
if ".csv" not in output_name: output_name = output_name + ".csv"
    
header = """VenomPred 2.0 is a chemoinformatics platform for toxicity predictions.

An endpoint score is provided for each compound if it was selected in the submission.
The score reports the probability of toxicity calculated by a machine learning consensus approach.  
For full details please read the paper "VenomPred: A Machine Learning Based Platform for Molecular Toxicity Predictions. Int. J. Mol. Sci. 2022;(23)2105. https://doi.org/10.3390/ijms23042105"

"""

#Function to calculate FP and convert it into a numpy array
def fp_as_array(mol, fp_choise):
    if fp_choise == "Morgan":
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)
    elif fp_choise == "RDKit":
        fp = Chem.RDKFingerprint(mol,fpSize=1024)
    elif fp_choise == "PubChem":
        fp = PubChemFingerprints.calcPubChemFingerAll(mol)
        bitstring = "".join([str(x) for x in fp])
        fp = DataStructs.CreateFromBitString(bitstring)
    arr = np.zeros((1,), int)
    DataStructs.ConvertToNumpyArray(fp, arr)
    return arr

def tanimoto_kernel(X, Y):
    xy = np.dot(X,Y.T)
    xx = (X*X).sum(1) #np.dot(X,X.T)
    m = xx.shape[0]
    xx = xx.reshape((m,1))
    yy = (Y*Y).sum(1) #np.dot(Y,Y.T)
    n = yy.shape[0]
    yy = yy.reshape((1,n))
    return xy/(xx+yy-xy)

#Calculate the endpoint single models predictions for a smile         
def get_predictions(smiles,endpoint):
    mol = Chem.MolFromSmiles(smiles)
    models_files = [x for x in os.scandir(f"Endpoints/{endpoint}/Models/") if "dump" in x.name]
    models = {x.name.replace(".dump",""):joblib.load(x) for x in models_files}
    preds,proba = [], []
    for m in models:
        desc_sel = m.split("_")[-1]
        model_ML = models[m]
        fp = fp_as_array(mol, desc_sel) #Calculate FP for each SMILES
        preds.append(model_ML.predict(np.stack(fp).reshape(1,-1))[0])
        prob = model_ML.predict_proba(np.stack(fp).reshape(1,-1))[0][1] #Probability of Toxic prediction
        proba.append(round(prob,2))
    return preds,proba

#Create the report
def make_report(data,selected_endpoints):
    data_csv = []
    for i,smiles in enumerate(data["SMILES"]): #Loop to iterate all the smiles
        name_mol = f"Entry_{i+1}"
        print("[INFO] Predictions for:",name_mol)
        #Chech if the molecule SMILE is correct or it is an empy space
        if Chem.MolFromSmiles(smiles) is None or smiles == "":
            print("[ERROR] in converting to RDKit mol object ")
            data_csv.append([name_mol,smiles] + ["ERROR" for x in selected_endpoints])
            continue
        row_csv = [name_mol,smiles]
        #Loop to predict the endpoints
        for endpoint in selected_endpoints:
            preds, perc_models = get_predictions(smiles,endpoint.replace(" ", "_"))
            row_csv.append(int(round(np.average(perc_models),2)*100))
        data_csv.append(row_csv)
    df = pd.DataFrame.from_records(data_csv, columns=["ID","SMILES"] + [f"Score_{endpoint}" for endpoint in selected_endpoints])
    df.to_csv(output_name,index=False)
    df_content = open(output_name,"r").read()
    df_file_final = open(output_name,"w")
    df_file_final.write(header+df_content)
    df_file_final.close()
    print(f"[INFO] Predictions results saved as {output_name}")

if __name__ == "__main__":
    from warnings import simplefilter
    # ignore all future warnings
    simplefilter(action='ignore', category=FutureWarning)
    try:
        df = pd.read_csv(args.input_file)
    except:
        print("[ERROR] Error opening input file. Make sure the formatting of the file is correct.")
        exit()
        
    if "SMILES" not in df.columns.tolist():
        print("[ERROR] Column SMILES not present in the input file.")
        exit()
        
    if args.endpoint == "All":
        selected_endpoints = ["Mutagenicity", "Carcinogenicity","Hepatotoxicity", "Estrogenicity", "Androgenicity", "Eye_Irritation",
                                                                 "Skin_Irritation", "Acute_Oral_Toxicity"]
    else:
        selected_endpoints = [args.endpoint]
        
    make_report(df,selected_endpoints)
