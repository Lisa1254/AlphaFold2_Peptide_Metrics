#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 16:14:17 2022

@author: lhoeg

Required input: path to folder with AlphaFold2 outputs, including pdb and scores file for each model, and a predicted_aligned_error file for each AF2 run. 

Output will be either a comma or tab delimited file with Run, Rank, Mean Interface PAE, Mean Interface pLDDT, pDockQ, and iptm for each model.

Default minimum distance for interaction is 9 Angstroms.

Use -h to access help and more options.

See https://www.nature.com/articles/s41467-022-28865-w and https://gitlab.com/ElofssonLab/FoldDock for more information about pDockQ.

"""


import sys
import argparse
import json
from statistics import mean
import glob
import re
import numpy as np
from biopandas.pdb import PandasPdb
from scipy.spatial import distance_matrix
from datetime import datetime

argparser = argparse.ArgumentParser(
        description="Takes folder of AlphaFold output, and returns csv of metrics on interface between chains")

argparser.add_argument("-F", "--folder", action="store", default=None, type=str,
    help="Provide path to folder with results for analysis.")
argparser.add_argument("-A", "--distance", action="store", default=9, type=int,
    help="Provide maximum distance (angstrom) for considering proximity between residues as interaction. Default=9.")
argparser.add_argument("-O", "--out_file", action="store", default="af2_interface_metrics", type=str,
    help="Optional, provide name for output file. Default is 'af2_interface_metrics'. Output file will be .txt if selecting tsv output filetype (default), or .csv if using csv as output filetype.")
argparser.add_argument("-C", "--out_file_csv", action="store_true",
    help="If selected, output will be 'csv' filetype. Otherwise 'tsv' filetype will be used")
argparser.add_argument("-L", "--logfile", action="store_true",
    help="If selected, a logfile will be produced detailing which models have been run")
argparser.add_argument("-V", "--verbose", action="store_true",
    help="If selected, current run being analysed will be printed to the terminal.")

args = argparser.parse_args(sys.argv[1:])
folder_string = args.folder
distance = args.distance
out_file = args.out_file
out_delim_type=args.out_file_csv
logfile = args.logfile
terminal_v = args.verbose


def dict_from_pdbATOM(pdb) :
    fpdb = PandasPdb()
    fpdb.read_pdb(pdb)
    ATOM_df=fpdb.df['ATOM']
    
    pdb_dict = {}
    chain_letters=ATOM_df.chain_id.unique()
    
    chain1_id = chain_letters[0]
    chain2_id = chain_letters[1]
    
    chain1_df = ATOM_df.query('chain_id == @chain1_id')
    chain2_df = ATOM_df.query('chain_id == @chain2_id')
    
    chain1_cb = chain1_df.query('(atom_name == "CB" & residue_name != "GLY") | (atom_name == "CA" & residue_name == "GLY")')
    chain2_cb = chain2_df.query('(atom_name == "CB" & residue_name != "GLY") | (atom_name == "CA" & residue_name == "GLY")')
    
    chain1_coords = chain1_cb[["x_coord", "y_coord", "z_coord"]]
    chain2_coords = chain2_cb[["x_coord", "y_coord", "z_coord"]]
    dist_c1_c2 = distance_matrix(chain1_coords, chain2_coords)
    chain1_plddt = chain1_cb[['b_factor']]
    chain2_plddt = chain2_cb[['b_factor']]
    pdb_dict['dist_mat'] = dist_c1_c2
    pdb_dict['chain1_plddt'] = chain1_plddt
    pdb_dict['chain2_plddt'] = chain2_plddt
    pdb_dict['chain1_id'] = chain1_id
    pdb_dict['chain2_id'] = chain2_id
    return pdb_dict

def mean_if_pae(len_ch1, paes, interfaces) :
    error = np.array(paes)
    rounded_errors = np.round(error.astype(np.float64), decimals=1)
    error_TR = rounded_errors[len_ch1:,0:len_ch1].transpose()
    error_BL = rounded_errors[0:len_ch1,len_ch1:]
    error_mat_mean = (error_TR+error_BL)/2
    if contacts.shape[0]<1:
        avg_if_pae = "NA"
    else:
        pae_keep = []
        for con in interfaces :
            pae_keep.append(error_mat_mean[con[0],con[1]])
        avg_if_pae = mean(pae_keep)
    return avg_if_pae

def calc_pdockq(mean_plddt_interfaces, interfaces):
    n_if_contacts = interfaces.shape[0]
    x = mean_plddt_interfaces*np.log10(n_if_contacts)
    pdockq = 0.724 / (1 + np.exp(-0.052*(x-152.611)))+0.018
    return pdockq

#Specify all files in provided folder as search parameters
if folder_string.endswith("/") :
    files_string = folder_string+"*"
else: 
    files_string = folder_string+"/*"

#Get files in folder
all_files = glob.glob(files_string)

#Separate files by type
pdb_files = [file for file in all_files if ".pdb" in file]
scores_files = [file for file in all_files if "scores.json" in file]
error_files = [file for file in all_files if "predicted_aligned_error" in file]

#get prefix for interaction. Should only be 1 PAE file per run, so take from prefix to the error_files
p = re.compile(r'_predicted_aligned_error_v1.json')
prefix_lst = [p.sub('', file) for file in error_files]
p = re.compile(files_string)
prefix_lst = [p.sub('', file) for file in prefix_lst]
prefix_lst.sort()
n_runs = len(prefix_lst)

#Initialize output file(s)
if out_delim_type :
    ods=","
    out_file = out_file+".csv"
else :
    ods="\t"
    out_file = out_file+".txt"

out_open = open(out_file, "w")
init_line="AF_Run"+ods+"Rank"+ods+"Mean_interface_PAE"+ods+"Mean_interface_pLDDT"+ods+"pDockQ"+ods+"iptm\n"
out_open.write(init_line)

now = datetime.now()
date_time = now.strftime("%Y/%m/%d, %H:%M:%S")
if logfile:
    log_outfile = "log_AF2Metrics_"+now.strftime("%m")+"_"+now.strftime("%d")+"_"+now.strftime("%H")+"_"+now.strftime("%M")+".txt"
    log_open = open(log_outfile, 'w')
    log_open.write("Running collected_interface_metrics.py script, using the following parameters:\n-F\t"+folder_string+"\n-A\t"+str(distance)+"\n-O\t"+out_file+"\n-C\t"+str(out_delim_type)+"\n\nStarting AF2 chain interface analysis at "+date_time+" for "+str(n_runs)+" detected runs.\n\n")

if terminal_v:
    print("Starting AF2 chain interface analysis at "+date_time+" for "+str(n_runs)+" detected runs.\n")

#For each prefix in list, get scores for each model. 
#Start with rank_1 for consistency
for run in prefix_lst :
    if logfile:
        log_open.write("Starting chain interface analysis for "+run+"\n")
    if terminal_v:
        print("Starting chain interface analysis for "+run+"\n")
    
    model_pdbs = [file for file in pdb_files if run in file]
    model_scores = [file for file in scores_files if run in file]
    for rank in range(1,len(model_pdbs)+1) :
        
        #Get pdb file & associated data
        current_pdb = [file for file in model_pdbs if "rank_"+str(rank) in file][0]
        current_pdb_dict = dict_from_pdbATOM(current_pdb)
        
        #Distance matrix
        dist_mat = current_pdb_dict['dist_mat']
        
        #Use distances and distance threshold to determine which interactions are proximal
        contacts = np.argwhere(dist_mat<=distance)
        
        #Use interface to determine interface plddt
        plddt1, plddt2 = np.array(current_pdb_dict['chain1_plddt']), np.array(current_pdb_dict['chain2_plddt'])
        if contacts.shape[0]<1:
            avg_if_plddt=0
            pdockq=0
        else:
            avg_if_plddt = np.average(np.concatenate([plddt1[np.unique(contacts[:,0])], plddt2[np.unique(contacts[:,1])]]))
            pdockq = calc_pdockq(avg_if_plddt, contacts)
        
        #Get iptm score
        current_scores = [file for file in model_scores if "rank_"+str(rank) in file][0]
        with open(current_scores, 'r') as f:
            scores = json.load(f)
        
        #Get mean interface PAE
        l_c1 = len(plddt1)
        avg_if_pae = mean_if_pae(l_c1, scores['pae'], contacts)
        
        #Get iptm
        iptm = scores['iptm']
        
        
        out_line=run+ods+str(rank)+ods+str(avg_if_pae)+ods+str(avg_if_plddt)+ods+str(pdockq)+ods+str(iptm)+"\n"
        out_open.write(out_line)

out_open.close()

now1 = datetime.now()
date_time1 = now1.strftime("%Y/%m/%d, %H:%M:%S")
if logfile:
    log_open.write("\nFinished at "+date_time1+".\n")
    log_open.close()

if terminal_v:
    print("Finished at "+date_time1+".\n")



