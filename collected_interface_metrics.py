#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 16:14:17 2022

@author: lhoeg
"""

#Input files: PDB, scores.json (is this the same as what alphafold buts into b_factor column of pdb?), pae.json
#Input parameters: Folder with above input files (F), interface distance (A)
#Alternative for future development, prefix of files to match, since for each input run, AF@ returns files with same prefix

#Output: csv of error metrics for each input model
#Metrics: Mean interface PAE, Mean interface pLDDT, pDockQ, iptm

#PAE should be summarized from both interaction quadrants (bottom-left and top-right)
#Get mean PAE from interacting residues, as in previous script, then take average for whole site of interaction
#Mean interface pLDDT: this is in the scores file, as well as b_factor column of pdb. Using pdb for convenience
#pDockQ is like in the NatureComm paper. Using code as described on gitlab for ElofssonLab
#iptm is that last field in the scores.json file

#Each model input will have one row output with each of these metrics as csv's
#Final output is a tsv or csv file with each row representing a different model input.


import sys
import argparse
import json
from math import sqrt
from statistics import mean
import glob
import re
import numpy as np
from biopandas.pdb import PandasPdb
from scipy.spatial import distance_matrix

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

args = argparser.parse_args(sys.argv[1:])
folder_string = args.folder
distance = args.distance
out_file = args.out_file
out_delim_type=args.out_file_csv


def dict_from_pdbATOM(ATOM_df) :
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
    #pdb_dict['chain2_len'] = dist_c1_c2.shape[0]
    pdb_dict['chain1_plddt'] = chain1_plddt
    pdb_dict['chain2_plddt'] = chain2_plddt
    pdb_dict['chain1_id'] = chain1_id
    pdb_dict['chain2_id'] = chain2_id
    return pdb_dict

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

#Initialize output file
if out_delim_type :
    ods=","
    out_file = out_file+".csv"
else :
    ods="\t"
    out_file = out_file+".txt"


out_open = open(out_file, "w")
init_line="AF_Run"+ods+"Rank"+ods+"Mean_interface_PAE"+ods+"Mean_interface_pLDDT"+ods+"pDockQ"+ods+"iptm\n"
out_open.write(init_line)

#For each prefix in list, get scores for each model. 
#Start with rank_1, since that is what is used in the PAE
for run in prefix_lst :
    model_pdbs = [file for file in pdb_files if run in file]
    model_scores = [file for file in scores_files if run in file]
    for rank in range(1,len(model_pdbs)+1) :
        #Get pdb file & associated data
        current_pdb = [file for file in model_pdbs if "rank_"+str(rank) in file][0]
        fpdb = PandasPdb()
        fpdb.read_pdb(current_pdb)
        fpdb_df=fpdb.df['ATOM']
        current_pdb_dict = dict_from_pdbATOM(fpdb_df)
        #Distance matrix
        dist_mat = current_pdb_dict['dist_mat']
        #Use distances and distance threshold to determine which interactions are proximal
        contacts = np.argwhere(dist_mat<=distance)
        #Use interface to determine if_plddt and pdockQ
        plddt1, plddt2 = np.array(current_pdb_dict['chain1_plddt']), np.array(current_pdb_dict['chain2_plddt'])
        if contacts.shape[0]<1:
            pdockq=0
            ppv=0
            avg_if_plddt=0
        else:
            avg_if_plddt = np.average(np.concatenate([plddt1[np.unique(contacts[:,0])], plddt2[np.unique(contacts[:,1])]]))
            n_if_contacts = contacts.shape[0]
            x = avg_if_plddt*np.log10(n_if_contacts)
            pdockq = 0.724 / (1 + np.exp(-0.052*(x-152.611)))+0.018
        #Get ptm score
        current_scores = [file for file in model_scores if "rank_"+str(rank) in file][0]
        with open(current_scores, 'r') as f:
            scores = json.load(f)
        
        iptm = scores['ptm']

        if rank == 1 :
            current_pae = [file for file in error_files if run in file][0]
            with open(current_pae, 'r') as f:
                pae_all = json.load(f)[0]
            sqr_mat = int(sqrt(len(pae_all['distance'])))
            error = np.reshape(np.array(pae_all['distance']), (sqr_mat,sqr_mat))
            l_c1 = len(plddt1)
            error_TR = error[l_c1:,0:l_c1].transpose()
            error_BL = error[0:l_c1,l_c1:]
            error_mat_mean = (error_TR+error_BL)/2
            if contacts.shape[0]<1:
                avg_if_pae = "NA"
            else:
                pae_keep = []
                for con in contacts :
                    pae_keep.append(error_mat_mean[con[0],con[1]])
                avg_if_pae = mean(pae_keep)
        else:
            avg_if_pae = "NA"
        
        out_line=run+ods+str(rank)+ods+str(avg_if_pae)+ods+str(avg_if_plddt)+ods+str(pdockq)+ods+str(iptm)+"\n"
        out_open.write(out_line)

out_open.close()


