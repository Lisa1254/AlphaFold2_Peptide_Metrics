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
#Mean interface pLDDT: this is in the scores file, but might also be found in the pdb??
#pDockQ is like in the NatureComm paper. I've saved their script as pdockq_ElofssonLab.py to see what they do
#iptm is that last field in the scores.json file

#Each model input will have one row output with each of these metrics as csv's
#Final output is a csv file with each row representing a different model input.


import sys
#import os
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
argparser.add_argument("-R", "-reverse_chains", action="store_true",
    help="Reverse base chain and comparison chain.")

args = argparser.parse_args(sys.argv[1:])
folder_string = args.folder
distance = args.distance
reverse_chains = args.reverse_chains

#Temp stuff
#Folder used for testing interactively (single run with 5 models): 
#folder_string = "/Users/lhoeg/Documents/Peptide_Interface/AF_result_GFLHVGG_c_zer1"
#folder_string = "/Users/localadmin/Documents/Reference/Example_AF2_Output"
#Folder used for testing interactively (two runs with 5 models each): 
#folder_string = "/Users/lhoeg/Documents/Peptide_Interface/Results_2runs/"
#reverse_chains=False
#distance = 9

#Function to make distance matrix
#Default treats shorter of chains as "Chain2", add argument later for reversing chains if flag is True
#Consider adding log to specify which chain is which
def dist_mat_from_pdbATOM(ATOM_df, Rev_Bool=False) :
    chain_letters=ATOM_df.chain_id.unique()
    
    top_chain_len=len(ATOM_df.query('chain_id == @chain_letters[0]'))
    bottom_chain_len=len(ATOM_df.query('chain_id == @chain_letters[1]'))
    
    if (top_chain_len > bottom_chain_len and not Rev_Bool) or (top_chain_len < bottom_chain_len and Rev_Bool):
        chain1_id = chain_letters[0]
        chain2_id = chain_letters[1]
    else:
        chain1_id = chain_letters[1]
        chain2_id = chain_letters[0]
    
    chain1_df = ATOM_df.query('chain_id == @chain1_id')
    chain2_df = ATOM_df.query('chain_id == @chain2_id')
    
    chain1_cb = chain1_df.query('(atom_name == "CB" & residue_name != "GLY") | (atom_name == "CA" & residue_name == "GLY")')
    chain2_cb = chain2_df.query('(atom_name == "CB" & residue_name != "GLY") | (atom_name == "CA" & residue_name == "GLY")')
    
    chain1_coords = chain1_cb[["x_coord", "y_coord", "z_coord"]]
    chain2_coords = chain2_cb[["x_coord", "y_coord", "z_coord"]]
    dist_c1_c2 = distance_matrix(chain1_coords, chain2_coords).transpose() #(len(chain2), len(chain1))
    return dist_c1_c2

def dict_from_pdbATOM(ATOM_df, Rev_Bool=False) :
    pdb_dict = {}

    chain_letters=ATOM_df.chain_id.unique()
    
    top_chain_len=len(ATOM_df.query('chain_id == @chain_letters[0]'))
    bottom_chain_len=len(ATOM_df.query('chain_id == @chain_letters[1]'))
    
    if (top_chain_len > bottom_chain_len and not Rev_Bool) or (top_chain_len < bottom_chain_len and Rev_Bool):
        chain1_id = chain_letters[0]
        chain2_id = chain_letters[1]
        pdb_dict['chain2'] = 1
    else:
        chain1_id = chain_letters[1]
        chain2_id = chain_letters[0]
        pdb_dict['chain2'] = 0
    
    chain1_df = ATOM_df.query('chain_id == @chain1_id')
    chain2_df = ATOM_df.query('chain_id == @chain2_id')
    
    chain1_cb = chain1_df.query('(atom_name == "CB" & residue_name != "GLY") | (atom_name == "CA" & residue_name == "GLY")')
    chain2_cb = chain2_df.query('(atom_name == "CB" & residue_name != "GLY") | (atom_name == "CA" & residue_name == "GLY")')
    
    chain1_coords = chain1_cb[["x_coord", "y_coord", "z_coord"]]
    chain2_coords = chain2_cb[["x_coord", "y_coord", "z_coord"]]
    dist_c1_c2 = distance_matrix(chain1_coords, chain2_coords).transpose() #(len(chain2), len(chain1))
    pdb_dict['dist_mat'] = dist_c1_c2
    pdb_dict['chain2_len'] = dist_c1_c2.shape[0]
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

#For each prefix in list, get scores for each model. 
#Start with rank_1, since that is what is used in the PAE
for run in prefix_lst :
    model_pdbs = [file for file in pdb_files if run in file]
    model_scores = [file for file in scores_files if run in file]
    for rank in range(1,len(model_pdbs)+1) :
        #Get pdb file
        current_pdb = [file for file in model_pdbs if "rank_"+str(rank) in file][0]
        fpdb = PandasPdb()
        fpdb.read_pdb(current_pdb)
        fpdb_df=fpdb.df['ATOM']
        current_pdb_dict = dict_from_pdbATOM(fpdb_df, Rev_Bool=reverse_chains)
        #Distance matrix
        #dist_mat = dist_mat_from_pdbATOM(ATOM_df=fpdb_df, Rev_Bool=reverse_chains)
        dist_mat = current_pdb_dict['dist_mat']
        #Use distances and distance threshold to determine which interactions are proximal
        prox_dict = {}
        prox_dict['all'] = set()
        for row in range(0,dist_mat.shape[0]):
            if np.any(dist_mat[row] < distance):
                col_ind_prox = np.where(dist_mat[row] < distance)[0]
                prox_dict[row] = col_ind_prox
                prox_dict['all'] = set(list(prox_dict['all'])+list(col_ind_prox))
        #Get scores
        current_scores = [file for file in model_scores if "rank_"+str(rank) in file][0]
        with open(current_scores, 'r') as f:
            scores = json.load(f)
        
        lddt_scores = np.array(scores['plddt'])
        prox_scores = []
        if pdb_dict['chain2'] == 0 :
            chain2_w_prox = list(prox_dict.keys())[1:]
            prox_scores+=list(lddt_scores[chain2_w_prox])
            chain1_w_prox = list(prox_dict['all'])
            chain1_prox_ind = [x+current_pdb_dict['chain2_len'] for x in chain1_w_prox]
            prox_scores+=list(lddt_scores[chain1_prox_ind])
        else:
            chain2_w_prox = list(prox_dict.keys())[1:]
            chain2_prox_ind = [x+dist_mat.shape[1] for x in chain2_w_prox]
            prox_scores+=list(lddt_scores[chain2_prox_ind])
            chain1_w_prox = list(prox_dict['all'])
            prox_scores+=list(lddt_scores[chain1_w_prox])
        
        mean_prox_lddt = mean(prox_scores)
        iptm = scores['ptm']

        if rank == 1 :
            current_pae = [file for file in error_files if run in file][0]
            with open(current_pae, 'r') as f:
                pae_all = json.load(f)[0]
            sqr_mat = int(sqrt(len(pae_all['distance'])))
            error = np.reshape(np.array(pae_all['distance']), (sqr_mat,sqr_mat))
            #Current position - HERE IN PROGRESS
        else:
            mean_prox_pae = "NA"




#Trying to figure out from colabfold and alphafold where the distance value in their pae file came from.
#Since there is a residue1, residue2, and distance, each with enough values for the chain1:chain2 x chain1:chain2 matrix
#The scores file also has a 'pae' key for the dictionary, which has length chain1:chain2, and each item has chain1:chain2 values, 
#So maybe this relates to the res1 or res2 or distance somehow? Maybe I can calculate for each of the models where it is not provided?
