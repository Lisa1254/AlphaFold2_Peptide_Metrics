#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  1 16:14:17 2022

@author: lhoeg
"""

#Input files: PDB, scores.json (is this the same as what alphafold buts into b_factor column of pdb?), pae.json
#Input parameters: interface distance (A)

#Output: csv of error metrics for each input model
#Metrics: Mean interface PAE, Mean interface pLDDT, pDockQ, iptm

#PAE should be summarized from both interaction quadrants (bottom-left and top-right)
#Get mean PAE from interacting residues, as in previous script, then take average for whole site of interaction
#Mean interface pLDDT: this is in the scores file, but might also be found in the pdb??
#pDockQ is like in the NatureComm paper. I've saved their script as pdockq_ElofssonLab.py to see what they do
#iptm is that last field in the scores.json file

#Each model input will have one row output with each of these metrics as csv's
#Final output is a csv file with each row representing a different model input.


