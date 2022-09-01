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


import sys
import os
import argparse

argparser = argparse.ArgumentParser(
        description="Takes folder of AlphaFold output, and returns csv of metrics on interface between chains")

argparser.add_argument("-F", "--folder", action="store", default=None, type=str,
    help="Provide path to folder with results for analysis.")
argparser.add_argument("-A", "--distance", action="store", default=9, type=int,
    help="Provide maximum distance (angstrom) for considering proximity between residues as interaction. Default=9.")

args = argparser.parse_args(sys.argv[1:])

#Temp stuff
folder_string = args.folder
print(folder_string)

distance = args.distance
print(distance)

