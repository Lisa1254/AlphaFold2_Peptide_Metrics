#Temp file of pdockQ stuff to add to main script
#From Eloffson Lab
#The distance and plddt calculations seem much easier than the ones that I made! Maybe I can use it.

chain_coords, chain_plddt = {}, {}
for line in range(0,fpdb_df.shape[0]):
     if fpdb_df.iloc[line]['atom_name']=='CB' or (fpdb_df.iloc[line]['atom_name']=='CA' and fpdb_df.iloc[line]['residue_name']=='GLY'):
        if fpdb_df.iloc[line]['chain_id'] in [*chain_coords.keys()]:
            chain_coords[fpdb_df.iloc[line]['chain_id']].append([fpdb_df.iloc[line]['x_coord'],fpdb_df.iloc[line]['y_coord'],fpdb_df.iloc[line]['z_coord']])
            chain_plddt[fpdb_df.iloc[line]['chain_id']].append(fpdb_df.iloc[line]['b_factor'])
        else:
            chain_coords[fpdb_df.iloc[line]['chain_id']] = [[fpdb_df.iloc[line]['x_coord'],fpdb_df.iloc[line]['y_coord'],fpdb_df.iloc[line]['z_coord']]]
            chain_plddt[fpdb_df.iloc[line]['chain_id']] = [fpdb_df.iloc[line]['b_factor']]


for chain in chain_coords:
    chain_coords[chain] = np.array(chain_coords[chain])
    chain_plddt[chain] = np.array(chain_plddt[chain])

ch1, ch2 = [*chain_coords.keys()]
coords1, coords2 = chain_coords[ch1], chain_coords[ch2]
plddt1, plddt2 = chain_plddt[ch1], chain_plddt[ch2]

mat = np.append(coords1, coords2,axis=0)
a_min_b = mat[:,np.newaxis,:] -mat[np.newaxis,:,:]
#Note: the transpose argument in this command makes the resulting array the transpose of the way that I've been using it in my scripts. If I omit the transpose here, then be sure to swap the order of plddt1 & plddt2 in the avg_if_plddt command below.
dists = np.sqrt(np.sum(a_min_b.T ** 2, axis=0)).T

l1 = len(coords1)
contact_dists = dists[:l1,l1:]
#Here manually input 9, fix for variable
contacts = np.argwhere(contact_dists<=9)

if contacts.shape[0]<1:
    pdockq=0
    ppv=0 #Might end up removing PPV Stuff
else:
    avg_if_plddt = np.average(np.concatenate([plddt1[np.unique(contacts[:,0])], plddt2[np.unique(contacts[:,1])]]))
    #avg_if_plddt = np.average(np.concatenate([plddt2[np.unique(contacts[:,0])], plddt1[np.unique(contacts[:,1])]])) #for if not transposing array above
    n_if_contacts = contacts.shape[0]
    x = avg_if_plddt*np.log10(n_if_contacts)
    pdockq = 0.724 / (1 + np.exp(-0.052*(x-152.611)))+0.018

#PPV Stuff from here, might remove
PPV = np.array([0.98128027, 0.96322524, 0.95333044, 0.9400192 ,0.93172991, 0.92420274, 0.91629946, 0.90952562, 0.90043139,0.8919553 , 0.88570037, 0.87822061, 0.87116417, 0.86040801,0.85453785, 0.84294946, 0.83367787, 0.82238224, 0.81190228,0.80223507, 0.78549007, 0.77766077, 0.75941223, 0.74006263,0.73044282, 0.71391784, 0.70615739, 0.68635536, 0.66728511,0.63555449, 0.55890174])
pdockq_thresholds = np.array([0.67333079, 0.65666073, 0.63254566, 0.62604391,0.60150931, 0.58313803, 0.5647381 , 0.54122438, 0.52314392,0.49659878, 0.4774676 , 0.44661346, 0.42628389, 0.39990988,0.38479715, 0.3649393 , 0.34526004, 0.3262589 , 0.31475668,0.29750023, 0.26673725, 0.24561247, 0.21882689, 0.19651314,0.17606258, 0.15398168, 0.13927677, 0.12024131, 0.09996019,0.06968505, 0.02946438])
inds = np.argwhere(pdockq_thresholds>=pdockq)

if len(inds)>0:
    ppv = PPV[inds[-1]][0]
else:
    ppv = PPV[0]