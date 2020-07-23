#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sys
import glob
import pandas as pd
import numpy as np
from collections import defaultdict



# ### Opening and storing the Xist sequence:

# In[2]:


sequencefile = open("../Reference_sequences/Xist.seq","r").readlines()
sequence = "".join([i.strip() for i in sequencefile])


# In[3]:


sequence[:10]


# ### Adding secondary structure information:

# In[4]:


secstructure = pd.read_csv("../Data/Xist_secondary_structure_info.txt",sep="\t",\
                          header=0,\
                          index_col=0,\
                          names=['nucleotide','base_paired_to'])



# In[5]:


secstructure.head()


# ### diffBUM_HMM data made by Toby on the 27th of May 2020:
# Two datasets, one normalized by denatured, the other without normalization.

# In[6]:


diffBUMHMMnorm = pd.read_csv("../Analysis/diffBUM-HMM/Xist_diff_BUM_HMM_analysis_withscaling_empty_first_line_last_line_removed.txt",\
                         sep="\t",\
                         header=0,\
                         index_col=0,\
                         names=["nucleotide","scaled_unmodified","scaled_ex_vivo","scaled_in_vivo","scaled_both"])

diffBUMHMM = pd.read_csv("../Analysis//diffBUM-HMM/Xist_diff_BUM_HMM_analysis_noscaling_empty_first_line_last_line_removed.txt",\
                         sep="\t",\
                         header=0,\
                         index_col=0,\
                         names=["nucleotide","unmodified","ex_vivo","in_vivo","both"])
# In[7]:


diffBUMHMMnorm.head()


# In[8]:


diffBUMHMM.head()


# ### Loading the dStruct analysis results:

# In[9]:

#WORK IN PROGRESS - CONTINUE FROM HERE

dStruct = pd.read_csv("../../../Data/Xist_dStruct_normal_vs_scaled_270520.txt",                        sep="\t",                        index_col=0,                        header=0)


# In[10]:


dStruct.head()


# ### Loading the SHAPE .map files for each replicate:

# In[11]:


ex_vivo_1 = pd.read_csv("../../../Data/Map_files_Xist/XIST_exvivo_1M7_rep1.map",                        sep="\t",                        index_col=0,                        header=None,                        names=["SHAPE_reactivity","Error","nucleotide"])
ex_vivo_2 = pd.read_csv("../../../Data/Map_files_Xist/XIST_exvivo_1M7_rep2.map",                        sep="\t",                        index_col=0,                        header=None,                        names=["SHAPE_reactivity","Error","nucleotide"])
in_cell_1 = pd.read_csv("../../../Data/Map_files_Xist/XIST_incell_1M7_rep1.map",                        sep="\t",                        index_col=0,                        header=None,                        names=["SHAPE_reactivity","Error","nucleotide"])
in_cell_2 = pd.read_csv("../../../Data/Map_files_Xist/XIST_incell_1M7_rep2.map",                        sep="\t",                        index_col=0,                        header=None,                        names=["SHAPE_reactivity","Error","nucleotide"])


# In[12]:


ex_vivo_1.head()


# ### Loading the RBP binding sites:

# In[13]:


proteins = ["CELF1","FUS","HuR","PTBP1","RBFOX2","TARDBP"]


# In[14]:


directory = "../../../Data/Bed_files_Xist/CLIPdb_transcript_coords"
proteindata = defaultdict(list)
for protein in proteins:
    ### First grabbing all the files with protein binding sites:
    files = glob.glob("%s/%s*" % (directory,protein))
    ### Now extracting the start and end positions and storing these
    ### in a list of tuples:
    for i in files:
        data = pd.read_csv(i,                           sep="\t",                           header=None,                           index_col=None)
        ### start is column 1 and end is column 2
        start = data[data.columns[1]].values
        end = data[data.columns[2]].values
        coordinates = (list(zip(start,end)))
        proteindata[protein].extend(coordinates)


# In[15]:


proteindata['CELF1']


# ### Loading the deltaSHAPE reactivities for each replicate:

# In[16]:


deltaSHAPE_1 = pd.read_csv("../../../Data/deltaSHAPE_confirmed_sites/deltaSHAPE_rep1.igv",                           sep="\t",                           header=0,
                           index_col=1,
                           names=["Chromosome","end","feature","deltaSHAPE"])
deltaSHAPE_2 = pd.read_csv("../../../Data/deltaSHAPE_confirmed_sites/deltaSHAPE_rep2.igv",                           sep="\t",                           header=0,
                           index_col=1,
                           names=["Chromosome","end","feature","deltaSHAPE"])


# ### Checking the number of differentially expressed nucleotides:

# In[17]:


len(deltaSHAPE_1[deltaSHAPE_1.deltaSHAPE > 0].index)


# In[18]:


len(deltaSHAPE_1[deltaSHAPE_1.deltaSHAPE < 0].index)


# In[19]:


len(deltaSHAPE_2[deltaSHAPE_2.deltaSHAPE > 0].index)


# In[20]:


len(deltaSHAPE_2[deltaSHAPE_2.deltaSHAPE < 0].index)


# ### Making the new DataFrame:

# In[21]:


columns = ["nucleotide"]
nucleotide_positions = np.arange(1,len(sequence)+1)
samples = ["SHAPE_reactivity_ex_vivo_1",           "SHAPE_reactivity_ex_vivo_2",           "SHAPE_reactivity_in_cell_1",           "SHAPE_reactivity_in_cell_2",           "deltaSHAPE_rep1",           "deltaSHAPE_rep2",           "dStruct",           "dStruct_scaled"]
columns.extend(proteins)
columns.extend(samples)
table = pd.DataFrame(0,index=nucleotide_positions,columns=columns)


# In[22]:


table.head()


# ### Adding the Xist sequence to the table:

# In[23]:


table.nucleotide = list(sequence)


# ### Adding the RBP-binding sites to the table:

# In[24]:


for protein in proteins:
    locations = proteindata[protein]
    for start,end in locations:
        table.loc[start:end+1,protein] = 1


# ### Adding the SHAPE reactivities from each replicate to the table:

# ### Check if the number of nucleotides in the data files is the same:
# If not then information for some nucleotides may be missing in the final table.

# In[25]:


assert len(table.index) == len(ex_vivo_1.index)
assert len(table.index) == len(ex_vivo_2.index)
assert len(table.index) == len(in_cell_1.index)
assert len(table.index) == len(in_cell_2.index)


# In[26]:


table.SHAPE_reactivity_ex_vivo_1 = ex_vivo_1.loc[table.index,"SHAPE_reactivity"]
table.SHAPE_reactivity_ex_vivo_2 = ex_vivo_2.loc[table.index,"SHAPE_reactivity"]
table.SHAPE_reactivity_in_cell_1 = in_cell_1.loc[table.index,"SHAPE_reactivity"]
table.SHAPE_reactivity_in_cell_2 = in_cell_2.loc[table.index,"SHAPE_reactivity"]


# In[27]:


table.head()


# ### Adding the deltaSHAPE values:

# In[28]:


positions = deltaSHAPE_1.index
table.loc[positions,"deltaSHAPE_rep1"] = deltaSHAPE_1.loc[positions,"deltaSHAPE"].values


# In[29]:


positions = deltaSHAPE_2.index
table.loc[positions,"deltaSHAPE_rep2"] = deltaSHAPE_2.loc[positions,"deltaSHAPE"].values


# In[30]:


table.head()


# ### Adding the dStruct data:

# In[31]:


table.dStruct = dStruct.dStruct
table.dStruct_scaled = dStruct.dStruct_scaled


# ### Some positions had no values. Replace with 0:

# In[32]:


len(table[table.deltaSHAPE_rep2 < 0].index)


# In[33]:


len(table[table.deltaSHAPE_rep2 > 0].index)


# In[34]:


table.replace(np.nan,0,inplace=True)


# ### Merging everyting together:

# In[46]:


table[diffBUMHMM.columns] = diffBUMHMM
table[diffBUMHMMnorm.columns] = diffBUMHMMnorm


# In[47]:


table.tail()


# ### Adding secondary structure information:

# In[48]:


table[secstructure.columns] = secstructure


# In[50]:


table.head()


# ### Writing output file:

# In[51]:


table.to_csv("../../../New_data_table.txt",sep="\t")


# ### How many positions do the first and second replicate of the Xist data have in common?

# In[52]:


len(table[(table.deltaSHAPE_rep1 != 0) & (table.deltaSHAPE_rep2 != 0)].index)


# In[ ]:

'''




