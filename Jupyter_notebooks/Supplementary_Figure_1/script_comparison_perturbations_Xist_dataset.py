import subprocess
#retcode = subprocess.call("Rscript --vanilla /Users/maran/Desktop/diff_BUM_HMM_Project/Github/diff_BUM_HMM/test.R", shell=True)
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
#import sys

original_df = pd.read_csv('../../Analysis/diffBUM-HMM/Xist_diff_BUM_HMM_analysis_withscaling_gaussiannoiseadded.txt', sep="\t", header=0)
original_df= original_df.rename_axis('nt').reset_index()

original_df['filteredUM'] = [0 if x<0.90 else x for x in original_df['UM']]
original_df['filteredMU'] = [0 if x<0.90 else x for x in original_df['MU']]
df= original_df[['filteredUM','filteredMU']]

'''
for i in range(0, len(df.index)):
    for j in range(0, len(df.columns)):
        if df.values[i,j] > 0:
            print(df.values[i,j])
            print('row ', i)
            print('column ', j)
'''
top_n = 1
A = pd.DataFrame({n: df.T[col].nlargest(top_n).index.tolist() for n, col in enumerate(df.T)}).T



#print(A)
#print(A.loc[16850:16890,:])

list_of_ratios = []
iterations_counter =0

for iter in range(150):
    print('Iteration: ' +str(iter))
    subprocess.call("Rscript --vanilla ../../Scripts/diffBUM_HMM_Xist_dataset_mismatch_rate.R", shell=True)
        
    temp_df = pd.read_csv('../../Analysis/diffBUM-HMM/Xist_diff_BUM_HMM_analysis_withscaling_gaussiannoiseadded.txt', sep="\t", header=0)
    
    
    
    
    if iter>0:
        print('This iter is larger than 0')
        check = temp_df.equals(ending_df)
        if check is True:
            print('Skipped this iteration')
            continue
        temporary_df=temp_df.copy()
        temporary_df= temporary_df.rename_axis('nt').reset_index()
        temporary_df['filteredUM'] = [0 if x<0.90 else x for x in temporary_df['UM']]
        temporary_df['filteredMU'] = [0 if x<0.90 else x for x in temporary_df['MU']]

        new_df= temporary_df[['filteredUM','filteredMU']]

        B = pd.DataFrame({n: new_df.T[col].nlargest(top_n).index.tolist() for n, col in enumerate(new_df.T)}).T

        counter=0
        other_counter=0
        for i in range(len(A.index)):
            if A.iloc[i:i+1,0:1].values[0][0] == B.iloc[i:i+1,0:1].values[0][0]:
                counter=counter+1
            else:
                other_counter=other_counter+1
 
    
        print(counter)
        print(other_counter)

        ratio = other_counter/counter
        print(ratio)
        list_of_ratios.append(ratio)
        del ratio
        del counter
        del other_counter
        del check
        del new_df

    ending_df=temp_df.copy()
    del temp_df
   
        #del B
        #del temporary_df
    
    
    iterations_counter=iterations_counter+1 
    










list_of_ratios_np = np.array(list_of_ratios).astype(np.float)
print('Mean is ' + str(np.mean(list_of_ratios_np)))
print('SD is ' + str(np.std(list_of_ratios_np)))
print('# of actual iterations' +str(iterations_counter))


#sys.stdout = open('file.txt', 'w')
