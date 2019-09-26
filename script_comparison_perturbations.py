import subprocess
#retcode = subprocess.call("Rscript --vanilla /Users/maran/Desktop/diff_BUM_HMM_Project/Github/diff_BUM_HMM/test.R", shell=True)
import matplotlib.pyplot as plt
import pandas as pd

original_df = pd.read_csv('C:/Users/maran/Desktop/diff_BUM_HMM_Project/Github/diff_BUM_HMM/Xist_in vivo_vs_ex vivo_diff_BUM_HMM_analysed.txt', sep="\t", header=0)
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

print(type(A))

#print(A)
#print(A.loc[16850:16890,:])


#for iter in range(5):
    #print(iter)
    #subprocess.call("Rscript --vanilla /Users/maran/Desktop/diff_BUM_HMM_Project/Github/diff_BUM_HMM/diffBUM_HMM_0.6.R", shell=True)

temporary_df = pd.read_csv('C:/Users/maran/Desktop/diff_BUM_HMM_Project/Github/diff_BUM_HMM/Xist_in vivo_vs_ex vivo_diff_BUM_HMM_analysed_gaussian_noise.txt', sep="\t", header=0)
temporary_df= temporary_df.rename_axis('nt').reset_index()
temporary_df['filteredUM'] = [0 if x<0.90 else x for x in temporary_df['UM']]
temporary_df['filteredMU'] = [0 if x<0.90 else x for x in temporary_df['MU']]

new_df= temporary_df[['filteredUM','filteredMU']]

B = pd.DataFrame({n: new_df.T[col].nlargest(top_n).index.tolist() for n, col in enumerate(new_df.T)}).T


counter=0
other_counter=0

#for i in range(len(A.index)):
    #for j in range(0, len(A.columns)):
        #print(A.values[i,0], B.values[i,0])
    

#print(A.iloc[16850:16851,:1].values[0][0])

#print(A.loc[16850:16851,:1].values[0][0])

#print(A)

#for index, row in A.iterrows():
#val = A['0'].values[0]
#print(val)   


for i in range(len(A.index)):
    if A.iloc[i:i+1,0:1].values[0][0] == B.iloc[i:i+1,0:1].values[0][0]:
        counter=counter+1
    else:
        other_counter=other_counter+1

print(counter)
print(other_counter)







        #print(new_df_data.values[i,j],new_df.value[i,j])
#del df




    




