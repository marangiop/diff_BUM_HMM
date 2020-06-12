#Extracting output
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

#35s MOLECULE
counter=0
list_output_values=[]
with open('output_35S_150_permutations') as topo_file:
    for index, line  in enumerate(topo_file):
        
        if index > 3304:
            split=line.split('.')
            if split[0] == '0':
                #print('ciaone')
                new_list =[]
                new_list.append('35S')
                new_list.append(float(line))
                counter=counter+1
                list_output_values.append(new_list)
                del new_list


df = pd.DataFrame(list_output_values, columns = ['Molecule', 'Prediction mismatch rate']) 
print(df)
print('Number of perturbations is' + str(counter))
del list_output_values
del counter

#----------------------------------------------------

counter=0
list_output_values=[]
with open('outputlatest') as topo_file:
    for index, line  in enumerate(topo_file):
        
        if index > 3276:
            split=line.split('.')
            if split[0] == '0':
                #print('ciaone')
                new_list =[]
                new_list.append('Xist')
                new_list.append(float(line))
                counter=counter+1
                list_output_values.append(new_list)
                del new_list

df_2 = pd.DataFrame(list_output_values, columns = ['Molecule', 'Prediction mismatch rate']) 
print('Number of perturbations is' + str(counter))

df_final=pd.concat([df,df_2], axis=0)

#ax = sns.boxplot(data=df_final)
ax = sns.boxplot(x="Molecule", y="Prediction mismatch rate", data=df_final)
plt.show()
ax.figure.savefig("boxplot_xist_35S_perturbations.png")
