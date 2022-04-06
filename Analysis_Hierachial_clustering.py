# -*- coding: utf-8 -*-
"""
Created on Fri Jan  8 12:57:45 2021

@author: akamnev
Based on original scrip from Michael Caldera "Analyse pre processed feature files" 
-> "2_AnalyseResultFiles" Jupiter script.
-> now split into individual script to plot specific plots only
"""

#%%
'''=====IMPORT LIBRARIES======================================='''
from matplotlib import  pylab as plt
from matplotlib import cm
import scipy.stats as stats
import numpy as np
from numpy import ma
import itertools
import seaborn as sns
import os
from scipy import spatial, linalg
from sklearn.preprocessing import StandardScaler
import seaborn as sns
import pandas as pd
import scipy.spatial as sp, scipy.cluster.hierarchy as hc

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib import cbook
from matplotlib import cm
from sklearn.decomposition import PCA

import umap

#time stamp
from datetime import datetime

#matplotlib module 'matplotlib.cbook' has no attribute 'iterable',
# apparently this is due to anaconda installing few matlab packages (>1, one by pip and one by anaconda)
# to check run: 
# $ conda list matplotlib
# if more then 1 instance found, deinstall:
# $ pip uninstall matplotlib
# then reinstall needed version (e.g. 3.2.2) as:
# $ conda install matplotlib=2.1.1

#%%
'''=====STANDARD FUNCTIONS & CLASSES======================================='''
class MidPointNorm(Normalize):    
    def __init__(self, midpoint=0, vmin=None, vmax=None, clip=False):
        Normalize.__init__(self,vmin, vmax, clip)
        self.midpoint = midpoint

    def __call__(self, value, clip=None):
        if clip is None:
            clip = self.clip

        result, is_scalar = self.process_value(value)

        self.autoscale_None(result)
        vmin, vmax, midpoint = self.vmin, self.vmax, self.midpoint

        if not (vmin < midpoint < vmax):
            raise ValueError("midpoint must be between maxvalue and minvalue.")       
        elif vmin == vmax:
            result.fill(0) # Or should it be all masked? Or 0.5?
        elif vmin > vmax:
            raise ValueError("maxvalue must be bigger than minvalue")
        else:
            vmin = float(vmin)
            vmax = float(vmax)
            if clip:
                mask = ma.getmask(result)
                result = ma.array(np.clip(result.filled(vmax), vmin, vmax),
                                  mask=mask)

            # ma division is very slow; we can take a shortcut
            resdat = result.data

            #First scale to -1 to 1 range, than to from 0 to 1.
            resdat -= midpoint            
            resdat[resdat>0] /= abs(vmax - midpoint)            
            resdat[resdat<0] /= abs(vmin - midpoint)

            resdat /= 2.
            resdat += 0.5
            result = ma.array(resdat, mask=result.mask, copy=False)                

        if is_scalar:
            result = result[0]            
        return result

    def inverse(self, value):
        if not self.scaled():
            raise ValueError("Not invertible until scaled")
        vmin, vmax, midpoint = self.vmin, self.vmax, self.midpoint

        if cbook.iterable(value):
            val = ma.asarray(value)
            val = 2 * (val-0.5)  
            val[val>0]  *= abs(vmax - midpoint)
            val[val<0] *= abs(vmin - midpoint)
            val += midpoint
            return val
        else:
            val = 2 * (value - 0.5)
            if val < 0: 
                return  val*abs(vmin-midpoint) + midpoint
            else:
                return  val*abs(vmax-midpoint) + midpoint
        

        
def reject_outliers_2(data, m = 3.):
    d = np.abs(data - np.median(data))
    mdev = np.median(d)
    s = d/(mdev if mdev else 1.)
    return [data[i] for i in range(0,len(data)) if s[i] < m]

def forcepath(path):
    """checks if directory exists at path <path> and, if not, makes one"""
    if not os.path.isdir(path):
        os.makedirs(path)

#%%
'''=====INPUTS & DEFINITIONS================================================'''
#PATHS------------------------------------------------------------------------
#Path to the raw data
DataPath = 'D:\\Temp\\APID\\APID_5_Screen2\\Morpho\\'
Folder = '7_Filtered_data\\Sample_wrap_liveCD8\\Plt_'
ResFolder = '8_SignatureAnalysis\\HierClust\\Plt_'
File =  'FinalresultFile_CD1.csv'

#path to core features:
# path_core= '4_Metadata\\PinnedFeatures.csv'

#Donors to remove
# unw_donors = ['PID_2599 (CD4/8 mix)', 'ND_TL115 (w7)']
unw_donors = [] #already removed at feature filtration

#plates to process
plate = 8

#data save to
ResPath = DataPath+ResFolder+str(plate)+'\\'
forcepath(ResPath)

#start of features
index_feature_start = 13

#do plots of indivudual features for each donor?
doIndPlots = False

## Define colors for the individual conditions-----------------------------
colors_donor = {'ND':'#525252',
                
                'ARPC1B':'#ff7f00',
                'DOCK11':'#ffff33',
                'HEM1':'#4daf4a',
                'MSN':'#377eb8',
                'PAK2':'#f781bf',
                'PIK3CG':'#984ea3',
                'PSTPIP1':'#00FFF7',
                'THEMIS2':'#33ff33',
                'WAS':'#a65628',
                'WDR1':'#e41a1c',         
               }

# feature_colors = {'AreaShape':'#008000', #green
#                 'Intensity':'#FFFF00', #yellow
#                 'RadialDistribution':'#00FFFF', #cyan
#                 'Texture':'#FF00FF', #magenta
#                 'Neighbors':'grey',
#                 'Number':'black',
#                 'number':'black',
#                 'Parent':'black',
#                   'Location':'black'}

# feature_colors = {
#                 #channels
#                 'Actin':'#008000', #green
#                 'DNA':'#00FFFF', #cyan
#                 'Marker':'#FF00FF', #magenta
#                 #geometry
#                 'AreaShape':'grey' #green
#                 }

feature_colors = {
                #channels
                'CellFoot':'#FF00FF', #magenta
                'Nuclei':'#00FFFF' #cyan
                }
#----------------------------------------------------------------------------
#%% Processing loop
print('')
print('')
print('Plate = '+str(plate))
# datetime object containing current date and time
now = datetime.now()
print("start =", now)
print('=========================================================================')
print('')
#%%
'''=====LOAD DATA================================================'''
#status
print ("Loading data...")
# #load pinned features
# path = DataPath + path_core
# core_features = pd.read_csv(path)
# core_features_list = core_features['Pinned_features'].tolist()

#get values for the specific group
#open csv file and save as dataframe
path = DataPath+Folder+str(plate)+'\\'+File
data_raw = pd.read_csv(path)
data_raw = data_raw.dropna()

#drop unwanted donors:
df = data_raw
don_ind = df['Patient_ID'].isin(unw_donors)
data = df[~don_ind]   

# Get all donor types
donors = list(set(data['Patient_ID']))
donors.sort()

# Get all mutation types
donor_groups = list(set(data['Genotype']))
donor_groups.sort()

#get all features
features = data.columns[index_feature_start:]
# features = core_features_list 

#%% Get core features only
# metadata = data.iloc[:,0:index_feature_start-1]
# core_data = data.loc[features]

# stop
#%%
'''====Create Dictonary of valid values ===================================
-> Remove clear wrong values'''

#status
print ("Making dictionaries for values...")

#get order of groups
donors =  list(set(data['Patient_ID']))
mutations = list(set(data['Genotype']))

valid_values_donors = {}
valid_values_mutations = {}
donor_to_mutation = {}
mutation_to_donor = {}

for m in mutations:
    valid_values_mutations[m] = {}
    mutation_to_donor[m] = []
    for f in features:
        valid_values_mutations[m][f] = []

for d in donors:
    valid_values_donors[d] = {}
    
    #get according mean
    mutation = data.loc[data['Patient_ID'] == d]['Genotype'].values[0]
    
    #save which donors belong to which mutation
    mutation_to_donor[mutation].append(d)
    donor_to_mutation[d] = mutation
    
    for f in features:      
        #Get all values
        values = data.loc[data['Patient_ID'] == d][f]
        
        #keep only good values --> this outlier detection is slightly different
        # than before (outlier for ALL NDs are calculated, and here per donor)
        valid_values = reject_outliers_2(list(values))
        
        #Calculate mean
        mean_val = np.median(valid_values)
        
        #add mean
        valid_values_donors[d][f] = mean_val
        
        #add to mutation mean
        valid_values_mutations[mutation][f].append(mean_val)
        
#Calculate median per mutation  
for m in mutations:
    for f in features:
        valid_values_mutations[m][f] = np.median(valid_values_mutations[m][f])
        
#Create column colors (= colors for the features)
col_colors = []
# stop

for f in features:
    for key in feature_colors.keys():
        if f.find(key) != -1:
            color = feature_colors[key]
            col_colors.append(color)

#Status  
print("Performed Outlier detection, created result dictionaries")
# stop
#%%
'''====CLUSTERMAPS========================================================='''

#status
print ("Making clustermaps of indivudual donors")

x_donors = []
labels_donors = []
row_colors_donors = []

for d in donors:
    
    #extract metadata
    labels_donors.append(d)
    mutation = data.loc[data['Patient_ID'] == d]['Genotype'].values[0]
    row_colors_donors.append(colors_donor[mutation])
    
    #Get values vor all features 
    values = []
    for f in features:
        values.append(valid_values_donors[d][f])
          
    x_donors.append(values)

max_val = min([max(np.array(x_donors).flatten()), 3])
min_val = max([min(np.array(x_donors).flatten()), -3])   

norm = MidPointNorm(midpoint=0.0)
# methods = ['complete','average','weighted','centroid','median','ward']
method = 'complete'
# metrics = ['braycurtis','canberra','chebyshev','cityblock','euclidean','minkowski']
metric = 'euclidean'
lut = 'seismic'
# for method in methods:
sns.clustermap(data=x_donors,vmin=min_val, vmax=max_val,method=method,metric=metric,
               cmap=lut,norm=norm, row_colors = row_colors_donors,
               col_colors=col_colors, yticklabels = labels_donors)
#sns.clustermap(data=x_donors,method='centroid',cmap='RdBu_r', row_colors = row_colors, col_colors=col_colors, yticklabels = labels_donors)

#save result
plt.savefig(ResPath+'HierClust_donors_plt'+str(plate)+'_'+method+'.png')
plt.close()


#%% Clustemap of individual donors
print ("Making clustermaps of donor groups")
x_mutations = []
labels_mutations = []
row_colors_mutations = []

for m in mutations:
    
    #extract metadata
    labels_mutations.append(m)
    row_colors_mutations.append(colors_donor[m])

    #Get values vor all features 
    values = []
    for f in features:
        values.append(valid_values_mutations[m][f])
    x_mutations.append(values)
    
    
max_val = min([max(np.array(x_mutations).flatten()), 3])
min_val = max([min(np.array(x_mutations).flatten()), -3])    

norm = MidPointNorm(midpoint=0)

# methods = ['complete','average','weighted','centroid','median','ward']
method = 'complete'
# metrics = ['braycurtis','canberra','chebyshev','cityblock','euclidean','minkowski']
metric = 'euclidean'
lut = 'seismic'
# lut = 'RdBu_r'

# for method in methods:
sns.clustermap(data=x_mutations,vmin=min_val, vmax=max_val,
               method=method, metric = metric,
               robust=True,#ignore outliers
               cmap=lut,norm=norm,
               row_colors = row_colors_mutations, col_colors=col_colors,
               yticklabels = labels_mutations)
#sns.clustermap(data=x_mutations,method='centroid',cmap='RdBu_r',  row_colors = row_colors, col_colors=col_colors, yticklabels = labels_mutations)
plt.savefig(ResPath+'HierClust_mutations_plt'+str(plate)+'_'+method+'.png')
plt.close()

# stop
#%%
'''====INDIVIDUAL FEATURES=================================================='''
########################################
#Make Feature Analysis (Barplot)
########################################
if doIndPlots:
    #status
    print ("Plotting individual features...")
    
    
    ensure_dir(ResPath+'IndividualFeaturesCore\\')
    mutations.sort()
    mutations.remove('ND')
    mutations.append('ND')
    
    #make csv file to write in results
    fp_out = open(ResPath+'IndividualFeatures\\Overview.csv','w')
    fp_out.write('Feature,'+','.join(donors)+'\n')
    
    #write in data
    for num,f in enumerate(features):
    
        #group values by Donor group (take mean and std)
        #y_values = data.groupby(['Patient_ID']).mean()[f]
        #errors = data.groupby(['Patient_ID']).std()[f]
    
        #get order of groups
        #groups =  list(data.groupby(['Patient_ID']).mean().index)
        y_values = []
        errors = []
        labels = []
        colors = []
        for m in mutations:
            for d in mutation_to_donor[m]:
                labels.append(d+'_'+m)
                colors.append(colors_donor[m])
                y_values.append(valid_values_donors[d][f])
        
        #make bar plot
        #plt.bar(range(0,len(y_values)),y_values, yerr = errors, color = 'grey')
        barlist = plt.bar(range(0,len(y_values)),y_values)
        
        for num,bar in enumerate(barlist):
            
            bar.set_color(colors[num])
        
        plt.axhline(0, ls='--', color='grey')
        plt.axhline(0, ls='--', color='grey')
        plt.xticks(range(0,len(y_values)), labels, rotation=90, size=2)
        plt.ylabel(f + '[+/- std]')
        plt.xlabel('Mutations')
        plt.savefig(ResPath+'IndividualFeatures\\'+f+'_plt'+str(plate)+'.png')
        plt.close()
        
        fp_out.write(f+','+','.join([str(x) for x in y_values])+'\n')
        
    fp_out.close()
#%%END












