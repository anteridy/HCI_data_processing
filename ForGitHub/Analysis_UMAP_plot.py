# -*- coding: utf-8 -*-
"""
Created on Tue Mar 16 17:43:07 2021

@author: akamnev
"""
'''=====IMPORT LIBRARIES======================================='''
from matplotlib import  pylab as plt
from matplotlib import cm
import numpy as np
from numpy import ma
import itertools
import seaborn as sns
import os
import pandas as pd
# import scipy.spatial as sp, scipy.cluster.hierarchy as hc

import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib import cbook
from matplotlib import cm
# from sklearn.decomposition import PCA

#matplotlib module 'matplotlib.cbook' has no attribute 'iterable',
# apparently this is due to anaconda installing few matlab packages (>1, one by pip and one by anaconda)
# to check run: 
# $ conda list matplotlib
# if more then 1 instance found, deinstall:
# $ pip uninstall matplotlib
# then reinstall needed version (e.g. 3.2.2) as:
# $ conda install matplotlib=2.1.1

#%%Functions
def forcepath(path):
    """checks if directory exists at path <path> and, if not, makes one"""
    if not os.path.isdir(path):
        os.makedirs(path)
        
#%%
'''=====INPUTS & DEFINITIONS================================================'''
#PATHS------------------------------------------------------------------------
#Path to the raw data
DataPath = 'D:\\Temp\\APID\\APID_5_Screen2\\Morpho\\'
Folder = '8_SignatureAnalysis\\UMAPs\\Plt_'
ResFolder = Folder

#plates to process
# plates = [2,5,8,11]
plates = [2]

for plate in plates:

    plt_id = 'filtered_core'
    
    #data save to
    ResPath = DataPath+ResFolder+str(plate)+'\\Plots\\Core\\'
    forcepath(ResPath)
    
    
    #File names
    FileUMAP =  'UMAP_data_core.csv'
    
    index_feature_start = 13
    
    #plot settings:
    mrks = 4
    alph = 0.3
    #----------------
    
    doIndividual =True
    
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
    
    #%% Processing loop
    print('Plate='+str(plate))
    #%% UMAP plot
    print ("Loading data...")
    #open csv file with data for plate i as pd dataframe
    path = DataPath+Folder+str(plate)+'\\'+FileUMAP 
    data = pd.read_csv(path)
    data_cols = list(data)
      
    # make color table
    labels = data['Genotype']
    # print(set(labels))
    colors = []
    for l in labels:
        colors.append(colors_donor[l])
    
    #save which donors belong to which mutation
    #get order of groups
    donors =  list(set(data['Patient_ID']))
    mutations = list(set(data['Genotype']))
    mutation_to_donor = {}
    for m in mutations:
        mutation_to_donor[m] = []    
    for d in donors:
        mutation = data.loc[data['Patient_ID'] == d]['Genotype'].values[0]
        mutation_to_donor[mutation].append(d)    
        
    #%%    
    ''' Plot Patient Groups ----------------------------------------------'''
    #make fugure object
    fig = plt.figure()
    ax = plt.subplot(111)
    
    #plot scatter
    ax.scatter(data['UMAP1'], data['UMAP2'], color=colors, alpha=alph, s=mrks)
    
    #Make legend based on the labels used-------------------
    # Get all mutation types
    donor_groups = list(set(data['Genotype']))
    donor_groups.sort()
    patches = []
    for key in donor_groups:
        patches.append(mpatches.Patch(color=colors_donor[key], label=key))
    #---------------------------------------------------------
        
    #Shrink current axis's height by 20% on the right side
    box = ax.get_position()
    ax.set_position([box.x0, box.y0,
                 box.width*0.8, box.height])
    
    #Modify plot
    plt.title('UMAP - '+plt_id+' - Plt_'+str(plate))
    plt.xlabel('UMAP1')
    plt.ylabel('UMAP2')
    
    #add ledend    
    plt.legend(handles=patches, fontsize=10, frameon=False,
               title='Genotype', bbox_to_anchor=(1, 1), loc='upper left')
    
    #save figure       
    plt.savefig(ResPath+'All.png')   
    plt.close()
    
     #%%
    '''====Overlay ind donors (per group) over UMAP (core) ========= '''
    if doIndividual:  
        #status
        print ("Making UMAPs (all + overlay of mutation)...") 
          
        #Get labels
        labels_mutations = data['Genotype']
        labels_donor = data['Patient_ID']
        
        #set all mutations
        all_mutations = list(set(labels_mutations)) \
            
        X = data.values[:,-2:]
        
        # Make standard scatter plot output for each mutation
        for m in all_mutations:
            print(m)
            donors_within_mutation = mutation_to_donor[m]
            
            #check donor N to pick color table
            if len(donors_within_mutation)>10:
                colormap = cm.get_cmap('tab20')
            else:
                colormap = cm.get_cmap('tab10')
            
            donor_colors = {}
            for num,d in enumerate(donors_within_mutation):
                donor_colors[d] = colormap(num)
            
            #make fugure object
            fig = plt.figure()
            ax = plt.subplot(111)
            
            for x_val,mut,do in zip(X,labels_mutations,labels_donor):
                if mut == m:
                    ax.scatter(x_val[0], x_val[1], color=donor_colors[do], alpha=0.6, s=mrks)
                else:
                    ax.scatter(x_val[0], x_val[1], facecolors='none', color='grey', alpha=0.3, s=mrks)
                
        
            #Make legend based on the labels usedJu
            patches = []
            for key in donor_colors:
                patches.append(mpatches.Patch(color=donor_colors[key], label=key))
                
            #shrink plot to fit legend
            # Shrink current axis's height by 20% on the right side
            box = ax.get_position()
            ax.set_position([box.x0, box.y0,
                         box.width*0.8, box.height])
            
            #add title
            # plt.title('Plt_'+str(plate)+'_'+m+' '+metric+'_nb'+str(nb)+'mind'+str(mind)+'nc'+str(ncomp)+'seed'+str(rnds))
            plt.title(m+' UMAP - '+plt_id+' - Plt'+str(plate))
            
            #add ledend    
            plt.legend(handles=patches, fontsize=10, frameon=False,
                       title='Donors', bbox_to_anchor=(1, 1), loc='upper left')
            
            #save fig
            plt.savefig(ResPath+m+'.png')
            #plt.show()
            plt.close()
    '''========================================================================'''
#%% END









