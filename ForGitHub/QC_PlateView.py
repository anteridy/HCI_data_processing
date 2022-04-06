# -*- coding: utf-8 -*-
"""
Created on Wed Aug 21 13:23:51 2019

@author: akamnev

Aim: to view QC params as a 384 w/pl layout heatmap
"""
#%%
#OS functions
import os #working with files and folders

#Dataframes and numberical operations
import numpy as np
from sklearn import preprocessing #for scaling dfs
import pickle #Create portable serialized representations of Python objects
import pandas as pd #pandas dataframe toolkit for large datasets

#plots
import matplotlib.pyplot as plt #main plotting engine
import seaborn as sns           # makes plots look nicer

#statistics
import scipy

#%%
"""==============================================
--------FUNCTIONS-------------------------
================================================"""
def forcepath(path):
    """checks if directory exists at path <path> and, if not, makes one"""
    if not os.path.isdir(path):
        os.makedirs(path)
"""=====END OF FUNCTIONS========================"""

#%%---PATHS---------------------------
#Path to the data
DataPath = 'D:\\Temp\\APID\\APID_5_Screen2\\Morpho\\4_QC\\'
DataName = 'SC_raw_pooled.pkl'
#Where to save plots
PltPath = 'D:\\Temp\\APID\\APID_5_Screen2\\Morpho\\4_QC\\PlateView\\' #for final plots
forcepath(PltPath)

#plates to process
#plates = ['01']
#plates = ['01','02','03','04','05','06','07','08','09','10','11','12']

#%%
'''Measurements to plot:'''
#imageQuality_Features =['CellFoot_AreaShape_Area',
#                        'Nuclei_AreaShape_Area',
#                        'CellFoot_Intensity_IntegratedIntensity_Actin',
#                        'Nuclei_Intensity_IntegratedIntensity_DNA']  
#imageQuality_Features =['CellFoot_AreaShape_Area']
#vmin = 500
#vmax = 3500
#imageQuality_Features =['Nuclei_AreaShape_Area']
#vmin = 300
#vmax = 2000
imageQuality_Features =['CellFoot_Intensity_MeanIntensity_Actin']
vmin = 0.001
vmax = 0.030
# imageQuality_Features =['Nuclei_Intensity_IntegratedIntensity_DNA']
# vmin = 1
# vmax = 30

#%% Load data
#status
print('Loading data...')
#open .pkl
path = DataPath+DataName
with open(path, 'rb') as f:
        df = pickle.load(f)

#show list of attributes in df
df_attr = list(df)

#%% Process plates
plates = df['Metadata_Plate'].unique()

for plate in plates:
    #status
    print('')
    print('Processing Plate '+str(plate)+'...')
    print('===================================')
    
    #get data for the plate
    df_plate = df.loc[df['Metadata_Plate']==plate]    
    data = df_plate
    
    #Full or part of the plate to plot:
    row_start = 1
    row_end = 16
    rows = range(row_start,row_end+1)
    
    col_start = 1
    col_end = 24
    cols = range(col_start,col_end+1)   
    
    for feature in imageQuality_Features:
        #print f    
        all_values = []
        for r in rows:
            tmp = []
            for c in cols:
                mean_well = data.loc[(data['Metadata_Row'] == r) & (data['Metadata_Column'] == c)][feature].astype('float').values.mean()
                tmp.append(mean_well)
            all_values.append(tmp)
        
        #Make figure window  
        fig=plt.figure() #?
        ax = plt.gca() #?
        #Plot
        im = ax.imshow(all_values, cmap='magma', interpolation='nearest'
                        ,vmin=vmin, vmax=vmax
                       )
        #title
        plt.title('Pt_'+str(plate)+' ave of '+feature)
        #Y axis
        plt.yticks(list(range(row_start-1,row_end)),list(map(str, list(range(row_end+1))[1:])))
        plt.ylabel('Row')
        #X axis
        plt.xticks(list(range(col_start-1,col_end)),list(map(str, list(range(col_end+1))[1:])))
        plt.xlabel('Column')
        #colorbar
        """# Create an axes for colorbar. The position of the axes is calculated based on the position of ax.
        # You can change 0.01 to adjust the distance between the main image and the colorbar.
        # You can change 0.02 to adjust the width of the colorbar.
        # This practice is universal for both subplots and GeoAxes."""
        cax = fig.add_axes([ax.get_position().x1+0.01,ax.get_position().y0,0.02,ax.get_position().height])
        plt.colorbar(im, cax=cax)
        #Layout
        plt.tight_layout(rect=[0,0,0.95,1])
        
#        stop
        
        #Save figure
        plt.savefig(PltPath+'Pt'+str(plate)+'_'+feature+'.png')
        plt.close()
#status
print('----------------------------------------------------------------')
#status
print('---> Finished. Plate heatmaps are done.')


        
    #%%
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    