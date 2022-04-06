# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 14:46:17 2020

@author: akamnev

==============================================
#------------INTRO--------------------------------
==================================================
Script to get images of cells within gate
"""
#%%Import 

#OS functions
import os #working with files and folders

#Dataframes and numberical operations
import numpy as np
import pickle #Create portable serialized representations of Python objects
import pandas as pd #pandas dataframe toolkit for large datasets
# import random  
from random import sample 


#plots
import matplotlib.pyplot as plt #main plotting engine
import seaborn as sns           # makes plots look nicer

#Image processing tools
from skimage import io, util

#%% Functions

def forcepath(path):
    """checks if directory exists at path <path> and, if not, makes one"""
    if not os.path.isdir(path):
        os.makedirs(path)
        
#%%
"""==============================================
--------INPUT - Raw dataset(s)-------------------------
================================================"""
DataPath = 'D:\\Temp\\APID\\APID_5_Screen2\\Morpho\\'
DataFolder = '3_Raw_data_mod\\' 
DataFile = 'SC_raw_plt'

plate = 2

#Columns to remove before saving-------------------------
FeatFolder = '4_Metadata\\'
# FeatName_cd8 = 'CD8_features.csv'
# FeatName_unwanted = 'UnwantedColumns.csv'
FeatName_unwanted_meta = 'UnwantedMetadata.csv'
FeatName_wanted = 'CoreFeatures.csv'
#------------------------------------------

#Other inputs:-----------------------------------------------
lastMetCol = 10 #last column of metadata in the table to split df
#------------------------------------------------------------

#path to images
ImgFolder = '1_Zproj\\Pt'
RawImgPath = DataPath+ImgFolder+str(plate).zfill(2)+'\\'

#Where to save plots
PltFolder = '6_ImgClassifier\\CellImages\\'
# PltPath = 'D:\\Temp\\APID\\APID_5_Screen2\\Morpho\\4_QC\\GatedCells\\RawImages\\' #for final plots
ResPath = DataPath+PltFolder
forcepath(ResPath)

#Cropping parameters
crop_width = 100 #pxl

#%% Load Morpho data
#status
print('Loading morpho data...')

#Load pkl file with morpho data
path = DataPath+DataFolder+DataFile+str(plate).zfill(2)+'.pkl'
with open(path, 'rb') as f:
        df_raw = pickle.load(f)
#show list of attributes in df
df_attr = list(df_raw)

#%%Load csv with special columns
foldpath = DataPath+FeatFolder

path = foldpath+FeatName_wanted
wanted_cols = pd.read_csv(path, index_col=None, header=0)['Feature'].astype(str)

path = foldpath+FeatName_unwanted_meta
unwanted_cols_meta = pd.read_csv(path, index_col=None, header=0)['Unwanted_metadata'].astype(str)
# unwanted_cols.replace('\n',' ', regex=True) 

#%% reduce dataset to wanted columns
#Move cell XY coordinatres to metadta------------------
df = df_raw
center_x = df['CellFoot_AreaShape_Center_X']
center_y = df['CellFoot_AreaShape_Center_Y']
df.drop(labels=['CellFoot_AreaShape_Center_X'], axis=1,inplace = True)
df.drop(labels=['CellFoot_AreaShape_Center_Y'], axis=1,inplace = True)
df.insert(lastMetCol+1, 'ObjectCenterX', center_x)
df.insert(lastMetCol+2, 'ObjectCenterY', center_y)
#------------------------------------------------------
    
#get metadata
df_meta = df.iloc[:,0:lastMetCol+3].drop(unwanted_cols_meta,axis=1)

#get features
df_feat = df[wanted_cols]

# merge back
df_final = pd.concat([df_meta,df_feat],axis=1)

#%% Get Sample
sampleN = 10 #number of samples
sample_size = 100 #number of objects per sample

df2sample = df_final

for i in range(1,sampleN):
    #get rnd sample
    df_sample = df2sample.sample(sample_size,replace=False, axis=0)
    #drop sample cases from the stock
    df2sample.drop(df_sample.index,inplace=True)
    #reser index in the sample
    df_sample.reset_index(drop=True)
    df_sample.index = np.arange(i*sample_size+1, (i+1)*sample_size+1)
    
    # stop
    #%% Find cells from sample
        
    #get XY of cell center------------------------------------------------------
    #---------------------------------------------------------------------------
    img_sample = []
    xc_array = []
    yc_array = []
    x_start_array = []
    y_start_array = []
    x_end_array = []
    y_end_array = []
    feature_val_array = []
    # gatedCells_df = pd.DataFrame()
    for obj in list(range(0,len(df_sample))):
        #get image ID
        r = str(df_sample['Metadata_Row'].iloc[obj])
        c = str(df_sample['Metadata_Column'].iloc[obj])
        f = str(df_sample['Metadata_FOV'].iloc[obj])
        
        img_id = 'r'+r.zfill(2)+'c'+c.zfill(2)+'f'+f.zfill(2)
        img_sample.append(img_id)
        
        #get cell XY center
        xc = df_sample['ObjectCenterX'].iloc[obj]
        yc = df_sample['ObjectCenterY'].iloc[obj]
        xc_array.append(xc)
        yc_array.append(yc)
        
        #get cropping window
        #X start
        if xc <= crop_width/2:
            x_start = 0
        else:
            x_start = xc - crop_width/2
        #X end
        if xc >= 1080 - crop_width/2:
            x_end = 1080
        else:
            x_end = xc + crop_width/2
        #Y start    
        if yc <= crop_width/2:
            y_start = 0
        else:
            y_start = yc - crop_width/2
        #Y end
        if yc >= 1080 - crop_width/2:
            y_end = 1080
        else:
            y_end = yc + crop_width/2
            
        x_start_array.append(x_start)
        y_start_array.append(y_start)
        x_end_array.append(x_end)
        y_end_array.append(y_end)
            
    # stop
    #%% Add data to resulting df
    df_sample.insert(loc=0, column='Image_ID', value=img_sample)
    df_sample.insert(loc=0, column='Ystart', value=y_start_array)
    df_sample.insert(loc=0, column='Xstart', value=x_start_array)
    df_sample.insert(loc=0, column='Xend', value=x_end_array)
    df_sample.insert(loc=0, column='Yend', value=y_end_array)
    
    #%%save result
    path = ResPath+'Data\\'
    forcepath(path)
    df_sample.to_csv(path+'Sample_'+str(i+1)+'.csv', index_label = 'Cell_Index')
    
    # stop
    #-------------------------------------------------------------------------------
    #%%
    # stop
    # Find and crop out cells from original image-------------------------------
    #---------------------------------------------------------------------------
    channels = [1,2,3,4]
    df = df_sample
    
    #Make dir for the donor
    #donor_path = ResPath+donor+'//'
    #os.mkdir(donor_path)
    
    for obj in list(range(0,len(df),1)):
        
        #Get img id
        img_id = df['Image_ID'].iloc[obj]
        
        #get cropping window
        x1 = int(df['Xstart'].iloc[obj])
        x2 = int(df['Xend'].iloc[obj])
        y1 = int(df['Ystart'].iloc[obj])
        y2 = int(df['Yend'].iloc[obj])
        
        for channel in channels:
            #make file name:
            img_name = img_id+'ch'+str(channel)+'pt'+str(plate).zfill(2)+'.tif'
            path = RawImgPath + img_name
            
            #load image
            img = io.imread(path)
            #show image
            #io.imshow(img)
            #io.show()
            
            # crop (for uncompressed or compressed with new version of Anaconda tiffs)
            crop = img[y1:y2,x1:x2] #rows = Y & cols = x in the dataframe
            
            # #crop (for compressed tiffs)
            # '''latest MIP script with compression produced multidimentional array as tif
            # image sits at level 4: img[0][0][0][1:1080,1:1080]'''
            # crop = img[0][0][0][y1:y2,x1:x2] #rows = Y & cols = x in the datafram
           
            # Save image
            img_name_fin = img_id +'obj'+str(obj)+'ch'+str(channel)+'_crop.tif'
            folder = ResPath+'Images//Sample_'+str(i+1)+'//'
            forcepath(folder)
            io.imsave(folder+img_name_fin, crop)
    # stop
    #-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
        
#%%status
print('----------------------------------------------------------------')
print('---> Finished. Results are saved.')






















