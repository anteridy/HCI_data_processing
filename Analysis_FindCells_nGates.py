# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 14:46:17 2020

@author: akamnev

==============================================
#------------INTRO--------------------------------
==================================================
Script to get images of cells within gate
"""
#%%
"""==============================================
--------IMPORT PACKAGES & MAKE FUNCTIONS-------------------------
================================================"""
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

#%%
'''=======FUNCTIONS==============================================
=================================================================
'''
#clean NaNs from pandas dataframe
def cleanNaNs(df):
    """<clearNaNs> first removes columns which contain only NaNs,
    then removes rows with at least one NaN
    """
    
    #Remove Columns with all NaNs:
    df_1 = df.dropna(axis='columns', how='all')
    
    #Remove rows with at least one NaN    
    df_2 = df_1.dropna()
    
    return df_2

#get a subset of the panda DF
def df_find_subset(data = None, groups='All', attr=None):
    """Find elements in the dataframe which satisfy group requirements
    outputs df filtered by attr and groups
    name of the chosen group as a string
    data   : pandas df with raw data
    attr   : list of attribute names to take out ['attr1', 'attr2']
    groups : list of lists with group name and values ['Name', [val1, val2]]
    """ 
    
    #get columns out
    df = data
    if not attr:
        df_tmp = df
    else:
        df_tmp = pd.DataFrame()
        for attribute in attr:
            df_tmp = pd.concat([df_tmp, df[attribute]], axis=1)
        
    #get rows
    if groups != 'All':
        for index in range(0,len(groups)):
            df_tmp = df_tmp.loc[df_tmp[groups[index][0]].isin(groups[index][1])]
        
    return df_tmp

def forcepath(path):
    """checks if directory exists at path <path> and, if not, makes one"""
    if not os.path.isdir(path):
        os.makedirs(path)
        
#%%
"""==============================================
--------INPUT - Raw dataset(s)-------------------------
================================================"""
#%%---PATHS---------------------------
#Path to the morpho data
DataPath = 'D:\\Temp\\APID\\APID_5_Screen2\\Morpho\\'
DataFolder = '6_PreProcessed_data\\ImgClassifier\\SC_liveCD8_ann\\'
DataFile = 'SC_liveCD8_ann_plt'

#Data to remove
# unw_pids = ['DOCK11','PAK2','THEMIS2']
unw_pids = ['DOCK11','PAK2','THEMIS2', 'PSTPIP1', 'WDR1','MSN' ]
#Donors to remove
unw_donors = ['PID_2599 (CD4/8 mix)', 'ND_TL115 (w7)']

#Cropping parameters
crop_width = 200 #pxl

#Attributes to take out before gating
wanted_attributes = ['Metadata_Row',
                     'Metadata_Column',
                     'Metadata_FOV',
                     'ObjectCenterX',
                     'ObjectCenterY']

#%% Data to process

#Gates-ICAM1-GM (Plt2)-----------------------
plate = 2
# #modify input data
# param1 = 'CellFoot_Intensity_MeanIntensity_Actin'
# param2 = 'CellFoot_AreaShape_EOP'
# param3 = 'CellFoot_AreaShape_Area'
# data_1[param1] = data_1[param1].apply(lambda x: x*0.106) # pxl^2 to um^2
# data_1[param2] = data_1[param2].apply(lambda x: x*0.325) # pxl to um

#gates
'unquote required block before running the script'
genotype = 'ND'
doAll = True
sel_donors =[]
gates = [[param1,[0.021,0.024]],
            [param2,[0.24, 0.28]],
            [param3, [2000, 2400]]
            ]

# #High-F-act
# genotype = 'WAS'
# doAll = False
# sel_donors = ['PID_WAS1', 'PID_WAS2']
# gates = [[param1,[0.019,0.023]],
#             [param2,[0.22, 0.25]]
#             # [param3, [1750, 2100]]
#             ]
#Low-F-act
# genotype = 'WAS'
# doAll = False
# sel_donors =['PID_WAST2']
# gates = [[param1,[0.014,0.016]],
#             [param2,[0.22, 0.25]]
#             # [param3, [1750, 2100]]
#             ]

# genotype = 'HEM1' 
# gates = [[param1,[0.014,0.016]],
#             [param2,[0.29, 0.33]],
#             [param3, [1750, 2100]]
#             ]

# genotype = 'ARPC1B' #pt1
# gates = [[param1,[0.011,0.014]],
#             [param2,[0.25, 0.30]],
#             [param3, [1750, 2100]]
#             ]
#-----------------------------------

#Gates-ICAM1-IL2 (Plt5)-----------------------
plate = 5
#modify input data
param1 = 'CellFoot_Intensity_MeanIntensity_Actin'
param2 = 'CellFoot_AreaShape_EOP'
param3 = 'CellFoot_AreaShape_Area'
# data_1[param1] = data_1[param1].apply(lambda x: x*0.106) # pxl^2 to um^2
# data_1[param2] = data_1[param2].apply(lambda x: x*0.325) # pxl to um

#gates
# genotype = 'ND'
# doAll = True
# # sel_donors =[]
# gates = [[param1,[0.021,0.025]],
#             [param2,[0.29, 0.34]]
#             # ,[param3, [2700, 2800]]
#             ]

# genotype = 'WAS'
# doAll = True
# # sel_donors =[]
# gates = [[param1,[0.020,0.022]],
#             [param2,[0.28, 0.33]]
#             # ,[param3, [2700, 2800]]
#             ]

# genotype = 'HEM1' 
# doAll = True
# # sel_donors =[]
# gates = [[param1,[0.018,0.021]],
#             [param2,[0.35, 0.40]]
#             # ,[param3,[2800, 2900]]
#             ]

# genotype = 'ARPC1B' #pt1
# doAll = True
# # sel_donors =[]
# gates = [[param1,[0.015,0.018]],
#             [param2,[0.31, 0.37]]
#             # ,[param3, [2500, 2600]]
#             ]

# genotype = 'PIK3CG'
# doAll = True
# # sel_donors =[]
# gates = [[param1,[0.017,0.021]],
#             [param2,[0.28, 0.31]]
#             # ,[param3, [2800, 2900]]
#             ]
#-----------------------------------

#%% IS-GM (plate 8)
# #modify input data
# plate = 8
# param1 = 'CellFoot_Intensity_MeanIntensity_Actin'
# # param2 = 'CellFoot_AreaShape_Eccentricity' #don't gate on this to increase cell number
# param3 = 'CellFoot_AreaShape_EOP'
# param4 = 'CellFoot_AreaShape_Area'

# #Gates-IS (Plt8)-----------------------
# genotype = 'ND'
# gates = [[param1,[0.018,0.021]],
# #             [param2,[0.6, 0.65]],
#             [param3, [0.21, 0.25]],
#             [param4,[1400,2400]]
#               ]
          
# # genotype = 'WAS'
# # gates = [[param1,[0.014,0.017]],
# # #             [param2,[0.7, 0.75]], #high range pt 1,2
# #             [param3, [0.25, 0.27]],
# #             [param4,[1400,2400]]
# #             ]

# # genotype = 'HEM1' 
# # gates = [[param1,[0.014,0.016]],
# # #             [param2,[0.625, 0.675]],
# #             [param3, [0.3, 0.32]],
# #             [param4,[1400,2400]]
# #             ]

# # genotype = 'ARPC1B' #low pts
# # gates = [[param1,[0.0115,0.0135]],
# # #             [param2,[0.545, 0.575]],
# #             [param3, [0.21, 0.25]],#middle pt 2
# #             [param4,[1400,2400]]
# #             ]

# genotype = 'PIK3CG'
# doAll = True
# # sel_donors =[]
# gates = [[param1,[0.0150,0.0175]],
# #             [param2,[0.55, 0.6]],
#             [param3, [0.19, 0.22]]
#             # ,[param4,[1400,2400]]
#             ]
# #-----------------------------------

#%%Gates-IS+IL2 (Plt11)-----------------------
# #modify input data
# plate = 11
# param1 = 'CellFoot_Intensity_MeanIntensity_Actin'
# # param2 = 'CellFoot_AreaShape_Eccentricity' #don't gate on this to increase cell number
# param3 = 'CellFoot_AreaShape_EOP'
# param4 = 'CellFoot_AreaShape_Area'

# genotype = 'ND'
# doAll = True
# # sel_donors =[]
# gates = [[param1,[0.018,0.021]],
#             [param3, [0.22, 0.275]],
#             [param4,[1200,2700]]
#               ]
          
# genotype = 'WAS'
# gates = [[param1,[0.016,0.018]],
#             [param3, [0.225, 0.275]], #high pts 1,2
#             [param4,[1200,2700]]
#             ]

# genotype = 'HEM1' 
# gates = [[param1,[0.0155,0.0175]],
#             [param3, [0.25, 0.3]],
#             [param4,[1200,2700]]
#             ]

# genotype = 'ARPC1B' #high EOCP
# doAll = False
# sel_donors =['PID_ARPC1B_Pt7']
# gates = [[param1,[0.011,0.0125]],
#             [param3, [0.19, 0.24]],#low, pt 1,2
#             [param4,[1200,2700]]
#             ]

# genotype = 'ARPC1B' #low EOCP
# doAll = False
# sel_donors =['PID_ARPC1B_Pt1','PID_2599 (CD8)']
# gates = [[param1,[0.011,0.0125]],
#             [param3, [0.19, 0.24]],#low, pt 1,2
#             [param4,[1200,2700]]
#             ]

# genotype = 'PIK3CG'
# doAll = True
# # sel_donors =[]
# gates = [[param1,[0.018,0.021]],
#             [param3, [0.16, 0.19]],
#             [param4,[1200,2700]]
#             ]
#-----------------------------------

#%% Load Morpho data
#status
print('Loading morpho data...')

#Load pkl file with morpho data
path = DataPath+DataFolder+DataFile+str(plate)+'.pkl'
with open(path, 'rb') as f:
        data_raw = pickle.load(f)
#show list of attributes in df
data_attr = list(data_raw)

#Clean data---------------------------------
#Drop unwanted unwanted PIDs
df = data_raw
pid_ind = df['Genotype'].isin(unw_pids)
data_raw_cln = df[~pid_ind]
df = data_raw_cln
don_ind = df['Patient_ID'].isin(unw_donors)
data_1 = df[~don_ind]  
#-------------------------------------------

#path to images
ImgFolder = '1_Zproj\\Pt'
RawImgPath = DataPath+ImgFolder+str(plate).zfill(2)+'\\'

#%%
wanted_group = [['Genotype',[genotype]]]
#%%
#Where to save plots
PltFolder = '8_SignatureAnalysis\\CellImages\\Plt_'+str(plate)+'\\'+genotype+'\\'
# PltPath = 'D:\\Temp\\APID\\APID_5_Screen2\\Morpho\\4_QC\\GatedCells\\RawImages\\' #for final plots
ResPath = DataPath+PltFolder
forcepath(ResPath)

#append chosen features
for i in range(0,len(gates)):
    wanted_attributes.append(gates[i][0])

     
#%%Get gated cells:
#Get genotype
df_gen = data_1.loc[data_1['Genotype']==genotype]
#Get patients
if doAll:
    df_donor = df_gen
else:
    df_donor = df_gen.loc[df_gen['Patient_ID'].isin(sel_donors)]

# stop

df2gate = df_donor[wanted_attributes]
df = df2gate
#gate dataset
for i in range(0,len(gates)):
    attribute = gates[i][0]
    g_min = gates[i][1][0]
    g_max = gates[i][1][1]
    df_tmp1 = df.loc[(df[attribute]>=g_min) & (df[attribute]<=g_max)]
    df = df_tmp1

df_gated = df

#Check if cells found
cellN = len(df_gated)
if cellN > 0:
    print("Found "+str(cellN)+' cells fitting gates')
else:
    print('')
    print('[!] Found no cells fitting into gate [!]')
    print('')
   
# Get XY coordinates of cells to crop out
#Get random 25 cells if possible
obj_n = len(df_gated)
#obj_sample = list(range(1,obj_n,1))
if obj_n >= 25:  
    obj_array = list(range(0,obj_n,1))
    obj_sample = sample(obj_array,25)
else:
    obj_sample = list(range(0,obj_n,1))


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
gatedCells_df = pd.DataFrame()
for obj in obj_sample:
    #get image ID
    r = str(df_gated['Metadata_Row'].iloc[obj])
    c = str(df_gated['Metadata_Column'].iloc[obj])
    f = str(df_gated['Metadata_FOV'].iloc[obj])
    img_id = 'r'+r.zfill(2)+'c'+c.zfill(2)+'f'+f.zfill(2)
    img_sample.append(img_id)
    
    #get cell XY center
    xc = df_gated['ObjectCenterX'].iloc[obj]
    yc = df_gated['ObjectCenterY'].iloc[obj]
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
    
    #Get gating value
#    feature_val_array.append(df_gated[gate_attribute].iloc[obj])
    
    
#Add data to resulting df
gatedCells_df['Image_ID'] = img_sample
#gatedCells_df[gate_attribute] = feature_val_array
gatedCells_df['Xc'] = xc_array
gatedCells_df['Yc'] = yc_array
gatedCells_df['Xstart'] = x_start_array
gatedCells_df['Ystart'] = y_start_array
gatedCells_df['Xend'] = x_end_array
gatedCells_df['Yend'] = y_end_array
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------



# Find and crop out cells from original image-------------------------------
#---------------------------------------------------------------------------
channels = [1,2,3,4]
df = gatedCells_df

#Make dir for the donor
#donor_path = ResPath+donor+'//'
#os.mkdir(donor_path)

#%%crop raw images
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
        
        #crop (for uncompressed tiffs)
        crop = img[y1:y2,x1:x2] #rows = Y & cols = x in the dataframe
        
        #crop (for compressed tiffs)
        '''latest MIP script with compression produced multidimentional array as tif
        image sits at level 4: img[0][0][0][1:1080,1:1080]
        '''
        # crop = img[0][0][0][y1:y2,x1:x2] #rows = Y & cols = x in the datafram
        
        # Save image
        path = ResPath+'Raw//'
        forcepath(path)
        img_name = img_id +'obj'+str(obj)+'ch'+str(channel)+'_crop.tif'
        io.imsave(path+img_name, crop)
        
#%%crop segmentation results
# 'Choose this if in addition to raw data segmentation results are also needed'
# for obj in list(range(0,len(df),1)):
    
#     #Get img id
#     img_id = df['Image_ID'].iloc[obj]
    
#     #get cropping window
#     x1 = int(df['Xstart'].iloc[obj])
#     x2 = int(df['Xend'].iloc[obj])
#     y1 = int(df['Ystart'].iloc[obj])
#     y2 = int(df['Yend'].iloc[obj])
    
#     for channel in [2]:
#         #make file name:
#         img_name = img_id+'ch'+str(channel)+'pt'+str(plate).zfill(2)+'_segm.png'
#         path = SegmImgPath + img_name
        
#         #load image
#         img = io.imread(path)
#         #show image
#         #io.imshow(img)
#         #io.show()
        
#         #crop (for uncompressed tiffs)
#         crop = img[y1:y2,x1:x2] #rows = Y & cols = x in the dataframe
        
#         #crop (for compressed tiffs)
#         '''latest MIP script with compression produced multidimentional array as tif
#         image sits at level 4: img[0][0][0][1:1080,1:1080]
#         '''
#         # crop = img[0][0][0][y1:y2,x1:x2] #rows = Y & cols = x in the datafram
        
#         # Save image
#         path = ResPath+'Segm//'
#         forcepath(path)
#         img_name = img_id +'obj'+str(obj)+'ch'+str(channel)+'_crop.png'
#         io.imsave(path+img_name, crop)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
        
#%%status
print('----------------------------------------------------------------')
print('---> Finished. Results are saved.')
print("Found "+str(cellN)+' cells fitting gates')























