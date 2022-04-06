# -*- coding: utf-8 -*-
"""
Created on Fri Mar  8 12:58:31 2019

@author: akamnev
"""

"""==============================================
#------------INTRO--------------------------------
==================================================
Script to calculate relative position of subcellular structures:
    nuclei center to cell center
    Cellmask center to actin center

#-------------------------------------------------
==================================================
"""

#%%
"""==============================================
--------IMPORT PACKAGES-------------------------
================================================"""
import os #working with files and folders
import pandas as pd #pca toolkit
import numpy as np
import pickle #Create portable serialized representations of Python objects
#import math #for square root calc etc, use numpy instead for lists and panda datasets

"""-------------------------------------------------
==================================================
"""
#%%
"""==============================================
--------FUNCTIONS-------------------------
================================================"""
def forcepath(path):
    """checks if directory exists at path <path> and, if not, makes one"""
    if not os.path.isdir(path):
        os.makedirs(path)
"""=====END OF FUNCTIONS========================"""


#%%
"""--------INPUTS--------------------"""

# #Common paths:
# PKL_name = 'Results_raw_pooled_plt'
# Res_name = 'Results_raw_pooled_mod_plt'
# RawMorphoPklPath = 'D:\\Temp\APID\\APID_5_Screen1\\3_Python_res\\Raw_pooled\\'
# ModMorphoPklPath = 'D:\\Temp\\APID\\APID_5_Screen1\\3_Python_res\\Raw_pooled_mod\\'

# #plates to process
# #plates = ['01','02','03','04','05','06','07','08','09']
# plates = ['10','11','12']

#Path to the cluster results (plate folder):
DataPath = 'D:\\Temp\\APID\\APID_5_Screen2\\Morpho\\3_Raw_data\\'
DataName = 'SC_raw_plt'
#Where to write results
OutPath = 'D:\\Temp\\APID\\APID_5_Screen2\\Morpho\\3_Raw_data_mod\\'
forcepath(OutPath)
#Name for single-cell raw dataset:
ResName = DataName

# plates to process
# plates = ['01','02']
plates = ['01','02','03','04','05','06','07','08','09','10','11','12']

#Other inputs
lastcol_met = 10 #to split morpho params from nuclei dataset
idx = lastcol_met+1 #where to insert new measurements

#Quick test?
test = False #Open only first x rows
rown = 1000

#"""-----------------------------------------------------------------------------"""

#%% Process plates
for plate in plates:
    #status
    print('')
    print('Processing Plate '+plate+'...')
    print('===================================')
    
    #build path to plate raw data
    FilePath = DataPath+DataName + plate +'.pkl'   
    
    # Load data
    #status
    print('Loading data...')
    
    path = FilePath
    with open(path, 'rb') as f:
        df_raw = pickle.load(f)
    
    #attributes in the dataframe
    features_all = list(df_raw)
    
    #Reduce data size for testing
    if test:
        df_morpho_raw = df_raw.iloc[:rown] #[!]for testing only
    else:
        df_morpho_raw = df_raw
    # stop
    
    # Calculate euclidian distance btw cell and nuclei center-------------------
    #status
    print('Computing euclidean distances...')
    
    #cell center
    xc = df_morpho_raw['CellFoot_AreaShape_Center_X']
    yc = df_morpho_raw['CellFoot_AreaShape_Center_Y']
    
    #Results df
    df_tmp = pd.DataFrame()
    
    #calculate displacement from cell center (footprint)
    for channel in ['DNA', 'Actin', 'Marker']:
        #get object mass center
        x = df_morpho_raw['CellFoot_Location_CenterMassIntensity_X_'+channel]
        y = df_morpho_raw['CellFoot_Location_CenterMassIntensity_Y_'+channel]
        
        #calculate distance & add to the result df
        new_col = np.sqrt(np.power(xc-x,2)+np.power(yc-y,2))
        #insert columns back to the dataframe
        col_name = 'CellFoot_Intensity_Displacement_'+channel
        df_morpho_raw.insert(loc=idx, column=col_name, value=new_col)

    # stop
    #---------------------------------------------------------------------------    
    
    #Calculate perimeter-area parameters
    #status
    print('Computing perimeter-area parameters:')
    
    """
    A) Nuclear lobularity: =P/A
    B) Nuclear contour ratio: =(4*3.1416)*(A/(P)^2)
    C) Excess of nuclear perimeter (EONP):
    Step 1: X =P/(2*3.1416*(A/3.1416)^0.5)
    Step 2: Y =1/X
         
    Step 3: Z =1-Y
    ----------------------
    P - perimeter
    A - area
    Z - EOP (excess of perimeter)
    """
    #Calculate lobularity-----------------------------
    print('Computing lobularity...')
    #Nuclei
    x = df_morpho_raw['Nuclei_AreaShape_Perimeter'].astype('float') #perimeter
    y = df_morpho_raw['Nuclei_AreaShape_Area'].astype('float') #area
    new_col = x/y
    #insert columns back to the dataframe
    col_name = 'Nuclei_AreaShape_Lobularity'
    df_morpho_raw.insert(loc=idx, column=col_name, value=new_col)
    
    #Cell
    x = df_morpho_raw['CellFoot_AreaShape_Perimeter'].astype('float') #perimeter
    y = df_morpho_raw['CellFoot_AreaShape_Area'].astype('float') #area
    new_col = x/y
    #insert columns back to the dataframe
    col_name = 'CellFoot_AreaShape_Lobularity'
    df_morpho_raw.insert(loc=idx, column=col_name, value=new_col)
    #------------------------------------------------

    # stop
    #Calculate contour ratio------------------------
    print('Computing contour ratio...')
    x = df_morpho_raw['Nuclei_AreaShape_Perimeter'].astype('float') #perimeter
    y = df_morpho_raw['Nuclei_AreaShape_Area'].astype('float') #area
    new_col = (4*np.pi)*(y/(np.power(x,2)))
    #insert columns back to the dataframe
    col_name = 'Nuclei_AreaShape_ContourRatio'
    df_morpho_raw.insert(loc=idx, column=col_name, value=new_col)
    
    x = df_morpho_raw['CellFoot_AreaShape_Perimeter'].astype('float') #perimeter
    y = df_morpho_raw['CellFoot_AreaShape_Area'].astype('float') #area
    new_col = (4*np.pi)*(y/(np.power(x,2)))
    #insert columns back to the dataframe
    col_name = 'CellFoot_AreaShape_ContourRatio'
    df_morpho_raw.insert(loc=idx, column=col_name, value=new_col)
    #------------------------------------------------
    
    #Calculate lobularity-----------------------------
    print('Computing EOP...')
    x = df_morpho_raw['Nuclei_AreaShape_Perimeter'].astype('float') #perimeter
    y = df_morpho_raw['Nuclei_AreaShape_Area'].astype('float') #area
    new_col = 1-((2*np.pi*np.power((y/np.pi),0.5))/x)
    #insert columns back to the dataframe
    col_name = 'Nuclei_AreaShape_EOP'
    df_morpho_raw.insert(loc=idx, column=col_name, value=new_col)
    
    x = df_morpho_raw['CellFoot_AreaShape_Perimeter'].astype('float') #perimeter
    y = df_morpho_raw['CellFoot_AreaShape_Area'].astype('float') #area
    new_col = 1-((2*np.pi*np.power((y/np.pi),0.5))/x)
    #insert columns back to the dataframe
    col_name = 'CellFoot_AreaShape_EOP'
    df_morpho_raw.insert(loc=idx, column=col_name, value=new_col)
    #------------------------------------------------
    
    # Calculate cell spread
    #status
    print('Computing cell spreading...')
    #computes difference in percent and total btw areas of cell nuclei and cytoplasm
    
    #get raw data
    x = df_morpho_raw['CellFoot_AreaShape_Area'].astype('float') # cell area
    y = df_morpho_raw['Nuclei_AreaShape_Area'].astype('float') # nuclei area
    
    #total spreading, pxl
    new_col = x-y
    #insert columns back to the dataframe
    col_name = 'CellFoot_AreaShape_CellSpreadTot'
    df_morpho_raw.insert(loc=idx, column=col_name, value=new_col)
    
    #Rel spreading, per
    new_col = (x-y)/x
    #insert columns back to the dataframe
    col_name = 'CellFoot_AreaShape_CellSpreadPer'
    df_morpho_raw.insert(loc=idx, column=col_name, value=new_col)
    
    #attributes in the dataframe
    df_morpho_raw_Attributes_new = list(df_morpho_raw)
    
    # SAVE RESULTS
    path_res = OutPath + ResName + plate +'.pkl'
    df_morpho_raw.to_pickle(path_res)
    
    #status
    print('----------------------------------------------------------------')
    print('---> Finished. Results are saved.')
#%%"""=====END OF CALCULATIONS========================="""

#Test area














