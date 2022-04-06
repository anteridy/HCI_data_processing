"""==============================================
#------------INTRO--------------------------------
==================================================
This script measures covenrts raw per-cell morpho-dataset to
    per-image dataset (average per image)

[19.06.25] - do not remove Nan col/rows before averaging, calculate mean ignoring NaNs instead

1. Works on SC annotated dataset (cropped to wanted features)
2. Projects on basis of donor to mean, median and std
3. merges plates to single pickle
4. saves single pkl file


#-------------------------------------------------
==================================================
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

#plots
import matplotlib.pyplot as plt #main plotting engine
import seaborn as sns           # makes plots look nicer



#------FUNCTIONS-----------------------------------
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

def forcepath(path):
    """checks if directory exists at path <path> and, if not, makes one"""
    if not os.path.isdir(path):
        os.makedirs(path)
#----END OF FUNCTIONS-----------------------------
"""-------------------------------------------------
==================================================
"""

#%%
"""==============================================
--------INPUT - Raw dataset(s)-------------------------
================================================"""
#PATHS-----------------------------------------------
DataPath = 'D:\\Temp\\APID\\APID_5_Screen2\\Morpho\\'
Folder = '7_Filtered_data\\Sample_wrap_liveCD8\\Plt_'
File = 'FinalResultFile_CD1.csv'

#Where to write results
plates = [4,5,6,10,11,12]

# #save res to:
# ResPath = DataPath+Folder+str(plate)+'\\'
#----------------------------------------------------------

#Other inputs:-----------------------------------------------
lastMetCol = 12 #last column of metadata in the table to split df
#------------------------------------------------------------

#Data to remove
unw_pids = ['DOCK11','PAK2','THEMIS2']
#Donors to remove
unw_donors = ['PID_2599 (CD4/8 mix)', 'ND_TL115 (w7)']

for plate in plates:    
    #save res to:
    ResPath = DataPath+Folder+str(plate)+'\\DonAve_filt\\'
    forcepath(ResPath)
    
    
    #%% Load data
    # Load data-------------------------------------
    #status
    print('Loading data...')
    
    #Load pkl file with morpho data
    path = DataPath+Folder+str(plate)+'\\'+File
    data_raw = pd.read_csv(path)
    data_cols = list(data_raw)
    
    #Clean data---------------------------------
    #Drop unwanted unwanted PIDs
    df = data_raw
    pid_ind = df['Genotype'].isin(unw_pids)
    data_raw_cln = df[~pid_ind]
    df = data_raw_cln
    don_ind = df['Patient_ID'].isin(unw_donors)
    data_fin = df[~don_ind]  
    #-------------------------------------------
    # stop
       
    #%% Donor level aggregation 
    merged = []
    # Calculate average per donor-------------------------------
    #status
    print('Calculating Donor average, median & stds ...')
    
    df = data_fin
    
    #Get donor IDs
    donors = np.sort(df['Donor_N'].unique())
    
    #location of metadata attributes
    metadata_cols = list(range(1,lastMetCol+1))
    feature_cols = list(range(lastMetCol+1,df.shape[1]))
    
    #location of numberical attributes to average
    # result_cols = features2keep
    
    #Prepare workspace
    means = []
    medians = []
    stds = []
    #cycle through donors range
    for donor in donors:
        #get subset which belongs to the image i
        df_donor = df.loc[(df['Donor_N']== donor)]
        
        #get metadata vector
        df_seed = df_donor.iloc[0,metadata_cols]
        # stop
        #extract numberical attributes
        # df_tmp = df_donor.iloc[:,result_cols]
        df_tmp = df_donor.iloc[:,feature_cols]
        # stop
        
        #clean NaNs & average along columns
        #df_tmp2 = cleanNaNs(df_tmp).mean(axis=0)
        
        #Average ingoring nan and null
        df_tmp2_mean = df_tmp.mean(axis=0, skipna=True)
        df_tmp2_median = df_tmp.median(axis=0, skipna=True)
        df_tmp2_std = df_tmp.std(axis=0, skipna=True)
        
        #merge sets
        df_fin_mean = pd.concat([df_seed, df_tmp2_mean], axis=0, ignore_index=False).transpose()
        df_fin_median = pd.concat([df_seed, df_tmp2_median], axis=0, ignore_index=False).transpose()
        df_fin_std = pd.concat([df_seed, df_tmp2_std], axis=0, ignore_index=False).transpose()
        
        #add number of cells in the FOV
        #df_fin['Cell_number'] = cleanNaNs(df_tmp).shape[0]
        # df_fin_mean['Cell_number'] = df_tmp.shape[0]
        # df_fin_median['Cell_number'] = df_tmp.shape[0]
        # df_fin_std['Cell_number'] = df_tmp.shape[0]
        
        #append
        means.append(df_fin_mean)
        medians.append(df_fin_median)
        stds.append(df_fin_std)
    
    #merge results to single dataframe & transpose
    df_means = pd.concat(means, axis=1, ignore_index=True).T
    df_medians = pd.concat(medians, axis=1, ignore_index=True).T
    df_stds = pd.concat(stds, axis=1, ignore_index=True).T
        
        
    """===============================================
    ------SAVE RESULT------------------------------------
    ==============================================="""
    #status
    #print('==========================================================')
    print('Saving results...')
    #save result - means
    path = ResPath+'DonorAve.csv'
    # df_means_merged.to_pickle(path)
    df_means.to_csv(path, index=False)
    
    #save result - medians
    path = ResPath +'DonorMed.csv'
    df_medians.to_csv(path, index=False)
    
    #save result - stds
    path = ResPath +'DonorStd.csv'
    df_stds.to_csv(path, index=False)
          
    #status
    print('==========================================================')
    print('---> DONE.')

#%%
"""===============================================
--------Testing area--------------------"""

# is_NaN = df_raw.isnull()
# row_has_NaN = is_NaN.any(axis=1)
# NaN_N_rows = row_has_NaN[0].value_counts()['True']
# col_has_NaN = is_NaN.any(axis=0)














