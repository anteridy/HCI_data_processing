# -*- coding: utf-8 -*-
"""
Created on Wed Dec 16 13:14:49 2020

@author: akamnev
"""
#%%
"""==============================================
--------IMPORT PACKAGES-------------------------
================================================"""
import os #working with files and folders
import pandas as pd #pca toolkit
import numpy as np
import pickle #Create portable serialized representations of Python objects


import seaborn as sns
import matplotlib.pyplot as plt

import math

plt.style.use('default')
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
#Path to the cluster results (plate folder):
DataPath = 'D:\\Temp\\APID\\APID_5_Screen2\\Morpho\\4_QC\\Wellave_pooled\\'

PltPath = 'D:\\Temp\\APID\\APID_5_Screen2\\Morpho\\4_QC\\Corr_matrices_pooled\\' #for final plots
forcepath(PltPath)

wrapmethods = ['ave','med','std']

datasets = ['raw','cln']

#%%the loop!
for dataset in datasets:
print('Dataset = '+dataset)    
    for wrapmethod in wrapmethods:
        print('Wrap method = '+wrapmethod)
    
        '''LOAD DATA===============================================================
        '''
        DataName = 'Well_'+dataset+'_qc_'+wrapmethod+'_pool.pkl'
        #open .pkl
        path = DataPath+DataName
        with open(path, 'rb') as f:
                df_raw = pickle.load(f)
        # df_columns = list(df)
        
        #clean duplicates
        df_raw = df_raw.loc[:,~df_raw.columns.duplicated()]
        df_raw_columns = list(df_raw)    
        '''========================================================================
        '''
        
        '''GET DATASET=============================================================
        '''     
        # get dataset for correlation donors vs wells
        df = df_raw
        donors = df['Donor_N'].unique()
        # donors = ['1']
        plates = df['Metadata_Plate'].unique()
        # plates = [1]
        # wells = df['Donor_N'].unique()
        # attr = 'CellFoot_AreaShape_Area'
        attributes = ['CellFoot_AreaShape_Area',
                      'Nuclei_AreaShape_Area',
                      'Nuclei_Intensity_IntegratedIntensity_DNA',
                      'CellFoot_Intensity_MeanIntensity_Actin',
                      'Cell_number']
        '''========================================================================
        '''
        
        '''CYCLE THROUGH ATTRIBUTES================================================
        '''
        for attr in attributes:
            print(attr)
            
            '''SELF-DONOR PAIRS====================================================
            '''
            #Make corr donor-self among all wells-------------------
            df2corr_dself = pd.DataFrame()
            # index1=1
            for plate in plates:
                for donor in donors:
                    #get df
                    df_tmp = df.loc[(df['Donor_N']==donor)&(df['Metadata_Plate']==plate)]
                    df_tmp.sort_values(by=['Metadata_Row','Metadata_Column'], inplace = True)
                    
                    #get wells for this donor
                    wells = df_tmp['Metadata_Well'].unique()
                    
                    #iterate over donor 2 quads to make an array
                    array = []
                    index2 = 1
                    for well in wells:
                        val = df_tmp.loc[(df_tmp['Metadata_Well']==well)][attr].iloc[0]
                        # welln = index1+index2-1
                        well_ind = 'p'+str(plate)+'w'+str(index2)
                        df2corr_dself.at[donor,well_ind]=val
                        index2 = index2+1
                # index1 = index1+8
            # stop
            
            #Correleate
            # df2corr_dself = df_replicates.astype('float')
            corrMatrix_dself = df2corr_dself.corr()
            
            #save corr matrix
            path = PltPath+wrapmethod+'//'+attr+'//'
            forcepath(path) # make folder for the histograms
            #save as .csv
            corrMatrix_dself.to_csv(path+'corrMatrix_dself_'+dataset+'.csv')
            
            
            #plot correlation map------------------
            # sns.color_palette("vlag", as_cmap=True)
            # ax = sns.heatmap(corrMatrix_dself, annot=True, vmin=0,vmax=1,cmap='bwr')
            ax = sns.heatmap(corrMatrix_dself, vmin=0,vmax=1,cmap='bwr')
            plt.show()
            #Plot title
            # plt.title('corrMatrix_dself_Plt-'+str(plate) + '_'+attr)
            plt.title('Dself_allPts_'+dataset+'_'+wrapmethod+'_'+attr)
            
            #Save figure
            path = PltPath+wrapmethod+'//'+attr+'//'
            forcepath(path) # make folder for the histograms
            # name = 'CorrMtrx_plt'+str(plate)+'.png'
            name = 'CorrMtrx_dself_allPts_'+dataset+'.png'
            plt.savefig(path+name)
            plt.close()
            #-------------------------------------
            '''==========================================================================
            '''
            
            '''RANDOM DONOR PAIRS========================================================
            '''
            #Get random pairs of donors (non-self, non-duplicated)---------------
            donorN = np.size(donors)
            #sample twice to make sure that Npairs >= Ndonors
            rnd_pairs1 = np.random.choice(donors,size=(donorN,2), replace = True)
            rnd_pairs2 = np.random.choice(donors,size=(donorN,2), replace = True)
            rnd_pairs = np.concatenate((rnd_pairs1,rnd_pairs2), axis=0)
            
            # remove duplicated pairs
            tmp = [tuple(row)for row in rnd_pairs]
            rnd_pairs_cln1 = np.unique(tmp,axis=0)
            
            # Find self-pairs
            mask = []
            for d1, d2 in rnd_pairs_cln1:
                if d1!=d2:
                    mask.append(True)
                else:
                    mask.append(False)
                    
            #remove self-paris
            rnd_paris_cln_fin = rnd_pairs_cln1[mask,...]
            #-----------------------------------------------------------------------
        
            #save random donor pairs as csv
            path = PltPath+wrapmethod+'//'+attr+'//'
            forcepath(path) # make folder for the histograms
            name = 'RandomDonorPairs.csv'
            path2save = path+name
            pd.DataFrame(rnd_pairs_cln1).to_csv(path2save, index = False)
            # stop
            
            #Make corr btw random pairs (non-self, non-duplicated)-------------------
            df2corr_d1 = pd.DataFrame()
            df2corr_d2 = pd.DataFrame()
            # index1=1
            for plate in plates:
                for d1,d2 in rnd_paris_cln_fin:
                    #get dfs for each donor
                    df_tmp1 = df.loc[(df['Donor_N']==d1)&(df['Metadata_Plate']==plate)]
                    df_tmp2 = df.loc[(df['Donor_N']==d2)&(df['Metadata_Plate']==plate)]
                    
                    #sort wells
                    df_tmp1.sort_values(by=['Metadata_Row','Metadata_Column'], inplace = True)
                    df_tmp2.sort_values(by=['Metadata_Row','Metadata_Column'], inplace = True)
                    
                    #get wells for this donor
                    wells1 = df_tmp1['Metadata_Well'].unique()
                    wells2 = df_tmp2['Metadata_Well'].unique()
        
                    # stop
                    #iterate over donor 2 quads to make an array
                    array = []
                    # index = 1
                    for index in list(range(0,8)):
                        val1 = df_tmp1[attr].iloc[index]
                        val2 = df_tmp2[attr].iloc[index]
                        # welln = index1+index2-1
                        well_ind = 'p'+str(plate)+'w'+str(index+1)
                        donor_ind = 'd'+d1+'_d'+d2
                        df2corr_d1.at[donor_ind,well_ind]=val1
                        df2corr_d2.at[donor_ind,well_ind]=val2
            
            #Correleate
            # df2corr_dself = df_replicates.astype('float')
            df2corr_drnd = pd.concat([df2corr_d1, df2corr_d2], axis=1, keys=['d1', 'd2'])
            corrMatrix_drnd = df2corr_drnd.corr().loc['d1', 'd2']
            
            #save corr matrix
            path = PltPath+wrapmethod+'//'+attr+'//'
            forcepath(path) # make folder for the histograms
            #save as .csv
            corrMatrix_drnd.to_csv(path+'corrMatrix_drnd_'+dataset+'.csv')
            
            #plot correlation map------------------
            # sns.color_palette("vlag", as_cmap=True)
            # ax = sns.heatmap(corrMatrix_dself, annot=True, vmin=0,vmax=1,cmap='bwr')
            df2plot = corrMatrix_drnd.abs() #make values as abosolute to see high negative correaltion
            ax = sns.heatmap(corrMatrix_drnd, vmin=0,vmax=1,cmap='bwr')
            plt.show()
            #Plot title
            # plt.title('corrMatrix_dself_Plt-'+str(plate) + '_'+attr)
            plt.title('Drnd_allPts_'+dataset+'_'+wrapmethod+'_'+attr)
            
            #Save figure
            path = PltPath+wrapmethod+'//'+attr+'//'
            forcepath(path) # make folder for the histograms
            # name = 'CorrMtrx_plt'+str(plate)+'.png'
            name = 'CorrMtrx_drnd_allPts_'+dataset+'.png'
            plt.savefig(path+name)
            plt.close()
            #-------------------------------------
            '''==========================================================================
            '''
            
            # STOP
            
            '''=DISTRIBUTION OF CORRELATINS IN ALL GROUPS===============================
            '''
            #Get unique paris of correlation (w/o self)------------------------
            df_pairs = corrMatrix_dself.stack() #get unique pairs
            df_pairs = df_pairs[df_pairs.index.get_level_values(0) != df_pairs.index.get_level_values(1)] #remove self pairs
            df_pairs.index = df_pairs.index.map('_'.join) #merge into single index with '_'
            df_pairs = df_pairs.to_frame().T #transpose 1 col -> 1 col per pair
            #-------------------------------------------------------------------
            
            #GeT distributions for groups ---------------------------------------------
            quad1 = [1,2,3,4]
            quad2 = [5,6,7,8]
            df_corr_groups = pd.DataFrame()
            
            #Quad self
            vals = []
            groups = []
            df_tmp = pd.DataFrame()
            for plate in plates:
                for well1 in quad1:
                    for well2 in quad1:
                        if well1 != well2:
                            well_index ='p'+str(plate)+'w'+str(well1)+'_p'+str(plate)+'w'+str(well2)
                            val = df_pairs[well_index].iloc[0]
                            vals.append(val)
                            groups.append('Quad_self')  
            df_tmp['Pearson'] = vals
            df_tmp['Group'] = groups
            df_corr_groups = df_tmp           
            
            #Btw quads
            vals = []
            groups = []
            df_tmp = pd.DataFrame()
            for plate in plates:
                for well1 in quad1:
                    for well2 in quad2:
                        if well1 != well2:
                            well_index ='p'+str(plate)+'w'+str(well1)+'_p'+str(plate)+'w'+str(well2)
                            val = df_pairs[well_index].iloc[0]
                            vals.append(val)
                            groups.append('Quad_Quad')  
            df_tmp['Pearson'] = vals
            df_tmp['Group'] = groups
            
            df_corr_groups = pd.concat([df_corr_groups,df_tmp], axis=0, ignore_index=True)
                   
            #Random pairs-------------------------------------------------
            #Get unique paris of correlation (w/o self)------------------------
            df_pairs_rnd = corrMatrix_drnd.stack() #get unique pairs
            df_pairs_rnd = df_pairs_rnd[df_pairs_rnd.index.get_level_values(0) != df_pairs_rnd.index.get_level_values(1)] #remove self pairs
            df_pairs_rnd.index = df_pairs_rnd.index.map('_'.join) #merge into single index with '_'
            df_pairs_rnd = df_pairs_rnd.to_frame().T #transpose 1 col -> 1 col per pair
            #-------------------------------------------------------------------
            #Get distributions for groups ---------------------------------------------  
            #Quad self
            vals = []
            groups = []
            df_tmp = pd.DataFrame()
            for plate in plates:
                for well1 in quad1:
                    for well2 in quad1:
                        if well1 != well2:
                            well_index ='p'+str(plate)+'w'+str(well1)+'_p'+str(plate)+'w'+str(well2)
                            val = df_pairs_rnd[well_index].iloc[0]
                            vals.append(val)
                            groups.append('Quad_self_Random')  
            df_tmp['Pearson'] = vals
            df_tmp['Group'] = groups
            df_corr_groups = pd.concat([df_corr_groups,df_tmp], axis=0, ignore_index=True)         
            
            #Btw quads
            vals = []
            groups = []
            df_tmp = pd.DataFrame()
            for plate in plates:
                for well1 in quad1:
                    for well2 in quad2:
                        if well1 != well2:
                            well_index ='p'+str(plate)+'w'+str(well1)+'_p'+str(plate)+'w'+str(well2)
                            val = df_pairs_rnd[well_index].iloc[0]
                            vals.append(val)
                            groups.append('Quad_Quad_Random')  
            df_tmp['Pearson'] = vals
            df_tmp['Group'] = groups
            
            df_corr_groups = pd.concat([df_corr_groups,df_tmp], axis=0, ignore_index=True)
            #-------------------------------------------------------------
                 
            # Plot
            df2plot = df_corr_groups
            x_group = 'Group'
            attribute = 'Pearson'
            plt.clf() #clear figure
            #make violing plot for the attribute
            ax = sns.violinplot(x = x_group, y = attribute, data = df2plot)
            plt.title('Corr_'+dataset+'_'+wrapmethod+'_'+attribute+'_allPts_'+attr)
            
            #save plot   
            path = PltPath+wrapmethod+'//'+attr+'//'
            forcepath(path) # make folder for the histograms
            name = 'Corr_'+dataset+'_'+attribute+'_allPts_'+attr+'.png'
            plt.savefig(path+name)
            plt.close()
            '''==========================================================================
            '''

#%%END
print('FINISHED--------------------------------------------------------')














