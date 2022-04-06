# -*- coding: utf-8 -*-
"""
Created on Wed Dec  2 13:40:21 2020

@author: akamnev

Based on original scrip from Michael Caldera "Extract data from (Anton) CellProfiler
 file" -> "1_Read_Data" Jupiter script.
"""
#%%
'''=====IMPORT======================================='''
#import all python libraries
import networkx as nx #code is build originally for v2.2, >=2.5 synthaxis changed a bit
import numpy as np
from os import listdir
from os.path import isfile, join
from random import choice
from matplotlib import pylab as plt
from collections import Counter
import os
import pandas as pd
import pickle
from sklearn.preprocessing import MinMaxScaler
from scipy import stats
#time stamp
from datetime import datetime
#to wrap nested loops
import itertools

#check networkx library version
nx_version = nx.__version__
print('nx version = '+ nx_version)
print('OK to proceed')
#abort if version is not right
if nx_version != '2.2':
    print('nx version = '+ nx_version+' must be 2.2 for script to work!')
    print ('install as <pip install networkx==2.2>')
    wrong_nx_version

#to change
#pip install networkx==2.2

#%%
'''=====STANDARD FUNCTIONS======================================='''
def mad(data, axis=None):
    return np.mean(np.absolute(data - np.mean(data, axis)), axis)


def ensure_dir(file_path):
    '''
    Function to ensure a file path exists, else creates the path

    :param file_path:
    :return:
    '''
    directory = os.path.dirname(file_path)
    if not os.path.exists(directory):
        os.makedirs(directory)


def reject_outliers_2(data, m = 3.):

    d = np.abs(data - np.median(data))

    mdev = np.median(d)
    s = d/(mdev if mdev else 1.)
    return [data[i] for i in range(0,len(data)) if s[i] < m]


#Effect size
def cohen_d(x, y):
    nx = len(x)
    ny = len(y)
    dof = nx + ny - 2
    return (np.mean(x) - np.mean(y)) / np.sqrt(((nx - 1) * np.std(x, ddof=1) ** 2 + (ny - 1) * np.std(y, ddof=1) ** 2) / dof)

def calculate_ZFactor(x,y):
    '''
    print '---'
    print np.std(drug)
    print np.std(dmso)
    print np.mean(drug)
    print np.mean(dmso)
    print '---'
    '''

    return 1-((3*np.std(x)+3*np.std(y))/(abs(np.mean(x) -np.mean(y))))

'''=====CORRELATION FUNCTIONS======================================='''
#Correlation Functions
def create_correlation_Network(data, effectSizes):

    #get all features to inspect
    features = list(data.columns)

    # create output graph
    G = nx.Graph()
    print ('\n\nCalculating PID Correlation...')
    #for f in feature_names:
    for f in features:
        #print f
        G.add_node(f)
        #for f2 in feature_names:
        for f2 in features:
            if f > f2:

                feature1 = data[f]
                feature2 = data[f2]

                if len(feature1) == 0 or len(feature2) == 0:
                    continue

                cor = np.corrcoef(feature1, feature2)[0, 1]

                if abs(cor) > 0.6:
                    G.add_edge(f,f2)
                    #Use G.nodes instead of G.node for never version of networks (>2.4)
                    G.node[f]['Effect'] = effectSizes[f]
                    G.node[f2]['Effect'] = effectSizes[f2]


    # Save correlation network for manual inspection
    #if saveNetwork:
    #    ensure_dir('../../results/Remove_Correlation/correlationnetwork'+filename+'.gml')
    #    nx.write_gml(G, '../../results/Remove_Correlation/correlationnetwork'+filename+'.gml')

    return G

def check_edge_between(min_nod,g):
    '''
    Find minimal amount of nodes to remove, so that min_nodes are not neighbors anymore

    :param min_nod: list of nodes to check
    :param g: connected component subgraph
    :return: list of nodes to remove so that nodes in min_nod are not connected
    '''

    #extract all pairwise edges between nodes in min_nod
    has_edge = []
    for min1 in min_nod:
        for min2 in min_nod:
            if min1 > min2:
                if g.has_edge(min1, min2):
                    has_edge.append(min1)
                    has_edge.append(min2)

    #create Counter instance
    data = Counter(has_edge)
    #get values of which node occured how often edges
    #e.g. [2,2,1,1,1,1] = two nodes are connected with two in min_nodes while 4 nodes are connected only with one
    freq_list = data.values()

    #if freq_list == 0, all nodes are separated from each other
    if len(freq_list) == 0:
        return []

    #find max edges of on of the nodes (e.g. example above would be 2)
    max_cnt = max(freq_list)

    #get all nodes that are involved in max_cnt (e.g. the first two nodes from the freq_list)
    total = list(freq_list).count(max_cnt)

    #if all nodes are equaly e.g. [1,1,1], it does'nt bother which one to remove, choose randomly one (so remove the other two)
    if total == len(freq_list):
        max_val = 0
        max_node = ''
        for node in min_nod:
            if g.node[node]['Effect'] > max_val:
                max_val = g.node[node]['Effect']
                max_node = node

        keep = max_node
        #keep = choice(min_nod)
        copylist = list(min_nod)
        copylist.remove(keep)
        return copylist

    #Return these nodes for removal (first two from example above)
    most_common = data.most_common(total)
    return [elem[0] for elem in most_common]

def analyse_component(g,draw=False):
    '''
    Main function for max fragmentation of subgraphs.
    Takes a graph (origins from bigger network as connected component)
    Slowly fragmentises it by removing best suited nodes

    :param g: connected component subgraph
    :param draw: True if there should be a step by step drawing output
    :return: number of
    '''
    #contains max number nodes for current component
    tmp_keep = []

    #fragmentise connected component subgraph until all nodes fragmented
    while len(g.nodes()) > 0:

        #draw option
        if draw == True:
            nx.draw_networkx(g, pos=nx.spring_layout(g), with_labels=False)
            plt.draw()  # pyplot draw()
            plt.show()

        #list for nodes that need to be removed in each iteration
        #Contains: Selected Nodes (find in tmp_keep), as well as their neighbors
        nodes_to_remove = set()

        #if (remaing) component is only two nodes ==>  A--B; take randomly one of the two
        if len(g.nodes()) == 2 and len(g.edges()) == 1:
            two_nodes = list(g.nodes())

            if g.node[two_nodes[0]]['Effect'] > g.node[two_nodes[1]]['Effect']:
                tmp_keep.append(two_nodes[0])
            else:
                tmp_keep.append(two_nodes[1])

            #purely random choice of which node to keep
            #rand_node = choice(g.nodes())
            #tmp_keep.append(rand_node)

            nodes_to_remove.add(list(g.nodes())[0])
            nodes_to_remove.add(list(g.nodes())[1])

        #if bigger than only two connected nodes
        else:

            #get node degrees
            degrees_tmp = g.degree()

            degrees = {}
            for d in degrees_tmp:
                degrees[d[0]] = d[1]


            #find terminal nodes (= degree 1)
            terminal_nodes = [x for x in degrees if degrees[x] == 1]

            #if subgraph still has terminal nodes, choose these
            if len(terminal_nodes) > 0:
                for tn in terminal_nodes:
                    tmp_keep.append(tn)
                    nodes_to_remove.add(list(g.edges(tn))[0][1])
                    nodes_to_remove.add(tn)

            #if no terminal nodes exist
            else:

                #Check if there are nodes with higher degree than other
                #if all degrees uniformly it's for example a triangle, rectangle etc. (circularity)
                if all(x == list(degrees.values())[0] for x in list(degrees.values())) == False:
                    #example for some nodes with lower/higher degree than others
                    # A-B
                    # |\|  ==> in this case the algorithm should pick B and C (other alternative would be only A or D)
                    # C-D

                    #extract smalles degree
                    min_degree = min(degrees.values())

                    #get nodes with this smallest degree
                    min_nodes = [x for x in degrees if degrees[x] == min_degree]

                    #check if these nodes with smallest degree are somehow neighbors
                    #e.g. two rectangles (4 nodes) connected by a middle node
                    # ==> always the three "outer" rectangle nodes would have degree 2 (togher with the middle one connecting the two rectangles)
                    while True:
                        #remove the minimum amount of nodes, so all selected "min_nodes" are no neighbors anymore
                        node_edge_remove = check_edge_between(min_nodes, g)
                        if len(node_edge_remove) == 0:
                            break
                        for node in node_edge_remove:
                            min_nodes.remove(node)

                    #Save the Min_nodes to tmp_keep and add them to nodes_to_remove (togher with their neighbors)
                    for mn in min_nodes:
                        tmp_keep.append(mn)
                        nodes_to_remove.add(mn)
                        edges = g.edges(mn)
                        for edge in edges:
                            nodes_to_remove.add(edge[1])

                #if all degrees are uniformly, meaning you have a triangle, rectangle, fully connected graph
                #e.g.:
                # A-B
                # | |  ==> e.g. first pick A (randomly); remove A + C and D (neighbors); in next
                # C-D           iteration there is only D left (will be picked)
                else:
                    #randomly choose a single node (all nodes equally anyway)
                    max_val = 0
                    max_node = ''
                    for node in g.nodes():
                        if g.node[node]['Effect'] > max_val:
                            max_val = g.node[node]['Effect']
                            max_node = node

                    rand_node = max_node

                    #rand_node = choice(g.nodes())
                    #add this random nood to tmp_keep and again remove him together with the neighbors
                    tmp_keep.append(rand_node)
                    nodes_to_remove.add(rand_node)
                    edges = g.edges(rand_node)
                    for edge in edges:
                        nodes_to_remove.add(edge[1])

        #Remove nodes from current subgraph
        for ntr in nodes_to_remove:
            g.remove_node(ntr)

        if draw:
            print (tmp_keep)

    return tmp_keep

def forcepath(path):
    """checks if directory exists at path <path> and, if not, makes one"""
    if not os.path.isdir(path):
        os.makedirs(path)

#%%
'''=====iNPUTS======================================='''
#Path to the raw data
DataPath = 'D:\\Temp\\APID\\APID_5_Screen2\\Morpho\\'
# Folder = '6_PreProcessed_data\\Sample_wrap_liveCD8\\'0
Folder = '6_PreProcessed_data\\ImgClassifier\\Sample_wrap_liveCD8\\'
File = 'AveSample_S20_plt_'

ResFolder ='7_Filtered_data\\Sample_wrap_liveCD8\\'
# ResFolder ='7_Filtered_data\\ImgClassifier\\Sample_wrap_liveCD8\\'
# ResFolder ='7_Filtered_data\\ImgClassifier\\Sample_wrap_liveCD8_forLoic\\'


#Metadata files (full path to the file)--------------------------
#Well metadata
# MetadataUnwPIDs_path = 'D:\\Temp\\APID\\APID_5_Screen2\\Morpho\\7_Filtered_data\\Unwanted_PIDs.csv'
unw_pids = ['DOCK11','PAK2','THEMIS2']
# unw_pids = []
unw_donors = ['PID_2599 (CD4/8 mix)', 'ND_TL115 (w7)']

numdatacol = 21 #first col with numerical morpho data

#Correlation removal
# D_min_array = [1,2] #min cohen's distance for feature to have to be kept
D_min = 1

#Quick test?
test = False #Open only first x rows
# test = True #Open only first x rows
rown = 300

plates = [1]


#path to pinned features:
path_pinned = '4_Metadata\\PinnedFeatures.csv'

#path to core features:
path_core = '4_Metadata\\CoreFeatures.csv'

   
for plate in plates:
    #%% Plate processing loop
    '''=====LOAD DATA======================================='''
    #status
    print('Loading data...')
    
    #load pinned features
    path = DataPath + path_pinned
    pinned_features = pd.read_csv(path)
    pinned_features_list = pinned_features['Pinned_features'].tolist()
    # data_cols = list(data)
    
    #load pinned features
    path = DataPath + path_core
    core_features = pd.read_csv(path)
    core_features_list = core_features['Feature'].tolist()
    # data_cols = list(data)
    
    # stop
    #load main data
    path = DataPath+Folder+File+str(plate)+'.pkl'
    # data = pandas.read_csv(DataPath+DataName)
    with open(path, 'rb') as f:
            data_raw = pickle.load(f)
    data_raw_columns = list(data_raw)
    # stop
    #Clean data-------------------------------------------------------------
    #Drop unwanted unwanted PIDs
    df = data_raw
    pid_ind = df['Genotype'].isin(unw_pids)
    data_tmp = df[~pid_ind]
    #Drop unwanted donors
    df = data_tmp
    don_ind = df['Patient_ID'].isin(unw_donors)
    data_raw_cln = df[~don_ind]
    
    #Drop unwanted features (Neighbors measurements)
    features2remove = [s for s in data_raw_columns if "Neighbors" in s]
    data_raw_cln.drop(features2remove,axis=1,inplace=True)
    features2remove = [s for s in data_raw_columns if "Cell_number" in s]
    data_raw_cln.drop(features2remove,axis=1,inplace=True)
    #-------------------------------------------------------------------------
    
    
    #Reduce data size for testing
    if test:
        data_pooled = data_raw_cln.iloc[:rown] #[!]for testing only
    else:
        data_pooled = data_raw_cln
    
    
    #make path to save data
    ResPath = DataPath+ResFolder+'\\Plt_'+str(plate)+'\\'
    forcepath(ResPath)
    
    #status
    print('')
    print('==============================================')
    # datetime object containing current date and time
    now = datetime.now()
    print("start =", now)
    print('Plate '+str(plate)+'...')
    
    #get df for the plate
    data = data_pooled.loc[data_pooled['Metadata_Plate']==plate]
    
    
    #Get list of morpho-feature names
    features = list(data.columns[numdatacol:])
    # print(features[0:4])
    # print('Number of features: %d' %len(features))
    features.sort()
    
    # This dataset has 9 types of mutations (including normal donor)
    # Each of the 9 mutations can have several donors (typically two)
    ##
    donor_groups = list(set(data['Patient_ID']))
    
    pid_groups = list(set(data['Genotype']))
    # print (donor_groups)
    # data.head(5)
    # stop
    #%%
    '''=====MAIN WORKFLOW=======================================
    1. Normalize everything to normal donor (mean of all normal donors)
    2. Filter features based on changes between ND and PID
    3. Filter correlating features
    4. Create Final ResultFile
    '''
    
    '''---1-Normalize-----------------------------------------'''
    ###
    # DECIDE NORMALIZATION METHOD
    # POC = Max/Min normalization using 0.5 and 99.5 percentiles
    # STD = calculates everything as Zscore given a reference (e.g. ND)
    ###
    
    #status
    print(' ')
    print('------------------------------------------------------')
    print('Normalizing to ND -> STD method ...')
    
    
    #make a normalized copy of the original dataframe
    data_normalized = data.copy()
    
    #below 4 digits only random noise
    '''AKM - Do I need this?'''
    data_normalized = data_normalized.round(5)
    
    ND_values = data_normalized.loc[(data_normalized['Genotype'] == 'ND')]
    
    #Go through all features
    for f in features:
        #f = 'Nuclei_Texture_Contrast_Actin_12_01'
    
        #get feature specific values; remove nan values
        ND_feature_val = [x for x in ND_values[f] if str(x) != 'nan']
        #print(ND_feature_val)
    
    
        #if less than 10 unique values exists == no real variability in the feature => feature is useless
        # These features just cause big problem --> there a plenty of fish in the sea ....
        if len(set(data_normalized[f])) < 10:
            data_normalized[f] = 0
    
        else:
            #Reject outlier
            '''AKM - why this is not used by Michi?'''
            #ND_feature_val = reject_outliers_2(list(ND_feature_val))
    
            #calculate mean and std
            mean_ND_val = np.mean(ND_feature_val)
            std_ND_val = np.std(ND_feature_val)
    
            #median_ND_val = np.median(ND_feature_val)
            #mad_ND_val = mad(ND_feature_val)
    
            data_normalized[f] = data_normalized[f].subtract(mean_ND_val)
            data_normalized[f] = data_normalized[f].div(std_ND_val)
    
            #ensure that columns are floats
            '''AKM: if numbers are object type then nx will get confused and throw
            error <float has no attribute shape>'''
            data_normalized[f] = pd.to_numeric(data_normalized[f])
    
    
    #drop columns that have only 'nan' values
    data_normalized = data_normalized.dropna(axis=1,how='all')
    
    # data_normalized.head()
    # stop
    #%%
    '''---2-Effect size-----------------------------------------
    -Get the ND values for a certain features
    -Calculate cohen's D (effect) size between ND and all types of mutations (for POC)
    -Use STD (ZScores) f
    '''
    #status
    print('------------------------------------------------------')
    print('Calculating effect sizes...')
    
    # CHECK FOR EFFECT SIZE
    ########
    #PID groups to check ND against
    PID_Donor_Groups = [x for x in donor_groups if 'ND' not in x]
    PID_Groups = [x for x in pid_groups if 'ND' not in x]
    # stop
    
    #result files
    feature_results_effect = {} #max effect size for all pairs
    feature_results_effects_per_Group = {} #effect size per group
    feature_results_effects_per_PID = {} #effect size per PID
    
    #some features got kicked out in the normalization step
    keep_features = [] #list of features processed by filter (and kept)
    
    #if STD you already normalized to ND
    
    # print("Effect Size: STD")
    
    #Go through all features
    for f in features:
        if f not in data_normalized:
            continue
    
        #remove useless features
        if sum(data_normalized[f]) == 0:
            continue
    
        #keep name of the selected feature
        keep_features.append(f)
    
        #create per Group results
        feature_results_effects_per_Group[f] = {}
        feature_results_effects_per_PID[f] = {}
    
        #tmp will save the separate donor results
        tmp = []
    
        #go through donor results
        for pid in PID_Donor_Groups:
    
            #get the PID values w/o NaNs
            pid_values = [x for x in data_normalized.loc[(data_normalized['Patient_ID'] == pid)][f].values if str(x) != 'nan']
    
            #remove clear outliers
            '''again michi doesn't use outliers removal - Why?'''
            #valid_values = reject_outliers_2(pid_values)
    
            #effect size = mean
            effect_size = abs(np.median(pid_values))
    
            #add to temporary results (use absolute as both directions equaly interesting)
            tmp.append(effect_size)
            feature_results_effects_per_Group[f][pid] = effect_size
    
    
        feature_results_effect[f] = max(tmp)
        
        #go through donor results
        tmp = []
        for pid in PID_Groups:
    
            #get the PID values w/o NaNs
            pid_values = [x for x in data_normalized.loc[(data_normalized['Genotype'] == pid)][f].values if str(x) != 'nan']
    
            #remove clear outliers
            '''again michi doesn't use outliers removal - Why?'''
            #valid_values = reject_outliers_2(pid_values)
    
            #effect size = mean
            effect_size = abs(np.median(pid_values))
    
            #add to temporary results (use absolute as both directions equaly interesting)
            tmp.append(effect_size)
            feature_results_effects_per_PID[f][pid] = effect_size
    
    
    features = keep_features
    
    
    # Save Output
    ########
    #status
    print('Saving effect size results...')
    
    path = ResPath
    ensure_dir(path)
    
    #make csv files for saving data
    fp = open(path+'EffectSize.csv','w')
    fp.write('Feature,Max_CohenD\n')
    
    fp2 = open(path+'EffectSize_Individual.csv','w')
    fp2.write('Feature,'+','.join(PID_Donor_Groups)+'\n')
    
    fp3 = open(path+'EffectSize_IndividualPIDs.csv','w')
    fp3.write('Feature,'+','.join(PID_Groups)+'\n')
    
    
    #Go through all features & write data to the result files
    for f in features:
        #write ouptut
        fp.write(f+','+str(feature_results_effect[f])+'\n')
    
        fp2.write(f)
        for pid in PID_Donor_Groups:
            fp2.write(','+str(feature_results_effects_per_Group[f][pid]))
        fp2.write('\n')
        
        fp3.write(f)
        for pid in PID_Groups:
            fp3.write(','+str(feature_results_effects_per_PID[f][pid]))
        fp3.write('\n')
    
    
    #close files
    fp.close()
    fp2.close()
    fp3.close()
    
    # stop
    #%%
    '''---3-REMOVE CORRELATING FEATURES-----------------------------------------
    '''
    print('')
    print('------------------------------------------------------')
    print('Removing correlating features (CohenD <'+str(D_min)+')...')
    
    #Dictionary to collect filtered features
    final_features_after_correlation = {}
    
    print('Total number of features: %d' %len(features))
    #extract those features that have at least a cohen's D of 1 (or as set)
    significant_features = []
    for f in features:
        if feature_results_effect[f] > D_min:
            significant_features.append(f)
    print('Number of significant features (by effect size): %d' %len(significant_features))
    
    final_features_after_correlation = []
    
    #create correlation network
    g = create_correlation_Network(data_normalized[significant_features],feature_results_effect)
    
    #save correlation network
    path = ResPath
    nx.write_gml(g, path+'CorrelationNetwork_CD'+str(D_min)+'.gml')
    
    #Extract connected components
    #v2.0
    components = nx.connected_component_subgraphs(g) #depreciated in nx. version >2.4
    #v2.5
    # components = (g.subgraph(c) for c in nx.connected_components(g))
    
    #keep contains the max amount of nodes that are never connected
    keep = []
    
    #Go through the components
    for comp in components:
        #if single node, add to keep (anyway not connected to anything)
        if len(comp.nodes()) == 1:
            keep.append(list(comp.nodes())[0])
        #if not single node check maximum fragmentation
        else:
            '''Here graph appears to be frozen (in v>2.4 by default?)'''
            keep = keep + analyse_component(comp,False)
    
    keep.sort()
    
    final_features_after_correlation = keep
    
    path = ResPath
    fp = open(path+'NoCorrelation_CD'+str(D_min)+'.csv','w')
    for feature in keep:
        fp.write(feature+'\n')
    fp.close()
    print ('Number of Features: %d' %len(keep))
    
    # stop
    #
    #%%
    '''---4-FINAL DATA FILE TO SAVE-----------------------------------------
    '''
    print ('Saving final results...')
    #define metadata columns
    metadata = ['Genotype','Patient_ID','Donor_N','Group','Metadata_Plate',
                         'Metadata_Column','Metadata_Row','Inplate_Replica','ImageNumber',
                         'Coating','IL-2','Marker']
    
    #Use only non correlating Features
    filtered_data = metadata + final_features_after_correlation
    # filtered_data.extend(final_features_after_correlation)
    
    #Use all feature that have some form of effect
    #extract those features that have at least a cohen's D of 1
    '''
    significant_features = []
    for f in features:
        if feature_results_effect[f] > 1:
            significant_features.append(f)
    important_columns.extend(significant_features)
    '''
    
    #get value specific for dataset
    #DS_Specific_Values = APID_normalized.loc[(APID_normalized['Metadata_aCD3_conc'] == ds)].copy()
    
    #get only nessecary columns
    DS_Specific_Values = data_normalized[filtered_data]
    #save as csv
    path = ResPath
    DS_Specific_Values.to_csv(path+'FinalResultFile_CD'+str(D_min)+'.csv')
    
    #%% Save core & pinned features
    metadata = ['Genotype','Patient_ID','Donor_N','Group','Metadata_Plate',
                         'Metadata_Column','Metadata_Row','Inplate_Replica','ImageNumber',
                         'Coating','IL-2','Marker']
    
    # pinned_data = metadata + pinned_features_list
    core_data = metadata + core_features_list
    
    # #get only nessecary columns - pinned features
    # DS_Specific_Values_pinned = data_normalized[pinned_data]
    # #save as csv
    # path = ResPath
    # DS_Specific_Values_pinned.to_csv(path+'FinalResultFile_pinned.csv')
    
    #get only nessecary columns - core features
    DS_Specific_Values_core = data_normalized[core_data]
    #save as csv
    path = ResPath
    DS_Specific_Values_core.to_csv(path+'FinalResultFile_core.csv')
    
    #Get average per PID
    df = DS_Specific_Values_core
    DS_Specific_Values_core_mean = df.groupby(['Genotype']).mean()
    #save as csv
    path = ResPath
    DS_Specific_Values_core_mean.to_csv(path+'FeatureMeanPerPID_core.csv')
    
    print("DONE with data preparation for plate "+str(plate))
    
    #%%end
    #status
    print('DONE----------------------------------')

#%% END


# a = pinned_features['Pinned_features'].tolist()
# b = metadata + a








