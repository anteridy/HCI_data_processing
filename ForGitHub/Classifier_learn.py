# -*- coding: utf-8 -*-
"""
Created on Thu Feb 25 16:31:03 2021

@author: akamnev

How to build ROC curves with Sklearn in python:
    https://stackabuse.com/understanding-roc-curves-with-python/
    
How to save learned model 
https://machinelearningmastery.com/save-load-machine-learning-models-python-scikit-learn/


"""
#%%Import 

#OS functions
import os #working with files and folders

#Dataframes and numberical operations
import numpy as np
import pickle #Create portable serialized representations of Python objects
import pandas as pd #pandas dataframe toolkit for large datasets

#machine learning
import sklearn

# roc curve and auc score
from sklearn.datasets import make_classification
from sklearn.neighbors import KNeighborsClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn import svm
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_curve
from sklearn.metrics import roc_auc_score

#plots
import matplotlib.pyplot as plt #main plotting engine
import seaborn as sns           # makes plots look nicer

#%% Functions

def forcepath(path):
    """checks if directory exists at path <path> and, if not, makes one"""
    if not os.path.isdir(path):
        os.makedirs(path)
        
def plot_roc_curve(fpr, tpr):
    plt.plot(fpr, tpr, color='orange', label='ROC')
    plt.plot([0, 1], [0, 1], color='darkblue', linestyle='--')
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('Receiver Operating Characteristic (ROC) Curve')
    plt.legend()
    plt.show()
        
#%% Inputs

DataPath = 'D:\\Temp\\APID\\APID_5_Screen2\\Morpho\\'
DataFolder = '6_ImgClassifier\\CellImages\\Data_ann\\'
ResFolder = '6_ImgClassifier\\Model\\'

training_sets = [1,2]
test_sets = range(1,11)

#Columns to remove before saving-------------------------
FeatFolder = '4_Metadata\\'
FeatName_wanted = 'CoreFeatures.csv'
#------------------------------------------


#%% Load data

#Load training sets
tmp = []
for setN in training_sets:
    path = DataPath+DataFolder+'Sample_'+str(setN)+'_ann.csv'
    data = pd.read_csv(path)
    tmp.append(data)
data_train = pd.concat(tmp,axis=0)

#Load testing sets
tmp = []
for setN in test_sets:
    path = DataPath+DataFolder+'Sample_'+str(setN)+'_ann.csv'
    data = pd.read_csv(path)
    tmp.append(data)
data_test = pd.concat(tmp,axis=0)    
    
data_attr = list(data_train)

#load wanted features
foldpath = DataPath+FeatFolder
path = foldpath+FeatName_wanted
wanted_cols = pd.read_csv(path, index_col=None, header=0)['Feature'].astype(str)

#%% Classifiation
# from sklearn import svm

#get classes
# X = data_train[wanted_cols]
# y = data_train['Trash']

#%% Plot values
# df = data['Trash','CellFoot_AreaShape_Area','CellFoot_Intensity_MeanIntensity_Actin']
# df = data_train[['Trash','CellFoot_AreaShape_Area','CellFoot_Intensity_MeanIntensity_Actin']]

# # stop
# # plot training data
# df.plot.scatter(x='CellFoot_AreaShape_Area',
#                       y='CellFoot_Intensity_MeanIntensity_Actin',
#                       c='Trash',
#                       colormap='viridis')
# plt.title('Training data')

#%%
#build classifier
# clf = svm.SVC()
# clf.fit(X,y)

# # Predict
# X_test = data_test[wanted_cols]
# y_test = data_test['Trash']

# a = clf.predict(X_test)

# # Compare
# mask = y_test==a

# b = mask.value_counts()


#%% ROC curve

# data_X, class_label = make_classification(n_samples=1000, n_classes=2, weights=[1,1], random_state=1)
# get data
data_X = data_test[wanted_cols]
class_label = data_test['Trash']

# Split dataset
trainX, testX, trainy, testy = train_test_split(data_X, class_label, test_size=0.3, random_state=1)

# Fit a model on the train data.
model = RandomForestClassifier()
# model = KNeighborsClassifier()
# model = svm.SVC(probability=True)
model.fit(trainX, trainy)

# redict probabilities for the test data.
probs = model.predict_proba(testX)

# Keep Probabilities of the positive class only.
probs = probs[:, 1]

# Compute the AUC Score.
auc = roc_auc_score(testy, probs)
print('AUC: %.2f' % auc)

# Get the ROC Curve.
fpr, tpr, thresholds = roc_curve(testy, probs)

# Plot ROC Curve using our defined function
plot_roc_curve(fpr, tpr)

#%% save the model to disk
# filename = 'finalized_model.sav'
path = DataPath + ResFolder
forcepath(path)
path1 = path + 'TrainedModel.sav'
pickle.dump(model, open(path1, 'wb'))

#%% load the model from disk
loaded_model = pickle.load(open(path1, 'rb'))
# result = loaded_model.score(X_test, Y_test)
# print(result)

#%% END

















