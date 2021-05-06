# -*- coding: utf-8 -*-
"""
Created on Mon Apr 26 11:42:08 2021

@author: Mudhurika
"""

#making data for glycosylation
a_file = open("protein_vector_positive.txt", "r")

list_of_lists = []
for line in a_file:
  stripped_line = line.strip()
  line_list = stripped_line.split()
  list_of_lists.append(line_list)

a_file.close()

import pandas as pd
dataframe = pd.DataFrame(list_of_lists)
dataframe['value'] =1
dataframe.head()

#making data for glycosylation
a_file = open("protein_vector_negative.txt", "r")

list_of_lists_n = []
for line in a_file:
  stripped_line = line.strip()
  line_list = stripped_line.split()
  list_of_lists_n.append(line_list)

a_file.close()

import pandas as pd
dataframe1 = pd.DataFrame(list_of_lists_n)
dataframe1['value'] =0
dataframe1.head()

frames = [dataframe,dataframe1]
result = pd.concat(frames,ignore_index= True)
result = result.sample(frac=1).reset_index(drop=True)


#SVM model train
x = result.iloc[:,:-1]
y = result.iloc[:,-1]

from sklearn.model_selection import train_test_split
X_train, X_test, y_train, y_test = train_test_split(x, y, test_size = 0.20)

from sklearn.svm import SVC
svclassifier = SVC(kernel='poly')
svclassifier.fit(X_train, y_train)

y_pred = svclassifier.predict(X_test)

from sklearn.metrics import classification_report, confusion_matrix
print(confusion_matrix(y_test,y_pred))
print(classification_report(y_test,y_pred))


