from os.path import abspath, join
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import OrdinalEncoder
from sklearn.compose import ColumnTransformer
from sklearn.ensemble import RandomForestClassifier
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import classification_report
from sklearn.model_selection import cross_val_score, cross_val_predict, GridSearchCV
from sklearn.metrics import accuracy_score, confusion_matrix, plot_confusion_matrix, auc, plot_roc_curve
from sklearn.preprocessing import MultiLabelBinarizer
import matplotlib.pyplot as plt
from collections import Counter
import numpy as np

def preprocess(features_combined,fs):
	features_combined.drop(['gene_id'], 1, inplace=True)

	#for phatscons/fs3
	if fs==3:
		val = {'cons_1':0, 'cons_2':0, 'cons_3':0, 'cons_4':0, 'cons_5':0, 'cons_6':0, 
		'cons_7':0, 'cons_8':0, 'cons_9':0, 'cons_10':0, 'cons_11':0, 'cons_12':0, 'cons_13':0, 
		'cons_14':0, 'cons_15':0, 'cons_16':0, 'cons_17':0, 'cons_18':0, 'cons_19':0, 'cons_20':0}
		features_combined.fillna(value=val, inplace=True)

	ncrna_train, ncrna_test = train_test_split(features_combined, test_size=0.2, random_state=42)
	ncrna_test_target=ncrna_test[["target"]]
	ncrna_test=ncrna_test.drop(["target"],1)
	print (f'Shape of y_test before encoding: {ncrna_test_target.shape}')
	print (f'Shape of X_test before encoding: {ncrna_test.shape}')

	ncrna_train_target=ncrna_train[["target"]]
	ncrna_train=ncrna_train.drop(["target"],1)
	print (f'Shape of y_train before encoding: {ncrna_train_target.shape}')
	print (f'Shape of X_train before encoding: {ncrna_train.shape}')

	#separating numerical features and categorical features
	num_feat=[]
	for x in ncrna_train.columns[1:]:
	    num_feat.append(x)
	cat_feat = ["seq"]

#encoding the sequence feature
#scaling the rest of the features
#and creating a ColumnTransformer object
	cct = ColumnTransformer([
	    ('oe', OrdinalEncoder(), cat_feat),
	    ('num', StandardScaler(), num_feat)
	    ], remainder='passthrough')

	#fitting the split training set and test set
	X_train = cct.fit_transform(ncrna_train)
	X_test = cct.fit_transform(ncrna_test)

	#encoding the target labels
	ordinal = OrdinalEncoder()
	y_train = ordinal.fit_transform(ncrna_train_target)
	y_test = ordinal.fit_transform(ncrna_test_target)
	y_train = np.ravel(y_train)

	print (f'Shape of X_train after encoding: {X_train.shape}')
	print (f'Shape of y_train after encoding: {y_train.shape}')
	print (f'Shape of X_test after encoding: {X_test.shape}')
	print (f'Shape of y_test after encoding: {y_test.shape}')
	
	print ("Training a random forest classifier with default number of trees (100)")
	rf_default (X_train, y_train, X_test, y_test)
	
def rf_default(X_train, y_train, X_test, y_test):
	#training the random forest classifier with default params
	rnd_clf = RandomForestClassifier(n_estimators=100, n_jobs=-1, random_state=42)
	rnd_clf.fit(X_train, y_train)
	y_pred = rnd_clf.predict(X_test)
	print (f'Accuracy on default params: {accuracy_score(y_test, y_pred)}')
	
	print ("Performance metrics on the original model")
	print(classification_report(y_test, y_pred, target_names=["NoHG","MIRHG","SNHG"]))
	disp = plot_confusion_matrix(rnd_clf, X_test, y_test,
                                     display_labels=["NoHG","MIRHG","SNHG"],
                                     cmap=plt.cm.Blues)
	disp.confusion_matrix
	
	grid_search_cv10(X_train, y_train, X_test, y_test)	
	
def grid_search_cv10(X_train, y_train, X_test, y_test):
	#grid search with 10-fold cv
	print ("10-fold cross validation with grid search. Number of trees [100,300,500,1000], bootstrap enabled, gini impurity or entropy, out-of-bagging")
	param_grid = {'oob_score': [True, False], 'bootstrap':[True],
              'criterion': ["gini", "entropy"], 'n_estimators':[100,300,500,1000]}
	grid_search_def10 = GridSearchCV(RandomForestClassifier(random_state=42),
	    param_grid, scoring="accuracy", n_jobs=-1, verbose=1, cv=10)
	grid_search_def10.fit(X_train, y_train)

	print (f"Grid search 10-fold cv best score: {grid_search_def10.best_score_}")
	best_params = grid_search_def10.best_params_
	print (f"Grid search accuracy best params: {best_params}")
	
	print ("Performance metrics on the best cross validated model")
	#y_pred_cv = grid_search_def10.best_estimator_.predict(X_test)
	y_pred_cv = grid_search_def10.predict(X_test)
	print (f'Accuracy: {accuracy_score(y_test, y_pred_cv)}')
	print(classification_report(y_test, y_pred_cv, target_names=["NoHG","MIRHG","SNHG"]))
	disp = plot_confusion_matrix(grid_search_def10.best_estimator_, X_test, y_test,
                                     display_labels=["NoHG","MIRHG","SNHG"],
                                     cmap=plt.cm.Blues)
	disp.confusion_matrix
