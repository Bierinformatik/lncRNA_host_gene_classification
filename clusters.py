from os.path import abspath, join
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import OrdinalEncoder
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import classification_report, adjusted_rand_score, f1_score
from sklearn.metrics import accuracy_score, confusion_matrix, plot_confusion_matrix, auc, plot_roc_curve
import matplotlib.pyplot as plt
from collections import Counter
import numpy as np
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.metrics.pairwise import pairwise_distances_argmin

def preprocess(features_combined,args):
	features_combined.drop(['gene_id'], 1, inplace=True)

	#for phatscons/fs3
	if args.fs==3:
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
	y_test = np.ravel(y_test)

	print (f'Shape of X_train after encoding: {X_train.shape}')
	print (f'Shape of y_train after encoding: {y_train.shape}')
	print (f'Shape of X_test after encoding: {X_test.shape}')
	print (f'Shape of y_test after encoding: {y_test.shape}')
	
	print ("PCA with 2 components")
	pca (features_combined, X_train, y_train, X_test, y_test, args)
	
def pca(features_combined, X_train, y_train, X_test, y_test, args):
	p = PCA(n_components=2, random_state=42)
	X= p.fit_transform(X_train)
	print(f'explained variance ratio: {p.explained_variance_ratio_}')
	print(f'n_components: {p.n_components_}')
	print(f'features: {p.n_features_}')
	print(f'samples: {p.n_samples_}')
	fig = plt.figure(figsize=(8,3))
	colours = ['navy', 'turquoise', 'darkorange']
	targets = ['NoHG', 'MIRHG', 'SNHG']
	lw=2
	for colour, i, target_name in zip(colours, [0, 1, 2], targets):
    		plt.scatter(X[y_train == i, 0], X[y_train == i, 1], color=colour, alpha=.8, #lw=lw,
                label=target_name)
	plt.legend(loc='best', shadow=False, scatterpoints=1)
	#plt.title('PCA of dataset4/feature set 1/wo_fickett')
	plt.text(200,30, f'PC1={p.explained_variance_ratio_[0]}\nPC2={p.explained_variance_ratio_[1]}',
		fontsize='small', bbox=dict(edgecolor='red', facecolor='white', alpha=0.3,linestyle=':'))
	fig.savefig("pca_ds"+str(args.t)+"_feat_"+str(args.fs)+"_fickett"+str(args.fi)+".svg", format="svg")
	
	print ("\n\nk-means with 3 clusters")
	kmeans(features_combined, X_train, y_train, X_test, y_test, args)

def kmeans(features_combined, X_train, y_train, X_test, y_test, args):
	
	km = KMeans (n_clusters=3, init="k-means++", random_state=42).fit(X_train)

	#print (f'kmeans labels : {km.labels_}')

	#print (f'number of iter : {km.n_iter_}')

	#print (f'score : {km.score(X_test)}')

	y = km.predict(X_test)

	#print(classification_report(y_test, y, target_names=['lnc', 'mir', 'sno']))

	print (confusion_matrix(y_test, y))

	print (f'Adjusted rand score: {adjusted_rand_score(y_test, y)}')

	print(f'X_pred: {Counter(km.labels_)}\n')
	print(f'X_labels: {Counter(y_train)}\n')
	print(f'accuracy training set: {accuracy_score(y_train, km.labels_)}\n')
	print(f'accuracy test set: {accuracy_score(y_test, y)}\n')
	print(f'num of correctly classified samples: {accuracy_score(y_train, km.labels_, normalize=False)}\n')
	
	k_means_labels = pairwise_distances_argmin(X_test, km.cluster_centers_)
	fig = plt.figure(figsize=(8, 3))	
	colors = ['#4EACC5', '#FF9C34', '#4E9A06']
	targets = ['NoHG', 'MIRHG', 'SNHG']
	for k, col, target in zip(range(3), colors, targets):
		my_members = k_means_labels == k
		cluster_center = km.cluster_centers_[k]
		plt.plot(X_test[my_members, 0], X_test[my_members, 1], 'w',
            markerfacecolor=col, marker='.')
		plt.plot(cluster_center[0], cluster_center[1], 'o', markerfacecolor=col,
            markeredgecolor='k', markersize=6, label=target)	
	#plt.title('ds2_fs3_wf_kmeans_on_test')
	plt.legend(loc='best', shadow=False, scatterpoints=1)
	fig.savefig("kmeans_ds"+str(args.t)+"_feat_"+str(args.fs)+"_fickett"+str(args.fi)+".svg", format="svg")

