## Initial steps

This readme provides instructions towards the usage of the included jupyter notebook and command line scripts to reproduce the findings 
of our study to identify lncRNA subtypes.

The folder *feature extraction scripts* contains Python scripts used to process the datasets, i.e. for data extraction and feature engineering. They are provided as an overview of the elaborate steps involved in data preprocessing.

The folder *machine learning steps* contains Python scripts that can be used to reproduce the training steps and generate the performance metrics in a terminal, should you choose not to use the iPython environment. See [Usage](#Usage) for more information.

*Prerequisites* are python >=3.6, jupyter, and anaconda.

To make life easier, all the necessary python packages are compiled into a conda environment and provided here.

To activate the provided conda environment, use

`conda env create -f lncRNA_classify_environment.yaml`

More information regarding managing environments can be found [here](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html).

_If you would rather build everything from scratch, we would recommend not to skip the next part, otherwise please head over to_ Usage.

To install anaconda, please follow the instructions on the [website](https://docs.anaconda.com/anaconda/install/).

To install jupyter, type:

`pip install notebook`

or use conda:

`conda install -c conda-forge notebook`

For more choices, please follow the instructions [here](https://jupyter.org/install).

The python packages required to successfully execute the training are:

* [pandas](https://pandas.pydata.org/)
* [numpy](https://numpy.org)
* [scikit-learn](https://scikit-learn.org/stable/) (release >=0.22 is recommended)
* [matplotlib](https://matplotlib.org/)


## Usage

### iPython environment

To open the notebook, navigate to the notebook's folder and execute

`jupyter notebook`

and choose **lncRNA_classify.ipynb** from the tree menu.

Each _index_ of _training\_\*_ folders stands for a _dataset_.

> training_1 stands for dataset 1

Each _training\_\*_ directory contains subdirectories labelled _fs\_\*_. Each _index_ stands for a _feature set_.

> fs1 stands for feature set 1

Every _fs\_\*_ directory two subdirectories, _w\_fickett_ and _wo\_fickett_, which stand for with and without _Fickett score_, respectively. These two directories contain the files necessary for the execution of the code.

Every _fs1_ folder contain a database of the _kmer weights_ transformed using a multi-label binarizer.

The whole notebook is divided into four parts for the four datasets and each part is further divided into a _supervised module_ and an _unsupervised module_.

The code to train the classifier using _feature set 1_ with _Fickett score_ is already provided. 

Only the directories _fs\[1-3\]_ need to be changed to proceed with the training using the features containing structural information and conservation scores. Inclusion of the _Fickett score_ can be regulated choosing the appropriate folder for every feature set.

### Command line

Two python scripts are included in the _machine learning steps_ folder, namely `supervised.py` and `unsupervised.py` to enable command line usage.

`supervised.py` trains a random forest classifier on the training sets as described above.

`unsupervised.py` performs PCA and k-means clustering on the training sets.

Each script can be run with the following switches:

`-t, -training_set` choose the training set for the classifier: \[1,2,3,4\]. The default value is 1.

`-fs, -feature_set` choose the feature set: \[1,2,3\]. The default value is 1.

`-fi, -fickett_score` should you choose to include fickett score, enter 1, else 0. Disabled by default.

Help with the switches can be accessed by:

`supervised.py -h`

OR

`unsupervised.py -h`

If no parameters are specified, both the files will run on default parameters.

`unsupervised.py` will save the plots generated for PCA and k-means clustering in the working directory.
