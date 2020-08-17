## Initial steps

This readme provides instructions towards the usage of the included jupyter notebook to reproduce the findings 
of our study to identify lncRNA subtypes.

Prerequisites are python >=3.6, jupyter, and anaconda.

To make life easier, all the necessary python packages are compiled into a conda environment and provided here.

To activate the provided conda environment, use

`conda env create -f lncRNA_classify_environment.yaml`

More information regarding managing environments can be found [here](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html).

_If you would rather build everything from scratch, we would recommend not to skip the next part._

To install anaconda, please follow the instructions on the [website](https://docs.anaconda.com/anaconda/install/).

To install jupyter, type:

`pip install notebook`

or use conda:

`conda install -c conda-forge notebook`

For more choices, please follow the instructions [here](https://jupyter.org/install) 

The python packages required to successfully execute the training are:

* [pandas](https://pandas.pydata.org/)
* [numpy](https://numpy.org)
* [scikit-learn](https://scikit-learn.org/stable/) (release >=0.22 is recommended)
* [matplotlib](https://matplotlib.org/)


## Usage

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
