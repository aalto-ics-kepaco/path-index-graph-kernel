path-index-graph-kernel
=======================

An implementation of an efficient graph kernel using a compressed path index. 

Contents
--------

* algorithm - The implementation of the algorithm to be used with .sdf or .mol files as input
* build.sh - A script that builds / compiles all the used scripts for the current system
* examples - Example use cases for the algorithm
    * keggReactionPrediction - Prediction of reaction ec numbers using the [KEGG dataset](http://www.genome.jp/kegg/)
* README.md - This readme

Usage
-----

** Prerequisites **

* python
* java 
* a c++ compiler
* matlab for some of the examples

** Building **

First you need to build / compile some of the files to be used on your system. The script `build.sh` does all this for you, just run it once and everything should be in place.

** Running the algorithm **

The algorithm currently only supports .sdf (e.g. PubChem) and .mol files (e.g. KEGG) as input graphs. Running the wrapper script (`python runAlgorithm.py`) in the `algorithm` directory  

Examples
--------


Old stuff below this line
________________________________



#### keggPreProcessing/1-reaction-list

Python script to create a condensed list of the reactions.

* Create KEGG reaction listing: `python extract-reactions.py path/to/kegg/latest/ path/to/kegg/latest/mol/ path/to/kegg/latest/ > kegg-reactions.txt`. The 3 arguments are path to "reaction" file, path to kegg mol's and path to kegg.

#### keggPreProcessing/2-feature-generator

Python code to create atom features from KEGG mol files

* Create atom features using `python generator2011.py path/to/kegg/latest/mol/* -k all --output-dir path/to/results/mol-features/`. The optional parameter k specifies the context size, default is all. There are other optional parameters for this script.

#### keggPreProcessing/3-atommapper

Java code to create reaction graphs out of mol files, feature files and the reaction-list

* Create reactiongraph mol files with atommapper: `java Mapper2000 -moldir  path/to/kegg/latest/mol/ -featdir path/to/atom/features/ -reacfile path/to/reactionlist/kegg-reactions.txt -output ../../../../results/computed/3-atommapper/`


Java code to extract all paths from each node as tries.

* Concatenate those files to a single file (e.g. cat-rgraphs.trees): `cat dir/with/seqs-files/*.seqs > cat-rgraphs.trees` 


C++ code to traverse all these trees and count the path frequencies of each graph


* create a listing of all .mol files in the reaction graph directory (results of keggPreProcessing/3-atommapper): `cd resultdir/of/3-atommapper` and `ls -1 *.mol > rgraphlist.txt`
* convert the .freqs result to a matlab sparse matrix file with `python freq2mtl.py path/to/result-kegg.freqs path/to/rgraphlist.txt`
* uncomment the kernels you want to use in `genkernelmatrices.m`, adjust the paths and run the script
*