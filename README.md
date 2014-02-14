path-index-graph-kernel
=======================

An implementation of an efficient graph kernel using a compressed path index.

Currently this is only working for the [KEGG dataset](http://www.genome.jp/kegg/) but the aim is to generalize this to any graph dataset. 

#### About the repository contents

* keggPreProcessing - The preprocessing steps to use the method with the KEGG dataset
* treeBorrowsWheelerTransform - The tree borrows-wheeler transform to convert path listings that were extracted from graphs into path frequencies
* kernel - a python script to convert the TBWT results into .mtl files and a matlab script to generate various kernels

Usage for KEGG
--------------

First make sure you have python, java and a c++ compiler up and running on your system. Then build keggPreProcessing/3-atommapper and keggPreProcessing/4-trie-generator for Java and treeBorrowsWheelerTransform with make (you'll need the builder and traverse).

#### keggPreProcessing/1-reaction-list

Python script to create a condensed list of the reactions.

* Create KEGG reaction listing: `python extract-reactions.py path/to/kegg/latest/ path/to/kegg/latest/mol/ path/to/kegg/latest/ > kegg-reactions.txt`. The 3 arguments are path to "reaction" file, path to kegg mol's and path to kegg.

#### keggPreProcessing/2-feature-generator

Python code to create atom features from KEGG mol files

* Create atom features using `python generator2011.py path/to/kegg/latest/mol/* -k all --output-dir path/to/results/mol-features/`. The optional parameter k specifies the context size, default is all. There are other optional parameters for this script.

#### keggPreProcessing/3-atommapper

Java code to create reaction graphs out of mol files, feature files and the reaction-list

* Create reactiongraph mol files with atommapper: `java Mapper2000 -moldir  path/to/kegg/latest/mol/ -featdir path/to/atom/features/ -reacfile path/to/reactionlist/kegg-reactions.txt -output ../../../../results/computed/3-atommapper/`
#### keggPreProcessing/4-trie-generator

Java code to extract all paths from each node as tries.
* `TrieGenerator.java` converts reactiongraph .mol files into .seqs files listing all the paths with a specified depth: `java TrieGenerator <maxdepth> dir/with/reactiongraph-mol-files/*.mol`
* Concatenate those files to a single file (e.g. cat-rgraphs.trees): `cat dir/with/seqs-files/*.seqs > cat-rgraphs.trees` 
#### treeBorrowsWheelerTransform

C++ code to traverse all these trees and count the path frequencies of each graph* convert the encoding of the graph to a unicode encoding via `./tconvert < dir/with/cat-rgraphs.trees > output.graph 2>encoding.txt`. This enables the use of edge labels and multi-character labels that are encoded with a single character. The encoding key is provided via stderr.* `builder` is used to create an index (output.graph.tbwt): `./builder output.graph`* `traverse` then creates resulting path frequencies out of the index generated by builder: `./traverse output.graph.tbwt > output.freqs`* and `resultconvert` translates the encoding back that tconvert produced: `./resultconvert cat-rgraphs.trees encoding.txt < output.freqs > result-kegg.freqs`
#### kernel

* create a listing of all .mol files in the reaction graph directory (results of keggPreProcessing/3-atommapper): `cd resultdir/of/3-atommapper` and `ls -1 *.mol > rgraphlist.txt`
* convert the .freqs result to a matlab sparse matrix file with `python freq2mtl.py path/to/result-kegg.freqs path/to/rgraphlist.txt`
* uncomment the kernels you want to use in `genkernelmatrices.m`, adjust the paths and run the script
*