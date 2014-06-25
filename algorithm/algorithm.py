# Runs the required steps for the algorithm
#
# @author Clemens Westrup

import argparse
import os
import errno
import shutil
import sys
import subprocess

 # default max depth value used in the paper
maxdepth = 50


# parse command line arguments
def parseargs():
    parser = argparse.ArgumentParser(
        description='Runs the indexing .'
        + ' Results are stored in folder' + os.path.abspath(outputdir))
    parser.add_argument(
        '-m', '--maxdepth', type=str, dest='maxdepth',
        required=True, help='maximal tree depth')
    parser.add_argument(
        '-v', '--verbose', dest='verbose', action='store_true',
        help='show verbose information')
    args = parser.parse_args()
    return args


# mkdir -p
def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise

# main function
if __name__ == '__main__':

    # parse arguments
    args = parseargs()
    maxdepth = args.maxdepth


#### trie-generator

# TrieGenerator.java converts reactiongraph .mol files into .seqs
# files listing all the paths with a specified depth
java TrieGenerator maxdepth dir/with/reactiongraph-mol-files/*.mol

# Concatenate those files to a single file (e.g. cat-rgraphs.trees)
cat dir/with/seqs-files/*.seqs > cat-rgraphs.trees

### treeBorrowsWheelerTransform

# convert the encoding of the graph to a unicode encoding via
./tconvert < dir/with/cat-rgraphs.trees > output.graph 2>encoding.txt

# builder is used to create an index (output.graph.tbwt)
./builder output.graph 

# traverse then creates resulting path frequencies out of the index generated 
# by builder
./traverse output.graph.tbwt > output.freqs

# and resultconvert translates the encoding back that tconvert produced
./resultconvert cat-rgraphs.trees encoding.txt < output.freqs > result-kegg.freqs

### kernels




