#!/bin/sh

# builds all needed files

curr_dir=`pwd`

dir=`dirname $0`
FILE_PATH=`cd  $dir;pwd`

echo "building atommapper"
ATOMMAPPERPATH=$FILE_PATH/examples/keggReactionPrediction/atommapper
javac -d $ATOMMAPPERPATH/bin -cp $ATOMMAPPERPATH/src/mapper/ \
-sourcepath $ATOMMAPPERPATH/src $ATOMMAPPERPATH/src/Mapper2000.java

echo "building trie-generator"
TRIEGENERATORPATH=$FILE_PATH/algorithm/trie-generator
javac -d $TRIEGENERATORPATH/bin -cp $TRIEGENERATORPATH/src/mechanism/ \
-sourcepath $TRIEGENERATORPATH/src $TRIEGENERATORPATH/src/TrieGenerator.java

echo "building treeBorrowsWheelerTransform"
TBWTPATH=$FILE_PATH/algorithm/treeBorrowsWheelerTransform
cd $TBWTPATH
make clean
make all
make tconvert
make resultconvert