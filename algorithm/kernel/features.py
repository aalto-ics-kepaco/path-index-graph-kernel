# Extracts features / path counts for graphs as .mtl sparse matrix files for
# use in matlab
#
# @author Clemens Westrup

import argparse
import logging
import os
import re
import sys
from scipy.sparse import *
from scipy.spatial.distance import *
from scipy import *


def parseargs():
    """Parse command line arguments.
    """
    parser = argparse.ArgumentParser(
        description='Extracts features (= path counts) for graphs as .mtl '
        'sparse matrix files for use in matlab. Note: Matrices are sorted '
        'in alphanumerical order of graph filenames so labels should be too.')
    parser.add_argument(
        '-g', '--graphpath', type=str, dest='graphpath', required=True,
        help='path to folder with input graphs for tbwt '
        '(.mol or .sdf files, hidden files are ignored) OR path to file '
        ' containing a listing of graphs in the format of '
        '\"ls -1 {graphdir}/*.mol > {outputdir}/graphlist.txt\"')
    parser.add_argument(
        '-t', '--tbwtresult', type=str, dest='tbwtresult',
        required=True, help='path to tbwt result file (e.g. result.freqs)')
    parser.add_argument(
        '-o', '--outpath', type=str, dest='outpath',
        required=True, help='path to store output files '
        '(\"kernels\" and if selected \"features\"')
    parser.add_argument(
        '-p', '--prefix', type=str, dest='prefix',
        help='add a prefix to the output files')
    parser.add_argument(
        '-c', '--common', type=int, dest='common',
        help='only take paths into consideration that are shared between more '
        'than {common} graphs')
    parser.add_argument(
        '-wf', '--writefeatures', dest='writefeatures',
        action='store_true', help='write out features in a file')
    parser.add_argument(
        '-v', '--verbose', dest='verbose', action='store_true',
        help='show verbose information')
    parser.add_argument(
        '-f', '--force', dest='force', action='store_true',
        help='force overwrite if output dir exists already')

    args = parser.parse_args()
    return args


def check_output(outputfile, force):
    """Check if file exists and if force is set to delete it.

    Keyword arguments:
    outputfile -- File to check
    force -- force overwrite / delete beforehand
    """
    if (os.path.isfile(outputfile)):
        if (force):
            try:
                os.remove(outputfile)
                logger.warn('Force parameter is set: Output file ' + outputfile
                            + ' exists and will be overwritten')

            except OSError:
                logger.error("Output file " + outputfile + "exists already "
                             + "but can't be removed. Aborting.")
                sys.exit(1)
        else:
            logger.error('Output file ' + outputfile
                         + ' exists already, aborting.')
            sys.exit(1)


def get_output_line_from_vector(vector, cast_integer):
    """Generate a string from a vector to be use a line to be written
    to an output file

    Keyword arguments:
    vector -- Vector to generate output from
    cast_integer -- whether values should be casted to integers
    """
    output_line = ""
    for entry in vector:
        if cast_integer:
            entry = int(entry)
        output_line += str(entry) + " "
    output_line = output_line[:-1] + "\n"
    return output_line


# main function
if __name__ == '__main__':

    # parse input arguments
    args = parseargs()
    outpath = args.outpath
    graphpath = args.graphpath

    # get logger and set it up
    logging.basicConfig(level=logging.INFO,
                        format="%(levelname)s - %(message)s")
    logger = logging.getLogger()
    # show debug level log messages if verbose is set
    if (args.verbose):
        logger.setLevel(logging.DEBUG)

    # add folder extension for output files
    if not outpath.endswith('/'):
        outpath = outpath + '/'
    # add prefix if set
    if args.prefix:
        outpath = outpath + args.prefix

    # check if graphpath is a file listing or a dir with the graphs
    if os.path.isfile(graphpath):
        logger.info("Reading graphs from file " + graphpath)
        with open(graphpath, 'r') as graph_listing:
            graph_list = graph_listing.read().splitlines()
    else:
        logger.info("Reading graphs from directory " + graphpath)
        # read all graphs and sort them in alphanumerical order
        graph_list = sorted(os.listdir(graphpath))

    # remove file extension for all graphs
    graph_list_no_ext = []
    for graph in graph_list:
        graph_list_no_ext.append(os.path.splitext(graph)[0])
    graph_list = graph_list_no_ext
    num_graphs = len(graph_list)

    # read tbwt results file and convert it to a list of lines
    with open(args.tbwtresult, 'r') as tbwtresultfile:
        tbwt_list = tbwtresultfile.read().splitlines()

    num_tbwt_entries = len(tbwt_list)

    featurefile = open(outpath + 'features', 'w')

    logger.info("Extracting features from tbwt result.")
    # pattern for filtering out paths with <= args.common listed graphs
    matching_pattern = r"^\S*(?:\s\S*){0," + str(args.common) + "}$"
    removed_paths = 0
    # extract features
    for path_i, line in enumerate(tbwt_list):
        # if common parameter and current path matches pattern
        # meaning it has <= args.common graphs listed it's discarded
        if (args.common and re.match(matching_pattern, line)):
            removed_paths += 1
        else:
            # extract graphs
            entries = line.split(' ')[1:]
            # add path to each of those graph entries in feature dictionary
            for entry in entries:
                [graph, frequency] = entry.split(':')
                try:
                    graph_id = graph_list.index(graph)
                    # add 1 to path and graph to avoid 0 index for matlab
                    featurefile.write(str(graph_id + 1) + " " + str(path_i + 1)
                                      + " " + frequency + "\n")
                except ValueError:
                    pass
            logger.debug("Extracted path " + str(path_i)
                         + " out of " + str(num_tbwt_entries))
    featurefile.close()
