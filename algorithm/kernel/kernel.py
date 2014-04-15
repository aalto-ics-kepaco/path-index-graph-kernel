# Extracts feature map Phi and the kernel matrix K from tbwt results 
# that can be used as input for MMCRF
#
# @author Clemens Westrup

import argparse, logging, os, re, sys
from scipy.sparse import *
from scipy.spatial.distance import *
from scipy import *
import numpy as np
from itertools import izip

# parse command line arguments
def parseargs():
    parser = argparse.ArgumentParser(description='Extracts feature map Phi and '
        'if specified by parameter writefeatures and computes the kernel '
        'matrix K from tbwt results. Note: Matrices are sorted '
        'in alphanumerical order of graph filenames so labels should be too.')
    parser.add_argument('-g', '--graphpath', type=str, dest='graphpath', 
        required=True, help='path to folder with input graphs for tbwt '
        '(.mol or .sdf files, hidden files are ignored) OR path to file '
        ' containing a listing of graphs in the format of '
        '\"ls -1 {graphdir}/*.mol > {outputdir}/graphlist.txt\"')
    parser.add_argument('-t', '--tbwtresult', type=str, dest='tbwtresult', 
        required=True, help='path to tbwt result file (e.g. result.freqs)')
    parser.add_argument('-o', '--outpath', type=str, dest='outpath', 
        required=True, help='path to store output files '
        '(\"kernels\" and if selected \"features\"')
    parser.add_argument('-p', '--prefix', type=str, dest='prefix', 
        help='add a prefix to the output files')    
    parser.add_argument('-c', '--common', type=int, dest='common', 
        help='only take paths into consideration that are shared between more '
        'than {common} graphs')
    parser.add_argument('-wf', '--writefeatures', dest='writefeatures', 
        action='store_true', help='write out features in a file')
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true',
        help='show verbose information')
    parser.add_argument('-f', '--force', dest='force', action='store_true',
        help='force overwrite if output dir exists already')

    args = parser.parse_args()
    return args

# check if file exists and if force is set to overwrite it
def check_output(outputfile, force):
    if (os.path.isfile(outputfile)):
        if (force):
            try:
                os.remove(outputfile)
                logger.warn('Force parameter is set: Output file ' + outputfile 
                    + ' exists and will be overwritten')

            except OSError:
                logger.error("Output file " + outputfile + "exists already "
                + "but can't be removed. Aborting.");
                sys.exit(1)
        else:
            logger.error('Output file ' + outputfile + 
                ' exists already, aborting.');
            sys.exit(1)

def get_output_line_from_vector(vector, cast_integer):
    output_line = ""
    for entry in vector:
        if cast_integer:
            entry = int(entry)
        output_line += str(entry) + " "
    output_line = output_line[:-1] + "\n"
    return output_line

def compute_feature_vector(graph, tbwt_list):
    """Compute a feature vector for a graph from its path frequencies.

    Keyword arguments:
    graph -- name of the graph 
    tbwt_list -- list of tbwt path entries read from the result file
    """
    # prepare empty array 
    feature_vector = zeros((len(tbwt_list),1), dtype=float)
    # for each reaction get the frequencies of each path
    for index,item in enumerate(tbwt_list):
        # if graph name is found in this line extract frequency
        matches = re.search(graph + r':\d', item)
        if matches:
            matches_list = matches.group(0)
            # get frequency by removing graph name + one char for ":"
            frequency = int(matches_list[len(graph)+1:])
            feature_vector[index] = frequency
    return feature_vector

def compute_kernel(phi_1, phi_2, kernel_type):
    if kernel_type == "linear":
        return phi_1.T.dot(phi_2)

# main function
if __name__ == '__main__':

    # parse input arguments 
    args = parseargs()
    outpath = args.outpath
    graphpath = args.graphpath

    # get logger and set it up
    logging.basicConfig(level=logging.INFO, format="%(levelname)s - %(message)s")
    logger = logging.getLogger()
    # show debug level log messages if verbose is set
    if (args.verbose):
        logger.setLevel(logging.DEBUG)

    # add folder extension for output files
    if not outpath.endswith('/'):
        outpath = outpath + '/'
    # add prefix if set
    if args.prefix:
        outpath  = outpath + args.prefix

    # check if graphpath is a file listing or a dir with the graphs
    if os.path.isfile(graphpath):
        logger.debug("reading graphs from file " + graphpath)
        with open(graphpath, 'r') as graph_listing:
            graph_list = graph_listing.read().splitlines()
    else:
        logger.debug("reading graphs from directory " + graphpath)
        # read all graphs and sort them in alphanumerical order
        graph_list = sorted(os.listdir(graphpath))

    # read tbwt results file and convert it to a list of lines
    with open(args.tbwtresult, 'r') as tbwtresultfile:
        tbwt_list = tbwtresultfile.read().splitlines()

    # if common parameter is set remove all paths with only one graph listed
    if (args.common):
        new_tbwt_list = []
        for item in tbwt_list:
            if (not re.match(r"^\S*(?:\s\S*){0," + str(args.common) + "}$", item)):
                new_tbwt_list.append(item)

        logger.debug(str(len(tbwt_list)-len(new_tbwt_list)) +
            " of " + str(len(tbwt_list)) + " paths removed " +
            "that are shared between <= " + str(args.common) + " graphs.")
        tbwt_list = new_tbwt_list

    # remove all file extensions from graph names 
    # and delete hidden file entries
    graph_list_copy = [];
    for idx, graph in enumerate(graph_list):
        # not hidden file
        if graph[0] != ".": 
            # remove file extension
            graph_list_copy.append(os.path.splitext(graph_list[idx])[0])
    graph_list = graph_list_copy

    # open feature output file to write feature map if requested
    if args.writefeatures:
        featurefile = open(outpath + 'features', 'w')

    # open kernel output file
    kernelfile = open(outpath + 'kernels', 'w')

    # generate diagonal kernels for normalization
    diagonal_kernels = np.zeros(len(graph_list))
    for i, graph_i in enumerate(graph_list): 
        phi_i = compute_feature_vector(graph_i, tbwt_list)
        diagonal_kernels[i] = compute_kernel(phi_i, phi_i, "linear")
        logger.debug("computing diagonal kernel #" + str(i))

    # iterate over all graphs to compute all kernels from feature vectors with 
    # the feature vector of the current graph (to avoid storage of full
    # feature matrix)
    for i, graph_i in enumerate(graph_list):

        # output graph 
        logger.debug("computing kernel row of graph i = " + str(i))

        # prepare array as line of kernel matrix
        kernel_matrix_row_i = np.zeros(len(graph_list))

        # feature vector phi for graph i
        phi_i = compute_feature_vector(graph_i, tbwt_list)

        # write out feature vector
        if args.writefeatures:
            output_line = get_output_line_from_vector(phi_i, True)
            featurefile.write(output_line)

        # iterate over all graphs again to generate the kernels with graph i
        for j, graph_j in enumerate(graph_list):

            # output graph 
            logger.debug("computing kernel value for graphs " 
                + str(i) + " and " + str(j))

            # feature vector phi for graph i
            phi_j = compute_feature_vector(graph_j, tbwt_list)

            # kernel for graphs i and j
            kernel_matrix_row_i[j] = compute_kernel(phi_i, phi_j, "linear")

        # normalize kernel row
        kernel_matrix_row_i_normalized = zeros(len(graph_list))
        for j, kernel in enumerate(kernel_matrix_row_i):
            kernel_matrix_row_i_normalized[j] = (kernel_matrix_row_i[j] 
                / sqrt(diagonal_kernels[i] * diagonal_kernels[i]))

        # write out normalized kernel row
        output_line = get_output_line_from_vector(kernel_matrix_row_i_normalized, False)
        kernelfile.write(output_line)

    # close output files 
    kernelfile.close()
    if args.writefeatures:
        featurefile.close()

