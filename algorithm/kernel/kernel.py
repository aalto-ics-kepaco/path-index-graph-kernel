# Extracts feature map Phi and the kernel matrix K from tbwt results 
# that can be used as input for MMCRF
#
# @author Clemens Westrup

import argparse, logging, os, re

# parse command line arguments
def parseargs():
    parser = argparse.ArgumentParser(description='Extracts feature map Phi and '
        'the kernel matrix K from tbwt results. Note: Matrices are sorted '
        'in alphanumerical order of graph filenames so labels should be too.')
    parser.add_argument('-g', '--graphpath', type=str, dest='graphpath', 
        required=True, help='path to folder with input graphs for tbwt '
        '(.mol or .sdf files, hidden files are ignored) OR path to file '
        ' containing a listing of graphs in the format of '
        '\"ls -1 {graphdir}/*.mol > {outputdir}/graphlist.txt\"')
    parser.add_argument('-t', '--tbwtresult', type=str, dest='tbwtresult', 
        required=True, help='path to tbwt result file (e.g. result.freqs)')
    parser.add_argument('-o', '--outputpath', type=str, dest='output', 
        required=True, help='path to store output files '
        '(\"features\" and \"kernels\"')
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true',
        help='show verbose information')
    parser.add_argument('-f', '--force', dest='force', action='store_true',
        help='force overwrite if output dir exists already')
    parser.add_argument('-c', '--common', type=int, dest='common', 
        help='only take paths into consideration that are shared between more '
        'than {mincommon} graphs')

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

# main function
if __name__ == '__main__':

    # get logger and set it up
    logging.basicConfig(level=logging.INFO, format="%(levelname)s - %(message)s")
    logger = logging.getLogger()

    # parse input arguments 
    args = parseargs()

    # show debug level log messages if verbose is set
    if (args.verbose):
        logger.setLevel(logging.DEBUG)

    # convert to absolute paths
    graphpath = os.path.abspath(args.graphpath)

    # check if graphpath is a file listing or a dir with the graphs
    if os.path.isfile(graphpath):
        logger.debug("reading graphs from file " + graphpath)
        with open(graphpath, 'r') as graph_listing:
            graph_list = graph_listing.read().splitlines()
    else:
        logger.debug("reading graphs from directory " + graphpath)
            # remove ending slash on path if existing
        if graphpath.endswith('/'):
            graphpath = graphpath[:-1]
        # read all graphs and sort them in alphanumerical order
        graph_list = sorted(os.listdir(graphpath))

    # read tbwt results file and convert it to a list of lines
    with open(args.tbwtresult, 'r') as tbwtresultfile:
        tbwt_list = tbwtresultfile.read().splitlines()

    # if common parameter is set remove all paths with only one graph listed
    if (args.mincommon):
        new_tbwt_list = []
        for item in tbwt_list:
            if (not re.match(r"^\S*(?:\s\S*){0," + str(args.mincommon) + "}$", item)):
                new_tbwt_list.append(item)

        logger.debug(str(len(tbwt_list)-len(new_tbwt_list)) + " paths removed "
            "that are shared between <= " + str(args.mincommon) + " graphs.")
        tbwt_list = new_tbwt_list

    # open output feature file to write feature map
    with open(args.output + "/pathkernel_features", 'w') as featurefile:

        # iterate over all graphs listed in the directory
        for i in range(len(graph_list)):

            # new line in features output file for the new graph
            output_line = ""

            # remove file name extension
            graph_list[i] = os.path.splitext(graph_list[i])[0]
            graph = graph_list[i]

            # ignore hidden files / files with names starting with .
            if graph[0] == ".":
                continue

            # for each reaction get the frequencies of each path
            for line in tbwt_list:
                # if graph name is found in this line  extract frequency
                matches = re.search(graph + r':\d', line)
                if matches:
                    matches_list = matches.group(0)
                    # get frequency by removing graph name + one char for ":"
                    frequency = int(matches_list[len(graph)+1:])
                    output_line += str(frequency) + " "
                else: 
                    output_line += "0 "

            # write line with path frequencies / features for graph i
            output_line = output_line[:-1] + "\n"
            featurefile.write(output_line)

