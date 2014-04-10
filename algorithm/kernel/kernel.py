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
        '(.mol or .sdf files). Watch out there\'s no other files in that '
        'folder, hidden files are ignored.')
    parser.add_argument('-t', '--tbwtresult', type=str, dest='tbwtresult', 
        required=True, help='path to tbwt result file (e.g. result.freqs)')
    parser.add_argument('-o', '--outputpath', type=str, dest='output', 
        required=True, help='path to store output files '
        '(\"features\" and \"kernels\"')
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

    # remove ending slash on path if existing
    if graphpath.endswith('/'):
        graphpath = graphpath[:-1]

    # read all graphs and sort them in alphanumerical order
    graph_list = sorted(os.listdir(graphpath))

    # read tbwt results file and convert it to a list of lines
    with open(args.tbwtresult, 'r') as tbwtresultfile:
        tbwt_list = tbwtresultfile.read().splitlines()

    # iterate over all graphs listed in the directory
    for i in range(len(graph_list)):

        # remove file name extension
        graph_list[i] = os.path.splitext(graph_list[i])[0]
        graph = graph_list[i]

        # ignore hidden files / files with names starting with .
        if graph[0] == ".":
            continue

        # for each reaction check the frequencies of each path
        for line in tbwt_list:
            # if graph name is found in this line extract frequency
            matches = re.search(graph + r':\d', line)
            if matches:
                matches_list = matches.group(0)
                # get frequency by removing graph name + one char for ":"
                frequency = int(matches_list[len(graph)+1:])
                logger.debug(graph + ": " + str(frequency))
            else: 
                logger.debug(graph + ": no match")



    logger.debug(graph_list)

