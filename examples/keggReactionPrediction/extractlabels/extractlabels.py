# Extracts the multilabel matrix Y for the MMCRF prediction task from KEGG
#
# @author Clemens Westrup

import argparse
import os
import re
import logging
import sys


# parse command line arguments
def parseargs():
    parser = argparse.ArgumentParser(
        description='Extracts the multilabel matrix Y for the MMCRF '
                    'prediction task from KEGG')
    parser.add_argument(
        '-k', '--keggpath', type=str, dest='keggpath',
        required=True, help='path to kegg ligand database')
    parser.add_argument(
        '-r', '--reactiongraphsdir', type=str, dest='rgrapgsdir',
        required=True, help='path to directory containing used kegg '
        'graphs as input for tbwt')
    parser.add_argument(
        '-o', '--output', type=str, dest='output',
        required=True, help='path to directory to store output result')
    parser.add_argument(
        '-i', '--ignoreincomplete', dest='ignoreincomplete',
        action='store_true', help='Ignore incomplete ec-labels like 1.3.2.-')
    parser.add_argument(
        '-d', '--depth', dest='depth', type=int, default=4,
        help='Depth of EC hierarchy level, 1 for only keeping the top level'
             + 'e.g. EC 1 up tp 4 for keeping the whole tree e.g. EC 1.4.2.23')
    parser.add_argument(
        '-v', '--verbose', dest='verbose', action='store_true',
        help='show verbose information')
    parser.add_argument(
        '-f', '--force', dest='force', action='store_true',
        help='force overwrite if output dir exists already')

    args = parser.parse_args()
    return args


# check if file exists and if force is set to overwrite it
def check_output(outputfile, force):
    if (os.path.isfile(outputfile)):
        if (force):
            try:
                os.remove(outputfile)
                logger.warn("Output file " + outputfile
                            + " exists and will be overwritten"
                            + " since parameter force was set.")
            except OSError:
                logger.error("Output file " + outputfile + "exists already"
                             + "but can't remove. Aborting.")
                sys.exit(1)
        else:
            logger.error('Output file ' + outputfile
                         + ' exists already, aborting.')
            sys.exit(1)

# main function
if __name__ == '__main__':

    # get logger and set it up
    logging.basicConfig(level=logging.INFO,
                        format="%(levelname)s - %(message)s")
    logger = logging.getLogger()

    # parse input arguments
    args = parseargs()

    # convert to absolute paths
    keggpath = os.path.abspath(args.keggpath)
    rgrapgsdir = os.path.abspath(args.rgrapgsdir)

    # remove ending slash on path if existing
    if keggpath.endswith('/'):
        keggpath = keggpath[:-1]

    # show debug level log messages if verbose is set
    if (args.verbose):
        logger.setLevel(logging.DEBUG)

    # prepare output
    check_output(os.path.abspath(args.output + "/labels"), args.force)
    check_output(os.path.abspath(args.output + "/graphlist.txt"), args.force)

    # dictionary to store output and set for unique ec-codes
    reactiondir = {}
    eccodes = set()

    # read all reactiongraphs
    reaction_list = sorted(os.listdir(rgrapgsdir))

    logger.debug("reactionlist: ")
    logger.debug(reaction_list)

    # ec-code search pattern
    if (args.ignoreincomplete):
        ecpattern = r"((?:(?:\d+)\.){3}(?:\d+))"
    else:
        ecpattern = r"((?:(?:\d+|-)\.){3}(?:\d+|-))"

    # process kegg reactions file
    with open(keggpath + '/reaction', 'r') as infile:
        reactions_found = 0
        reactions_total = 0
        reaction_read = False
        ecnumbers_read = False
        reading_ecnumbers = False

        for line in infile:

            # new reaction entry
            if re.match('ENTRY.*Reaction', line):
                reactions_total += 1
                try:
                    reaction = re.search('R[0-9]{5}', line).group(0)
                    reaction_read = True
                except AttributeError:
                    pass

            # check if reading ec-number lines
            elif re.match('ENZYME', line):
                reading_ecnumbers = True
                ecnumbers = []
            elif re.match('(ORTHOLOGY)|(///)|(REFERENCE)', line):
                reading_ecnumbers = False
                ecnumbers_read = True
            # extract ec numbers if currently reading those
            if reading_ecnumbers:
                try:
                    new_ecnumbers = re.findall(ecpattern, line)
                    for ec in new_ecnumbers:
                        ecnumbers.append(".".join(ec.split(".")[0:args.depth]))
                except AttributeError:
                    pass

            # if an entry was read (meaning reaction number and ec codes)
            # then store the reaction entry with eccodes
            if reaction_read and ecnumbers_read:
                reaction_read = False
                ecnumbers_read = False

                # iterate over reactions in graphs directory
                for listed_reaction in reaction_list:
                    # remove file extension
                    listed_reaction_basename = os.path.splitext(
                        listed_reaction)[0]
                    # if the reaction exists then store it with the filename
                    # e.g. reaction could be R00258 and a matching
                    # listed_reaction would be R00258_0_b.mol
                    # where R00258_0_b is stored
                    if (re.match(reaction, listed_reaction_basename)):
                        reactions_found += 1
                        reactiondir[listed_reaction_basename] = ecnumbers
                        # add ec numbers to set of unique ec codes
                        for ec in ecnumbers:
                            eccodes.add(ec)

     # convert unique ec-codes to sorted list
    sorted_ecs = sorted(list(eccodes))

    logger.debug(sorted_ecs)

    logger.info("found " + str(reactions_found)
                + " reactions out of " + str(reactions_total))
    logger.info("extracted " + str(len(sorted_ecs))
                + " unique ec-codes from found reactions")

    # write output file
    with open(args.output + "/labels", 'w') as outfile:
        with open(args.output + "/graphlist.txt", 'w') as graphlist_output:
            for reaction in sorted(reactiondir):
                logger.debug("found reaction " + reaction + " with ec-codes:")

                # prepare line of output for a new graph with its labels
                output_line = ""

                # write line to reaction listing
                graphlist_output.write(reaction + ".mol\n")

                # write line with binary vector for graph
                # with eccodes as labels
                for ec in sorted_ecs:
                    if ec in reactiondir[reaction]:
                        output_line += "1 "
                        logger.debug(ec)
                    else:
                        output_line += "0 "
                output_line = output_line[:-1] + "\n"
                outfile.write(output_line)
