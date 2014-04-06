# Extracts the multilabel matrix Y for the MMCRF prediction task
#
# @author Clemens Westrup

import argparse, os, time, re, logging, sys

# parse command line arguments
def parseargs():
    parser = argparse.ArgumentParser(description='Extracts the multilabel'
        'matrix Y for the MMCRF prediction task')
    parser.add_argument('-k', '--keggpath', type=str, dest='keggpath', 
        required=True, help='path to kegg ligand database')
    parser.add_argument('-r', '--reactiongraphsdir', type=str, dest='rgrapgsdir', 
        required=True, help='file with listing of used kegg reactions')
    parser.add_argument('-o', '--output', type=str, dest='output', 
        required=True, help='path to store output result')
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
                logger.warn('Output file ' + outputfile + ' exists and will be overwritten'
                    ' since parameter force was set.')
            except OSError:
                logger.error("Output file " + outputfile + "exists already"
                + "but can't remove. Aborting.");
                sys.exit(1)
        else:
            logger.error('Output file ' + outputfile + ' exists already, aborting.');
            sys.exit(1)

# main function
if __name__ == '__main__':

    logging.basicConfig(level=logging.INFO, format="%(levelname)s - %(message)s")
    logger = logging.getLogger()

    args = parseargs()

    # convert paths
    keggpath = os.path.abspath(args.keggpath)
    rgrapgsdir = os.path.abspath(args.rgrapgsdir)

    # remove ending slash on path if existing
    if keggpath.endswith('/'):
        keggpath = keggpath[:-1]

    # show more output
    if (args.verbose):
        logger.setLevel(logging.DEBUG)

    # prepare output
    check_output(os.path.abspath(args.output), args.force)

    # dictionary to store output and set for unique ec-codes
    reactiondir = {}
    eccodes = set()

    # read all reactiongraphs 
    reactionlist_content = " ".join(line.strip() for line in os.listdir(rgrapgsdir))

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
                    new_ecnumbers = re.findall(r"((?:(?:\d+|-)\.){3}(?:\d+|-))", line)
                    ecnumbers.extend(new_ecnumbers)
                except AttributeError:
                    pass

            # if an entry was read (meaning reaction number and ec codes)
            if line == "///\n" and reaction_read and ecnumbers_read:
                
                reaction_read = False
                ecnumbers_read = False
                # and if the entry exists in the reactionslist file
                # then store it
                if re.search(reaction, reactionlist_content):
                    reactions_found += 1
                    reactiondir[reaction] = ecnumbers
                    # add ec numbers to set of unique ec codes
                    for ec in ecnumbers: eccodes.add(ec)

    # convert unique ec-codes to sorted list
    sorted_ecs = sorted(list(eccodes))

    logger.info("found " + str(reactions_found)
        + " reactions out of " + str(reactions_total))
    logger.info("extracted " + str(len(sorted_ecs))
        + " unique ec-codes from found reactions")

    # write output file
    with open(args.output, 'w') as outfile:
        for reaction in reactiondir:
            output_line = ""
            logger.debug("found reaction " + reaction + " with ec-codes:")
            
            for ec in sorted_ecs: 
                if ec in reactiondir[reaction]:
                    output_line += "1 "
                    logger.debug(ec)
                else:
                    output_line += "0 "
            output_line = output_line[:-1] + "\n"
            outfile.write(output_line)
