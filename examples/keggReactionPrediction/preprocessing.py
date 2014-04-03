# Runs all preprocessing steps for preprocessing mol files from KEGG for 
# reaction prediction task
#
# @author Clemens Westrup

import os, errno, shutil, sys, argparse, time, subprocess

# parse command line arguments
def parseargs():
	parser = argparse.ArgumentParser(description='Preprocess kegg files'
		+ ' for use with path-index-graph-kernel algorithm.')
	parser.add_argument('-k', '--keggpath', type=str, dest='keggpath', 
		required=True, help='path to kegg ligand database')
	parser.add_argument('-v', '--verbose', dest='verbose', action='store_true',
		help='show verbose information')
	parser.add_argument('-o', '--output', type=str, dest='output', 
		required=True, help='path to store output results')
	parser.add_argument('-f', '--force', dest='force', action='store_true',
		help='force overwrite if output dir exists already')

	args = parser.parse_args()
	return args

# mkdir -p
def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc:
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else: raise

# main function
if __name__ == '__main__':
	args = parseargs()

	# read arguments
	verbose = args.verbose
	force = args.force
	keggpath = os.path.abspath(args.keggpath)
	outputdir = os.path.abspath(args.output)

	molfeaturesdir = outputdir + "/mol-features"
	reactiongraphsdir = outputdir + "/reaction-graphs"

	print ""
	print "--- preprocessing kegg data for path index graph kernel ---"
	print ""

	time.sleep(1) # delay for stout output

	# remove ending slash on path if existing
	if keggpath.endswith('/'):
		keggpath = keggpath[:-1]

	# show or hide subprocess output
	if (verbose):
		OUTPUT = None
	else:
		OUTPUT = open(os.devnull, 'w')

	# check if output dir exists
	if (os.path.isdir(outputdir)):
		if (force):
			print ('Directory "' + outputdir 
				+ '/ for storing results exists already and will be '
				+ 'overwritten since force overwrite was set.')
			shutil.rmtree(outputdir)
		else: 
			print ('Directory "' + outputdir 
				+ '" for storing results exists already, aborting.')
			sys.exit(1)

	# create temp dir
	print "Writing results in dir: " + outputdir
	mkdir_p(outputdir)

	### reaction list

	# Extract reaction list
	print "Extracting reaction list from kegg, this might take a while..."
	p = subprocess.Popen(
		'python reaction-list/extract-reactions.py -k ' + keggpath 
		+ ' > ' + outputdir + '/kegg-reactions.txt',
		shell=True, stdout=OUTPUT)
	p.wait()
	if p.returncode != 0:
		sys.exit(1)

	### mol-features

	# create output dir for mol-features
	print "Creating dir " + molfeaturesdir
	mkdir_p(molfeaturesdir)
	# create atom features from KEGG mol files
	print "Creating atom features from KEGG mol files ..."
	p = subprocess.Popen(
		'python feature-generator/generator2011.py ' + keggpath + '/mol/*'
		+ ' -k all --output-dir ' + molfeaturesdir + '/',
		shell=True, stdout=OUTPUT)
	p.wait()
	if p.returncode != 0:
		sys.exit(1)

	### reaction graphs

	# create output dir reaction graphs
	print "Creating temp dir " + reactiongraphsdir
	mkdir_p(reactiongraphsdir)
	# Create reactiongraph mol files with atommapper
	print "Create reactiongraph mol files with atommapper ..."
	os.chdir("./atommapper/bin");
	p = subprocess.Popen(
		'java Mapper2000 -rgraphs -moldir ' + keggpath + '/mol/ '
		+ '-featdir ' + molfeaturesdir + '/ '
		+ '-reacfile ' + outputdir + '/kegg-reactions.txt ' 
		+ '-output ' + reactiongraphsdir + '/',
		shell=True, stdout=OUTPUT)
	p.wait()
	if p.returncode != 0:
		sys.exit(1)

	### extract label matrix

	# Create reactiongraph mol files with atommapper
	print "Extracting label matrix into file \"labels\""
	os.chdir("../../extractlabels");
	p = subprocess.Popen(
		'python extractlabels.py -k ' + keggpath + ' -r ' + reactiongraphsdir 
		+ ' -o ' + outputdir + '/labels -f', shell=True, stdout=OUTPUT)
	p.wait()
	if p.returncode != 0:
		sys.exit(1)

	print("All done. Reaction graph .mol files can be "
		+ "found in " + os.path.abspath(reactiongraphsdir))
